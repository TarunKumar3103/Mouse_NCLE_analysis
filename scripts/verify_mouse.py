#!/usr/bin/env python3
"""Independent reconstruction from original sources and cached public responses.

No functions from preparation or analysis scripts are imported.
Fits logistic-score equations with scipy.optimize.root rather than the
analysis code's trust-exact likelihood minimization.
"""
import argparse,csv,io,json,re,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import optimize,special,stats

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--results-dir',type=Path);p.add_argument('--verification-dir',type=Path)
    a=p.parse_args();root=a.root
    raw=root/'raw_sources';results=a.results_dir or root/'results';out=a.verification_dir or root/'verification';out.mkdir(parents=True,exist_ok=True)
    native={}
    with zipfile.ZipFile(raw/'Mouse_Native_HQ_GE.zip') as z:
        for name in z.namelist():
            if name.endswith('/'):continue
            m=re.fullmatch(r'AF-(.+)-F1-model_v6\.csv',Path(name).name)
            if not m:raise AssertionError('Unexpected native model')
            accession=m.group(1);rows=list(csv.DictReader(io.StringIO(z.read(name).decode()),delimiter='|'))
            covalent=any(r['CCbond']=='True' for r in rows)
            detected=any(r['CCbond']=='False' and r['ENT']=='True' for r in rows)
            high=any(r['CCbond']=='False' and r['ENT']=='True' and r['Quality']=='High' for r in rows)
            native[accession]={'NCLE_any':None if covalent else int(detected),
                'NCLE_High_only':None if covalent or (detected and not high) else int(high),
                'Native_covalent':int(covalent),'Model_ID':Path(name).stem}
    with (raw/'knot_all.txt').open() as f:
        reader=csv.DictReader(f,delimiter=';');knots={r['Uniprot ID'].strip() for r in reader}
    lassos=set()
    with zipfile.ZipFile(raw/'AlphaLasso.zip') as z:
        for name in z.namelist():
            if not name.endswith('.tsv') or '__MACOSX' in name:continue
            with z.open(name) as f:
                for chunk in pd.read_csv(f,sep='\t',dtype=str,chunksize=100000):
                    mouse=chunk.loc[chunk.Organism.str.strip().eq('Mouse'),'Uniprot']
                    lassos.update(re.sub(r'-F\d+$','',str(x).strip()) for x in mouse)
    keyword={}
    with (raw/'uniprot_mouse_keyword_snapshot_2026_10_10.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['Organism (ID)']!='10090':raise AssertionError('Wrong keyword organism')
            kw={s.strip() for s in r['Keyword ID'].split(';')}
            keyword[r['Entry']]={'Membrane':'KW-0472' in kw,'Secreted':'KW-0964' in kw}
    source=pd.read_excel(raw/'McShane_mouse_mmc1.xlsx')
    source=source.rename(columns={'Protein IDs (Uniprot)':'UniProt_ID','Degradation profile':'Degradation_profile'})
    rebuilt=[]
    for r in source.itertuples(index=False):
        acc=r.UniProt_ID.strip();n=native.get(acc)
        if not n:continue
        mp=raw/'raw_prediction_json'/f'{acc}.json'
        if not mp.exists():continue
        data=json.loads(mp.read_text());candidates=[x for x in data if x.get('entryId',x.get('modelEntityId'))==f'AF-{acc}-F1']
        if len(candidates)!=1:continue
        m=candidates[0]
        if m['latestVersion']!=6 or m.get('taxId')!=10090 or not m['pdbUrl'].endswith(n['Model_ID']+'.pdb'):continue
        seq=m.get('sequence',m.get('uniprotSequence'))
        if m.get('sequenceStart',m.get('uniprotStart'))!=1 or m.get('sequenceEnd',m.get('uniprotEnd'))!=len(seq):continue
        cp=raw/'raw_confidence_json'/f'{acc}_confidence_v6.json'
        if not cp.exists():continue
        confidence=json.loads(cp.read_text());values=confidence['confidenceScore']
        if confidence['residueNumber']!=list(range(1,len(values)+1)) or len(values)!=len(seq):continue
        if not values or not np.isfinite(values).all() or min(values)<0 or max(values)>100:continue
        mean=float(np.mean(values))
        if mean<70 or n['NCLE_any'] is None or acc in knots or n['Native_covalent'] or acc in lassos:continue
        kw=keyword.get(acc)
        if kw is None or kw['Membrane'] or kw['Secreted']:continue
        rebuilt.append({'UniProt_ID':acc,'Degradation_profile':r.Degradation_profile,'Protein_length':len(seq),'Mean_pLDDT':mean,
            'NCLE_any':n['NCLE_any'],'NCLE_High_only':n['NCLE_High_only']})
    rebuilt=pd.DataFrame(rebuilt).sort_values('UniProt_ID')
    reference=pd.read_csv(results/'association_results.csv');checks=[];memberships=[]
    for analysis,tag,labels in [('Primary: NED vs ED','primary',{'ED','NED'}),('Sensitivity: NED vs ED + UN','sensitivity',{'ED','NED','UN'})]:
        base_cohort=rebuilt.loc[rebuilt.Degradation_profile.isin(labels)].copy()
        for definition,col,filetag in [('Any detected NCLE','NCLE_any','any_ncle'),('High-quality detection only','NCLE_High_only','high_only')]:
            cohort=base_cohort.loc[base_cohort[col].notna()].copy()
            shipped=pd.read_csv(results/f'{tag}_{filetag}_cohort.csv').sort_values('UniProt_ID')
            if set(cohort.UniProt_ID)!=set(shipped.UniProt_ID):raise AssertionError(f'Cohort IDs differ: {analysis}')
            merged=cohort.merge(shipped,on='UniProt_ID',validate='one_to_one',suffixes=('_independent','_shipped'))
            for k in ['Degradation_profile','Protein_length','NCLE_any','NCLE_High_only']:
                if not merged[k+'_independent'].fillna(-1).eq(merged[k+'_shipped'].fillna(-1)).all():raise AssertionError(f'Field differs: {k}')
            if not np.allclose(merged.Mean_pLDDT_independent,merged.Mean_pLDDT_shipped,atol=1e-10,rtol=0):raise AssertionError('Mean pLDDT differs')
            memberships.append({'Analysis':analysis,'NCLE_definition':definition,'N':len(cohort),'Exact_ID_match':True,'Field_disagreements':0})
            nc=cohort[col].to_numpy(dtype=float);y=cohort.Degradation_profile.eq('NED').to_numpy(dtype=float)
            length=cohort.Protein_length.to_numpy(dtype=float);z=(length-np.mean(length))/np.std(length,ddof=1)
            x=np.column_stack([np.ones(len(y)),nc,z])
            def score(beta):return x.T@(special.expit(x@beta)-y)
            def info(beta):
                prob=special.expit(x@beta);return x.T@((prob*(1-prob))[:,None]*x)
            fit=optimize.root(score,[special.logit(y.mean()),0,0],jac=info,tol=1e-11)
            if np.max(np.abs(score(fit.x)))>1e-7:raise AssertionError('Independent score fit failed')
            errors=np.sqrt(np.diag(np.linalg.inv(info(fit.x))))
            target=reference.loc[reference.Analysis.eq(analysis)&reference.NCLE_definition.eq(definition)]
            for term,beta,se in zip(['Intercept','NCLE','z_Protein_length'],fit.x,errors):
                ref=target.loc[target.Method.str.startswith('Model 1')&target.Term.eq(term)].iloc[0]
                values={'Coefficient':beta,'SE':se,'OR':np.exp(beta),'CI_lower':np.exp(beta-stats.norm.ppf(.975)*se),
                    'CI_upper':np.exp(beta+stats.norm.ppf(.975)*se),'P_value':2*stats.norm.sf(abs(beta/se))}
                for field,value in values.items():
                    match=bool(np.isclose(value,ref[field],rtol=1e-6,atol=1e-8))
                    checks.append({'Analysis':analysis,'NCLE_definition':definition,'Method':'Model 1','Term':term,'Field':field,
                        'Independent':float(value),'Shipped':float(ref[field]),'Pass':match})
            a1=int(((nc==1)&(y==1)).sum());b=int(((nc==1)&(y==0)).sum());c=int(((nc==0)&(y==1)).sum());d=int(((nc==0)&(y==0)).sum())
            odds=a1*d/(b*c);se=np.sqrt(1/a1+1/b+1/c+1/d)
            # Enumerate the fixed-margin hypergeometric distribution for a two-sided Fisher test.
            total=a1+b+c+d;events=a1+c;exposed=a1+b
            support=np.arange(max(0,exposed-(total-events)),min(exposed,events)+1)
            probs=stats.hypergeom.pmf(support,total,events,exposed);obs=stats.hypergeom.pmf(a1,total,events,exposed)
            fisher=float(probs[probs<=obs*(1+1e-12)].sum())
            values={'OR':odds,'CI_lower':np.exp(np.log(odds)-stats.norm.ppf(.975)*se),'CI_upper':np.exp(np.log(odds)+stats.norm.ppf(.975)*se),'P_value':fisher}
            ref=target.loc[target.Method.eq('Crude Fisher exact')].iloc[0]
            for field,value in values.items():
                checks.append({'Analysis':analysis,'NCLE_definition':definition,'Method':'Crude','Term':'NCLE','Field':field,
                    'Independent':float(value),'Shipped':float(ref[field]),'Pass':bool(np.isclose(value,ref[field],rtol=1e-6,atol=1e-8))})
    checks=pd.DataFrame(checks);checks.to_csv(out/'independent_numerical_checks.csv',index=False)
    pd.DataFrame(memberships).to_csv(out/'independent_cohort_checks.csv',index=False)
    summary={'Independent_cohort_IDs_and_fields':'PASS','Numerical_fields_checked':len(checks),'Numerical_failures':int((~checks.Pass).sum()),
        'Independent_regression_method':'Solve logistic score equations; no analysis-script functions imported',
        'Independent_crude_P_method':'Enumerate fixed-margin hypergeometric probabilities'}
    (out/'verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not checks.Pass.all():raise AssertionError('Independent numerical check failed')

if __name__=='__main__':main()
