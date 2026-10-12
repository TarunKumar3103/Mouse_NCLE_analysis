#!/usr/bin/env python3
"""Mouse NCLE association with McShane NED, crude and length-adjusted Model 1.

Main NCLE rule matches the human analysis: any native noncovalent detection.
The additional local-quality check requires a detection with Quality=High.
Mean protein pLDDT>=70 is required in every analysis. Disorder is not a filter.
"""
import argparse, hashlib, json, platform, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from stats_model1 import crude, fit_model1

def fingerprint(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return {'filename':p.name,'bytes':p.stat().st_size,'sha256':h.hexdigest()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--prepared',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--confidence',type=Path,required=True)
    p.add_argument('--knots',type=Path,required=True);p.add_argument('--alphalasso',type=Path,required=True)
    p.add_argument('--keywords',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):raise ValueError('Use a new/empty result directory')
    a.out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.prepared,dtype={'UniProt_ID':str})
    metadata=pd.read_csv(a.metadata,dtype={'UniProt_ID':str})
    if not d.UniProt_ID.is_unique or not metadata.UniProt_ID.is_unique:raise ValueError('Duplicate merge IDs')
    if metadata.Metadata_status.str.startswith('request_failed').any():raise ValueError('Retry transient metadata-query failures before constructing cohorts')
    expected=set(d.loc[d.Model_ID.notna(),'UniProt_ID'])
    if set(metadata.UniProt_ID)!=expected:raise ValueError('Metadata query must cover every exact native match')
    metadata=metadata.rename(columns={'Mean_pLDDT':'AFDB_reported_mean_pLDDT','pLDDT_method':'AFDB_reported_pLDDT_method'})
    d=d.merge(metadata.drop(columns=['Sequence'],errors='ignore'),on='UniProt_ID',how='left',validate='one_to_one')
    confidence=pd.read_csv(a.confidence,dtype={'UniProt_ID':str})
    if not confidence.UniProt_ID.is_unique or not set(confidence.UniProt_ID).issubset(expected):raise ValueError('Invalid residue-confidence accession coverage')
    d=d.merge(confidence.rename(columns={'Residue_mean_pLDDT':'Mean_pLDDT'}),on='UniProt_ID',how='left',validate='one_to_one')
    valid=d.Metadata_status.eq('ok')&d.Confidence_status.eq('ok')&d.Protein_length.eq(d.Confidence_length)
    if not d.loc[valid,'AF_model_ID'].eq(d.loc[valid,'Model_ID']).all():raise ValueError('Metadata model mismatch')
    if not d.loc[valid,'AF_model_version'].eq(6).all():raise ValueError('Wrong AlphaFold version')
    if not d.loc[valid,'Confidence_model_ID'].eq(d.loc[valid,'Model_ID']).all():raise ValueError('Confidence model mismatch')
    d['pLDDT_summary_difference']=d.Mean_pLDDT-d.AFDB_reported_mean_pLDDT
    d['pLDDT_cutoff_call_difference']=(d.Mean_pLDDT.notna()&d.AFDB_reported_mean_pLDDT.notna()&d.Mean_pLDDT.ge(70).ne(d.AFDB_reported_mean_pLDDT.ge(70))).astype(int)
    k=pd.read_csv(a.knots,sep=';',usecols=[0],dtype=str)
    knots=set(k.iloc[:,0].dropna().str.strip())
    d['Knotted']=d.UniProt_ID.isin(knots).astype(int)
    # Exact source-organism rule was confirmed by inspecting matched records.
    evidence=[];mouse_ids=set()
    with zipfile.ZipFile(a.alphalasso) as z:
        for name in z.namelist():
            if not name.endswith('.tsv') or '__MACOSX' in name:continue
            with z.open(name) as f:
                for c in pd.read_csv(f,sep='\t',dtype=str,chunksize=100000):
                    c=c.loc[c.Organism.str.strip().eq('Mouse')].copy()
                    c['UniProt_ID']=c.Uniprot.str.strip().str.replace(r'-F\d+$','',regex=True)
                    mouse_ids.update(c.UniProt_ID)
                    e=c.loc[c.UniProt_ID.isin(d.UniProt_ID)].copy()
                    if len(e):e['AlphaLasso_table']=name;evidence.append(e)
    evidence=pd.concat(evidence,ignore_index=True) if evidence else pd.DataFrame(columns=['UniProt_ID'])
    d['AlphaLasso_covalent']=d.UniProt_ID.isin(mouse_ids).astype(int)
    d['Covalent_any_source']=(d.Covalent_native_records.eq(1)|d.AlphaLasso_covalent.eq(1)).astype(int)
    kw=pd.read_csv(a.keywords,sep='\t',dtype=str)
    if not kw.Entry.is_unique or not kw['Organism (ID)'].eq('10090').all():raise ValueError('Invalid keyword snapshot')
    kw=kw.rename(columns={'Entry':'UniProt_ID','Length':'Current_UniProt_length','Organism (ID)':'Current_taxonomy_ID',
        'Organism':'Current_organism','Keyword ID':'UniProt_keyword_IDs','Keywords':'UniProt_keywords'})
    d=d.merge(kw,on='UniProt_ID',how='left',validate='one_to_one')
    d['Keyword_source_covered']=d.Current_taxonomy_ID.notna()
    sets=d.UniProt_keyword_IDs.fillna('').str.split(';').map(lambda x:{s.strip() for s in x})
    for column,keyword in [('Membrane_keyword','KW-0472'),('Secreted_keyword','KW-0964')]:
        d[column]=pd.array([int(keyword in s) if covered else pd.NA for s,covered in zip(sets,d.Keyword_source_covered)],dtype='Int64')
    d['Current_UniProt_length']=pd.to_numeric(d.Current_UniProt_length,errors='coerce')
    d['AF_UniProt_length_disagreement']=(d.Protein_length.notna()&d.Current_UniProt_length.notna()&d.Protein_length.ne(d.Current_UniProt_length)).astype(int)
    prequality=d.Model_ID.notna()&d.NCLE_any.notna()&d.Knotted.eq(0)&d.Covalent_any_source.eq(0)&d.Keyword_source_covered&d.Membrane_keyword.eq(0)&d.Secreted_keyword.eq(0)
    if set(confidence.UniProt_ID)!=set(d.loc[prequality,'UniProt_ID']):raise ValueError('Residue-confidence queries must cover all and only pre-quality candidates')
    if confidence.Confidence_status.str.startswith('request_failed').any():raise ValueError('Retry transient confidence-query failures; do not exclude proteins based on a recoverable request failure')
    filters=[('Require exact native accession match',d.Model_ID.notna()),
        ('Require binary native NCLE status',d.NCLE_any.notna()),('Exclude supplied knotted proteins',d.Knotted.eq(0)),
        ('Exclude supplied or detected covalent lassos',d.Covalent_any_source.eq(0)),
        ('Require mouse UniProt keyword coverage',d.Keyword_source_covered),
        ('Exclude UniProt Membrane keyword',d.Membrane_keyword.eq(0)),('Exclude UniProt Secreted keyword',d.Secreted_keyword.eq(0)),
        ('Require matching v6 metadata and residue confidence',valid),('Retain protein mean pLDDT >=70',d.Mean_pLDDT.ge(70))]
    flows=[];result_rows=[];counts=[];diagnostics=[];cohort_sizes={}
    d['Core_eligible']=False;d['First_exclusion']=''
    remaining=pd.Series(True,index=d.index)
    for step,mask in filters:
        mask=mask.fillna(False)
        rejected=remaining&~mask;d.loc[rejected,'First_exclusion']=step;remaining&=mask
    d['Core_eligible']=remaining
    for analysis,profiles in [('Primary: NED vs ED',{'ED','NED'}),('Sensitivity: NED vs ED + UN',{'ED','NED','UN'})]:
        selected=d.Degradation_profile.isin(profiles)
        def flow(step,old=None):
            subset=d.loc[selected];pcounts=subset.Degradation_profile.value_counts()
            flows.append({'Analysis':analysis,'Step':step,'Total':len(subset),'Removed':None if old is None else old-len(subset),
                **{k:int(pcounts.get(k,0)) for k in ['ED','NED','UN']}})
        flow('Select experimental degradation categories')
        for step,mask in filters:
            old=int(selected.sum());selected&=mask.fillna(False);flow(step,old)
        cohort=d.loc[selected].copy()
        if cohort.empty:raise ValueError('Empty final cohort')
        if not cohort.UniProt_ID.is_unique:raise ValueError('Duplicate cohort IDs')
        cohort['NED_binary']=cohort.Degradation_profile.eq('NED').astype(int)
        cohort_sizes[analysis]=len(cohort)
        for definition,col in [('Any detected NCLE','NCLE_any'),('High-quality detection only','NCLE_High_only')]:
            # Low-only detections remain positive in the main human-consistent
            # rule, but are uncertain (excluded) in the conservative High-only check.
            q=cohort.loc[cohort[col].notna()].copy();q['NCLE']=q[col].astype(int)
            estimate,ct=crude(q)
            result_rows.append({'Analysis':analysis,'NCLE_definition':definition,'Method':'Crude Fisher exact','Term':'NCLE',**estimate})
            for r in ct:counts.append({'Analysis':analysis,'NCLE_definition':definition,**r,'NED_frequency':r['NED']/r['Total']})
            rows,diag=fit_model1(q)
            if diag['Status']!='ok':raise ValueError(f'Model1 diagnostic failed: {diag}')
            result_rows.extend({'Analysis':analysis,'NCLE_definition':definition,'Method':'Model 1: NED ~ NCLE + z(length)',**r} for r in rows)
            diagnostics.append({'Analysis':analysis,'NCLE_definition':definition,**diag})
            fname=('primary' if analysis.startswith('Primary') else 'sensitivity')+('_any_ncle' if col=='NCLE_any' else '_high_only')+'_cohort.csv'
            q.sort_values('UniProt_ID').to_csv(a.out/fname,index=False)
    d.sort_values('UniProt_ID').to_csv(a.out/'all_experimental_record_audit.csv',index=False)
    pd.DataFrame(flows).to_csv(a.out/'cohort_flow.csv',index=False)
    res=pd.DataFrame(result_rows);res.to_csv(a.out/'association_results.csv',index=False)
    pd.DataFrame(counts).to_csv(a.out/'contingency_counts.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(a.out/'model_diagnostics.csv',index=False)
    evidence.to_csv(a.out/'covalent_lasso_evidence.csv',index=False)
    d.loc[d.NCLE_any.ne(d.NCLE_High_only)&d.NCLE_any.notna()].to_csv(a.out/'local_quality_definition_differences.csv',index=False)
    d.loc[d.AF_UniProt_length_disagreement.eq(1)].to_csv(a.out/'sequence_length_version_differences.csv',index=False)
    d.loc[d.Mean_pLDDT.notna(),['UniProt_ID','AFDB_reported_mean_pLDDT','Mean_pLDDT','pLDDT_summary_difference','pLDDT_cutoff_call_difference']].to_csv(a.out/'pLDDT_summary_comparison.csv',index=False)
    summary={'cohort_sizes':cohort_sizes,'mean_pLDDT_cutoff':'>=70','protein_length_source':'Matching AlphaFold model v6 sequence length',
        'main_NCLE':'Any ENT=True and CCbond=False row; any native covalent row excluded',
        'local_quality_check':'Require a Quality=High noncovalent detection for positives; exclude Low-only detections; empty files remain no detected NCLE',
        'pLDDT_method':'Arithmetic mean of exact-v6 confidenceScore values, one per residue; 1..N coordinates and sequence length verified',
        'mouse_AlphaLasso_IDs':len(mouse_ids),'metadata_statuses':metadata.Metadata_status.value_counts().to_dict(),
        'eligible_AF_UniProt_length_disagreements':int(d.loc[d.Core_eligible,'AF_UniProt_length_disagreement'].sum()),
        'eligible_NCLE_definition_differences':int((d.Core_eligible&d.NCLE_any.ne(d.NCLE_High_only)).sum()),
        'confidence_statuses':confidence.Confidence_status.value_counts().to_dict(),
        'AFDB_vs_residue_mean_cutoff_disagreements':int(d.pLDDT_cutoff_call_difference.sum()),
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
        'sources':[fingerprint(x) for x in [a.prepared,a.metadata,a.confidence,a.knots,a.alphalasso,a.keywords]],
        'notes':['MobiDB disorder and complex membership are not used.','UN is included only in the sensitivity analysis.',
            'Exact accession matching preserves isoform suffixes.','Empty native detection CSVs are treated as no detected NCLE; completion logs were not supplied.',
            'Any-source covalent and supplied knot exclusions apply regardless of detection quality or source model version.',
            'UniProt keywords are a current frozen snapshot, not historical McShane annotations.',
            'The High-only check excludes Low-only detections rather than classifying them as confirmed negatives.']}
    (a.out/'analysis_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print(res.loc[res.Term.eq('NCLE')].to_string(index=False))

if __name__=='__main__':main()
