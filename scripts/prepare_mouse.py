#!/usr/bin/env python3
"""Audit supplied mouse native records. Preserve all-record and High-only calls."""
import argparse, csv, hashlib, io, json, re, zipfile
from collections import Counter
from pathlib import Path
import pandas as pd

def fingerprint(p):
    return {"name":p.name,"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--mcshane',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    annotations=[];quality=Counter();reasons=Counter();columns=None
    with zipfile.ZipFile(a.native) as z:
        bad=z.testzip()
        if bad:raise ValueError(f'Bad ZIP CRC: {bad}')
        for name in z.namelist():
            if name.endswith('/'):continue
            m=re.fullmatch(r'AF-(.+)-F(\d+)-model_v(\d+)\.csv',Path(name).name)
            if not m:raise ValueError(f'Unexpected file: {name}')
            accession,fragment,version=m.groups()
            reader=csv.DictReader(io.StringIO(z.read(name).decode()),delimiter='|')
            if columns is None:columns=reader.fieldnames
            if reader.fieldnames!=columns:raise ValueError(f'Unexpected schema: {name}')
            rows=list(reader)
            for r in rows:
                if r['ID']!=Path(name).stem:raise ValueError(f'ID disagreement: {name}')
                if r['CCbond'] not in ['True','False'] or r['ENT']!='True' or r['Quality'] not in ['High','Low']:
                    raise ValueError(f'Unexpected annotation: {name}: {r}')
                quality[r['Quality']]+=1;reasons[r['Reason']]+=1
            covalent=any(r['CCbond']=='True' for r in rows)
            any_ncle=any(r['CCbond']=='False' and r['ENT']=='True' for r in rows)
            high_ncle=any(r['CCbond']=='False' and r['ENT']=='True' and r['Quality']=='High' for r in rows)
            annotations.append({'UniProt_ID':accession,'Model_ID':Path(name).stem,
                'Archive_file':name,'Fragment':int(fragment),'Model_version':int(version),
                'Detected_rows':len(rows),'High_rows':sum(r['Quality']=='High' for r in rows),
                'Low_rows':sum(r['Quality']=='Low' for r in rows),
                'Covalent_native_records':int(covalent),
                'NCLE_any':None if covalent else int(any_ncle),
                'NCLE_High_only':None if covalent or (any_ncle and not high_ncle) else int(high_ncle),
                'Low_only_noncovalent_detection':int(any_ncle and not high_ncle and not covalent),
                'Archive_status':'covalent_exclusion' if covalent else 'detected_noncovalent' if any_ncle else 'empty_detection_file'})
    n=pd.DataFrame(annotations).sort_values('UniProt_ID')
    if not n.UniProt_ID.is_unique or not n.Fragment.eq(1).all() or not n.Model_version.eq(6).all():
        raise ValueError('Multiple/non-F1 or non-v6 models need an explicit aggregation policy')
    d=pd.read_excel(a.mcshane,sheet_name='FinalTable_PostNorm_3T3_')
    d=d.rename(columns={'Protein IDs (Uniprot)':'UniProt_ID','Degradation profile':'Degradation_profile'})
    d.UniProt_ID=d.UniProt_ID.astype('string').str.strip()
    d.Degradation_profile=d.Degradation_profile.astype('string').str.strip()
    d['Source_Excel_row']=range(2,len(d)+2)
    if d.UniProt_ID.isna().any() or not d.UniProt_ID.is_unique or not d.Degradation_profile.isin(['ED','NED','UN']).all():
        raise ValueError('Invalid McShane IDs/profiles')
    d['Isoform_ID']=d.UniProt_ID.str.contains(r'-\d+$').astype(int)
    merged=d.merge(n,on='UniProt_ID',how='left',validate='one_to_one')
    if len(merged)!=len(d):raise ValueError('Merge expanded rows')
    n.to_csv(a.out/'native_protein_annotations.csv',index=False)
    merged.to_csv(a.out/'mcshane_native_exact_merge.csv',index=False)
    matched=merged.loc[merged.Model_ID.notna()]
    matched[['UniProt_ID','Model_ID','Model_version']].to_csv(a.out/'metadata_query_ids.csv',index=False)
    d.loc[~d.UniProt_ID.isin(n.UniProt_ID)].to_csv(a.out/'unmatched_degradation_records.csv',index=False)
    overlap=pd.crosstab(merged.Degradation_profile,merged.Archive_status.fillna('no_archive_file'))
    overlap.to_csv(a.out/'archive_profile_overlap.csv')
    summary={'native_files':len(n),'native_record_rows':int(n.Detected_rows.sum()),'empty_files':int(n.Detected_rows.eq(0).sum()),
        'record_quality_counts':dict(quality),'record_reason_counts':dict(reasons),'columns':columns,
        'native_covalent_proteins':int(n.Covalent_native_records.sum()),
        'any_ncle_proteins':int(n.NCLE_any.eq(1).sum()),'any_ncle_negative_proteins':int(n.NCLE_any.eq(0).sum()),
        'high_only_ncle_proteins':int(n.NCLE_High_only.eq(1).sum()),'high_only_negative_proteins':int(n.NCLE_High_only.eq(0).sum()),
        'noncovalent_low_only_proteins':int(n.Low_only_noncovalent_detection.sum()),
        'mcshane_rows':len(d),'mcshane_profiles':d.Degradation_profile.value_counts().to_dict(),
        'mcshane_isoform_ids':int(d.Isoform_ID.sum()),'exact_native_matches':len(matched),
        'unmatched_records':int(merged.Model_ID.isna().sum()),'sources':[fingerprint(a.native),fingerprint(a.mcshane)]}
    (a.out/'preparation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
