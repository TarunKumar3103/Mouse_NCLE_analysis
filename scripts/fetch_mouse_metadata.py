#!/usr/bin/env python3
"""Freeze official AlphaFold metadata, requiring an exact v6/F1 model match."""
import argparse, concurrent.futures, datetime, hashlib, json, time, urllib.error, urllib.request
from pathlib import Path
import pandas as pd

def fetch(row,cache):
    accession=row.UniProt_ID; expected=row.Model_ID
    output=cache/f'{accession}.json'
    url=f'https://alphafold.ebi.ac.uk/api/prediction/{accession}'
    if output.exists():
        body=output.read_bytes()
    else:
        last=None
        for attempt in range(3):
            try:
                req=urllib.request.Request(url,headers={'User-Agent':'NCLE-reproducibility-analysis/1.0'})
                with urllib.request.urlopen(req,timeout=25) as r:body=r.read()
                json.loads(body)
                output.write_bytes(body)
                break
            except urllib.error.HTTPError as e:
                if e.code==404:return {'UniProt_ID':accession,'Metadata_status':'http_404','Metadata_URL':url}
                last=f'HTTP_{e.code}'
                time.sleep(2*(attempt+1))
            except Exception as e:
                last=type(e).__name__;time.sleep(attempt+1)
        else:return {'UniProt_ID':accession,'Metadata_status':f'request_failed:{last}','Metadata_URL':url}
    data=json.loads(body)
    candidates=[r for r in data if r.get('modelEntityId',r.get('entryId'))==expected.replace('-model_v6','')]
    if len(candidates)!=1:return {'UniProt_ID':accession,'Metadata_status':'ambiguous_or_missing_exact_model','Metadata_URL':url}
    r=candidates[0]
    if r.get('taxId')!=10090:
        return {'UniProt_ID':accession,'Metadata_status':'wrong_species','Metadata_URL':url}
    if r.get('latestVersion')!=6 or not r.get('pdbUrl','').endswith(expected+'.pdb'):
        return {'UniProt_ID':accession,'Metadata_status':'model_version_mismatch','Metadata_URL':url}
    sequence=r.get('sequence',r.get('uniprotSequence',''))
    start=r.get('sequenceStart',r.get('uniprotStart'));end=r.get('sequenceEnd',r.get('uniprotEnd'))
    if start!=1 or end!=len(sequence) or not sequence:
        return {'UniProt_ID':accession,'Metadata_status':'partial_or_inconsistent_sequence','Metadata_URL':url}
    mean=float(r['globalMetricValue']);method='AFDB globalMetricValue (retained for comparison; not the final cutoff field)'
    if not 0<=mean<=100:raise ValueError('Invalid mean pLDDT')
    return {'UniProt_ID':accession,'Metadata_status':'ok','Protein_length':len(sequence),'Mean_pLDDT':mean,
        'pLDDT_method':method,'AF_model_version':6,'AF_model_ID':expected,'AF_model_created_date':r.get('modelCreatedDate'),
        'AF_sequence_SHA256':hashlib.sha256(sequence.encode()).hexdigest(),'Metadata_URL':url,
        'Metadata_raw_SHA256':hashlib.sha256(body).hexdigest(),'Sequence':sequence}

def main():
    p=argparse.ArgumentParser();p.add_argument('--ids',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    cache=a.out/'raw_prediction_json';cache.mkdir(exist_ok=True)
    ids=pd.read_csv(a.ids,dtype={'UniProt_ID':str});results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs=[pool.submit(fetch,r,cache) for r in ids.itertuples(index=False)]
        for job in concurrent.futures.as_completed(jobs):
            results.append(job.result())
            if len(results)%50==0 or len(results)==len(ids):
                pd.DataFrame(results).sort_values('UniProt_ID').to_csv(a.out/'alphafold_metadata.csv',index=False)
                print(f'AlphaFold metadata: {len(results)}/{len(ids)}; {sum(r["Metadata_status"]=="ok" for r in results)} matched v6',flush=True)
    d=pd.DataFrame(results).sort_values('UniProt_ID');d.to_csv(a.out/'alphafold_metadata.csv',index=False)
    manifest={'retrieved_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':'Official AlphaFold DB API',
        'required_model_version':6,'queries':len(ids),'statuses':d.Metadata_status.value_counts().to_dict(),
        'threshold':'Final cutoff uses separately retrieved per-residue means, not these summary values'}
    (a.out/'metadata_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':main()
