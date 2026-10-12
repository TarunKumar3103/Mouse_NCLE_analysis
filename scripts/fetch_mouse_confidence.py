#!/usr/bin/env python3
"""Freeze exact-v6 residue confidence and calculate the protein-level arithmetic mean."""
import argparse,concurrent.futures,datetime,email.utils,hashlib,json,threading,time,urllib.error,urllib.request
from pathlib import Path
import pandas as pd

class RateGate:
    def __init__(self,interval):self.interval=interval;self.next_request=0.;self.lock=threading.Lock()
    def wait(self):
        with self.lock:
            now=time.monotonic();wait=max(0,self.next_request-now);self.next_request=max(now,self.next_request)+self.interval
        if wait:time.sleep(wait)
    def pause(self,delay):
        with self.lock:self.next_request=max(self.next_request,time.monotonic()+delay)

def fetch(row,cache,gate):
    acc=row.UniProt_ID;model=row.Model_ID
    url=f'https://alphafold.ebi.ac.uk/files/{model.replace("-model_v6", "")}-confidence_v6.json'
    cp=cache/f'{acc}_confidence_v6.json'
    if cp.exists():body=cp.read_bytes()
    else:
        last=None
        for attempt in range(5):
            try:
                gate.wait()
                with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'NCLE-reproducibility-analysis/1.0'}),timeout=25) as r:body=r.read()
                json.loads(body);cp.write_bytes(body);break
            except urllib.error.HTTPError as e:
                if e.code==404:return {'UniProt_ID':acc,'Confidence_status':'http_404','Confidence_URL':url}
                last=f'HTTP_{e.code}'
                if e.code==429:
                    retry=e.headers.get('Retry-After','30')
                    try:delay=float(retry)
                    except ValueError:delay=max(0,(email.utils.parsedate_to_datetime(retry)-datetime.datetime.now(datetime.timezone.utc)).total_seconds())
                    gate.pause(max(1,delay))
                else:time.sleep(attempt+1)
            except Exception as e:last=type(e).__name__;time.sleep(attempt+1)
        else:return {'UniProt_ID':acc,'Confidence_status':f'request_failed:{last}','Confidence_URL':url}
    d=json.loads(body);scores=d.get('confidenceScore');numbers=d.get('residueNumber')
    if not isinstance(scores,list) or not scores or numbers!=list(range(1,len(scores)+1)):
        return {'UniProt_ID':acc,'Confidence_status':'invalid_residue_coordinates','Confidence_URL':url}
    if not all(isinstance(s,(int,float)) and 0<=s<=100 for s in scores):raise ValueError('Invalid residue confidence')
    return {'UniProt_ID':acc,'Confidence_status':'ok','Residue_mean_pLDDT':sum(scores)/len(scores),
        'Confidence_length':len(scores),'Confidence_URL':url,'Confidence_raw_SHA256':hashlib.sha256(body).hexdigest(),
        'Confidence_model_ID':model}

def main():
    p=argparse.ArgumentParser();p.add_argument('--ids',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=6);p.add_argument('--request-interval',type=float,default=.7)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);gate=RateGate(a.request_interval)
    cache=a.out/'raw_confidence_json';cache.mkdir(exist_ok=True)
    ids=pd.read_csv(a.ids,dtype={'UniProt_ID':str});result=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs=[pool.submit(fetch,r,cache,gate) for r in ids.itertuples(index=False)]
        for job in concurrent.futures.as_completed(jobs):
            result.append(job.result())
            if len(result)%100==0 or len(result)==len(ids):
                pd.DataFrame(result).sort_values('UniProt_ID').to_csv(a.out/'residue_confidence_means.csv',index=False)
                print(f'Residue means: {len(result)}/{len(ids)}; {sum(r["Confidence_status"]=="ok" for r in result)} valid',flush=True)
    d=pd.DataFrame(result).sort_values('UniProt_ID');d.to_csv(a.out/'residue_confidence_means.csv',index=False)
    manifest={'retrieved_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'Source':'Exact version 6 AlphaFold residue confidence JSON',
        'statuses':d.Confidence_status.value_counts().to_dict(),'Mean':'Arithmetic mean of one score per residue','Coordinates':'Consecutive 1..N'}
    (a.out/'confidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
