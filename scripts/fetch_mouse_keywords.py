#!/usr/bin/env python3
"""Freeze a mouse UniProt keyword snapshot; absence of a matched entry stays missing."""
import argparse, datetime, hashlib, json, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    fields='accession,length,organism_id,organism_name,keywordid,keyword'
    url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode({'query':'(organism_id:10090)','format':'tsv','fields':fields})
    req=urllib.request.Request(url,headers={'User-Agent':'NCLE-reproducibility-analysis/1.0'})
    with urllib.request.urlopen(req,timeout=60) as response:
        body=response.read();headers=dict(response.headers)
    path=a.out/'uniprot_mouse_keyword_snapshot_2026_10_10.tsv';path.write_bytes(body)
    d=pd.read_csv(path,sep='\t',dtype=str)
    if d.Entry.isna().any() or not d.Entry.is_unique or not d['Organism (ID)'].eq('10090').all():raise ValueError('Invalid keyword snapshot')
    manifest={'URL':url,'retrieved_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'entry_count':len(d),'SHA256':hashlib.sha256(body).hexdigest(),
        'response_headers':{k:v for k,v in headers.items() if k.lower() in ['x-uniprot-release','x-uniprot-release-date','x-total-results','last-modified']},
        'Membrane_keyword':'KW-0472','Secreted_keyword':'KW-0964'}
    (a.out/'uniprot_keyword_metadata.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':main()
