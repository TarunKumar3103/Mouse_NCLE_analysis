import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';

const root=process.argv[2];
if(!root)throw new Error('Pass the absolute Mouse_NCLE_analysis directory.');
const p=JSON.parse(await fs.readFile(path.join(root,'results/workbook_payload.json'),'utf8'));
const wb=Workbook.create();
const res=wb.worksheets.add('Results'),flow=wb.worksheets.add('Cohort flow'),records=wb.worksheets.add('Protein records');
const navy='#243D59',ink='#172B40',pale='#EAF0F6';
for(const s of [res,flow,records])s.showGridLines=false;
res.tabColor=navy;flow.tabColor='#63809B';
function style(s,range){s.getRange(range).format.font={name:'Arial',size:11,color:ink};s.getRange(range).format.rowHeight=23;s.getRange(range).format.verticalAlignment='center';}
function header(s,range){const r=s.getRange(range);r.format.fill=navy;r.format.font={name:'Arial',size:11,bold:true,color:'#FFFFFF'};r.format.horizontalAlignment='center';r.format.wrapText=true;r.format.rowHeight=36;r.format.borders={insideVertical:{style:'thin',color:'#FFFFFF'}};}
function title(s,text,end='G'){s.getRange('A2').values=[[text]];s.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:ink};s.getRange('A2').format.rowHeight=29;s.getRange(`A3:${end}3`).format.borders={bottom:{style:'thin',color:'#A9B9CA'}};}
const keys=['UniProt_ID','Degradation_profile','NCLE_any','NCLE_High_only','Mean_pLDDT','Protein_length','Knotted','Covalent_any_source','Membrane_keyword','Secreted_keyword','Core_eligible','NED_binary','Primary_included','Sensitivity_included','First_exclusion','AF_UniProt_length_disagreement','Model_ID','Metadata_status','Source_Excel_row','Current_UniProt_length'];
const matrix=[keys,...p.audit.map(r=>keys.map(k=>r[k]??null))],last=matrix.length;
style(records,`A1:T${last}`);records.getRangeByIndexes(0,0,last,20).values=matrix;
records.tables.add(`A1:T${last}`,true,'MouseProteinRecords');header(records,'A1:T1');records.getRange('A1:T1').format.rowHeight=49;
records.freezePanes.freezeRows(1);records.freezePanes.freezeColumns(1);
for(let c=0;c<20;c++)records.getRangeByIndexes(0,c,last,1).format.columnWidth=c===0?19:c===14?67:c===16?38:c===17?30:17;
records.getRange(`E2:E${last}`).setNumberFormat('0.00');records.getRange(`F2:T${last}`).setNumberFormat('0');
for(const c of ['O','Q','R'])records.getRange(`${c}2:${c}${last}`).setNumberFormat('@');
records.getRange('L2').formulas=[['=IF(B2="NED",1,IF(OR(B2="ED",B2="UN"),0,""))']];records.getRange(`L2:L${last}`).fillDown();
records.getRange('M2').formulas=[['=IF(AND(K2=1,OR(B2="ED",B2="NED")),1,0)']];records.getRange(`M2:M${last}`).fillDown();
records.getRange('N2').formulas=[['=IF(K2=1,1,0)']];records.getRange(`N2:N${last}`).fillDown();
const sources=[['Source','Fields','Frozen source','Interpretation'],
 ['McShane 2016 mouse supplement','UniProt accession and ED/NED/UN','McShane_mouse_mmc1.xlsx; FinalTable_PostNorm_3T3_','ED=0; NED=1; UN only in NED-vs-non-NED sensitivity'],
 ['Supplied native archive','Native detection and covalent flags','Mouse_Native_HQ_GE.zip; model_v6.csv files','Any CCbond=True excludes; noncovalent ENT=True establishes the main positive call'],
 ['Official AlphaFold DB API','Model identity, sequence length','api/prediction/{UniProt_ID}; exact F1/v6','Full 1..N sequence required; model version and species verified'],
 ['Official AlphaFold confidence','Average protein pLDDT','AF-{ID}-F1-confidence_v6.json','Arithmetic mean of one score per residue; retain mean >=70'],
 ['Supplied knots','Knotted flag','knot_all.txt','Exact accession membership; retain isoform suffixes'],
 ['Supplied covalent lassos','Covalent flag','AlphaLasso.zip; Organism=Mouse','Remove terminal -F<number>; exclude matches in all four tables'],
 ['Frozen UniProt mouse query','Membrane and Secreted keyword flags','organism_id:10090; release 2026_03; retrieved 2026-10-10','KW-0472 and KW-0964; missing entry remains missing'],
 ['Statistical method','Crude and Model 1 estimates','NED ~ NCLE + z(protein length)','MLE logistic; sample SD (ddof=1); Fisher crude P; Wald intervals and regression P'],
 ['Conservative local-quality check','NCLE_High_only','Require Quality=High for a positive detection','Low-only detections are excluded from this check, not assigned to the negative group']];
style(records,'V1:Y10');records.getRange('V1:Y10').values=sources;header(records,'V1:Y1');
for(const [c,width]of [['V',34],['W',37],['X',62],['Y',89]])records.getRange(`${c}1:${c}10`).format.columnWidth=width;
style(res,'A1:G51');title(res,'Mouse NCLE–NED association');
res.getRange('A4').values=[['Mean pLDDT ≥70. Membrane, secreted, knotted and covalent-lasso proteins excluded.']];
const main=p.results.filter(r=>r.Term==='NCLE'&&r.NCLE_definition==='Any detected NCLE');
const high=p.results.filter(r=>r.Term==='NCLE'&&r.NCLE_definition==='High-quality detection only');
const label=a=>a.startsWith('Primary')?'NED vs ED':'NED vs ED + UN';
const tableHeader=[['Comparison','Method','N','NCLE OR','95% CI lower','95% CI upper','P-value']];
for(const [start,rr]of [[6,main],[14,high]]){
 res.getRange(`A${start}:G${start}`).values=tableHeader;header(res,`A${start}:G${start}`);
 res.getRange(`A${start+1}:G${start+4}`).values=rr.map(r=>[label(r.Analysis),r.Method==='Crude Fisher exact'?'Crude':'Length-adjusted',r.N,r.OR,r.CI_lower,r.CI_upper,r.P_value]);
 res.getRange(`C${start+1}:C${start+4}`).setNumberFormat('#,##0');res.getRange(`D${start+1}:G${start+4}`).setNumberFormat('0.000');
 for(const r of [start+1,start+3])res.getRange(`A${r}:G${r}`).format.fill=pale;
}
res.getRange('A12').values=[['Conservative check: require High-quality detection; exclude Low-only proteins']];
res.getRange('A12').format.font={name:'Arial',size:12,bold:true,color:ink};
res.getRange('A21').values=[['Counts supporting each analysis']];res.getRange('A21').format.font={name:'Arial',size:12,bold:true,color:ink};
res.getRange('A23:G23').values=[['Comparison','Detection rule','NCLE','NED','Non-NED','Total','NED fraction']];header(res,'A23:G23');
const countrows=p.counts;
res.getRange('A24:G31').values=countrows.map(r=>[label(r.Analysis),r.NCLE_definition==='Any detected NCLE'?'Any detection':'High only',r.NCLE,r.NED,r.Non_NED,r.Total,r.NED_frequency]);
for(let i=0;i<countrows.length;i++){
 const r=24+i,v=countrows[i],nc=v.NCLE_definition==='Any detected NCLE'?'C':'D',include=v.Analysis.startsWith('Primary')?'M':'N';
 const refs=`'Protein records'!$${include}$2:$${include}$${last},1,'Protein records'!$${nc}$2:$${nc}$${last},${v.NCLE},'Protein records'!$L$2:$L$${last}`;
 const guard=nc==='D'?`,'Protein records'!$C$2:$C$${last},${v.NCLE}`:'';
 res.getRange(`D${r}`).formulas=[[`=COUNTIFS(${refs},1${guard})`]];res.getRange(`E${r}`).formulas=[[`=COUNTIFS(${refs},0${guard})`]];
 res.getRange(`F${r}`).formulas=[[`=SUM(D${r}:E${r})`]];res.getRange(`G${r}`).formulas=[[`=D${r}/F${r}`]];
}
res.getRange('C24:F31').setNumberFormat('#,##0');res.getRange('G24:G31').setNumberFormat('0.0%');
for(const [resultRow,countRow]of [[7,24],[9,28],[15,26],[17,30]]){
 const a=countRow,b=countRow+1;
 res.getRange(`C${resultRow}`).formulas=[[`=SUM(F${a}:F${b})`]];
 res.getRange(`D${resultRow}`).formulas=[[`=D${a}*E${b}/(E${a}*D${b})`]];
 const se=`SQRT(1/D${a}+1/E${a}+1/D${b}+1/E${b})`;
 res.getRange(`E${resultRow}`).formulas=[[`=EXP(LN(D${resultRow})-1.95996398454005*${se})`]];
 res.getRange(`F${resultRow}`).formulas=[[`=EXP(LN(D${resultRow})+1.95996398454005*${se})`]];
}
res.getRange('A34').values=[['Protein-length effect in Model 1 (OR per one SD)']];res.getRange('A34').format.font={name:'Arial',size:12,bold:true,color:ink};
res.getRange('A36:G36').values=[['Comparison','Detection rule','N','Length OR','95% CI lower','95% CI upper','P-value']];header(res,'A36:G36');
res.getRange('A37:G40').values=p.results.filter(r=>r.Term==='z_Protein_length').map(r=>[label(r.Analysis),r.NCLE_definition==='Any detected NCLE'?'Any detection':'High only',r.N,r.OR,r.CI_lower,r.CI_upper,r.P_value]);
res.getRange('C37:C40').setNumberFormat('#,##0');res.getRange('D37:G40').setNumberFormat('0.000');
const notes=['Main NCLE rule matches the human coding: any detected native non-covalent entanglement.',
 'Model 1 controls only z-standardized protein length. Disorder and complex membership are not included.',
 'The High-only check excludes Low-only detections. Empty files denote no detected native NCLE.',
 'Crude P-values are two-sided Fisher exact. Regression P-values and 95% CIs use Wald tests.',
 'Regression and P-values are saved Python fits. Rerun scripts after changing source records.',
 'UN means unclassified. Including it as non-NED is a sensitivity analysis, not an ED assignment.'];
notes.forEach((n,i)=>res.getRange(`A${43+i}`).values=[[n]]);
for(const [c,w]of [['A',29],['B',26],['C',12],['D',16],['E',17],['F',17],['G',15]])res.getRange(`${c}1:${c}51`).format.columnWidth=w;
style(flow,'A1:F33');title(flow,'Mouse cohort construction','F');
flow.getRange('A4').values=[['Exact accession matching preserves isoform suffixes; no disorder filter is applied.']];
for(const [analysis,start]of [['Primary: NED vs ED',6],['Sensitivity: NED vs ED + UN',20]]){
 flow.getRange(`A${start}`).values=[[analysis]];flow.getRange(`A${start}`).format.font={name:'Arial',size:12,bold:true,color:ink};
 const rr=p.flow.filter(r=>r.Analysis===analysis);
 flow.getRange(`A${start+1}:F${start+1}`).values=[['Sequential step','Remaining','Removed','ED','NED','UN']];header(flow,`A${start+1}:F${start+1}`);
 flow.getRange(`A${start+2}:F${start+1+rr.length}`).values=rr.map(r=>[r.Step,r.Total,r.Removed,r.ED,r.NED,r.UN]);
 flow.getRange(`B${start+2}:F${start+1+rr.length}`).setNumberFormat('#,##0');
 flow.getRange(`A${start+1+rr.length}:F${start+1+rr.length}`).format.fill=pale;
}
flow.getRange('A1:A33').format.columnWidth=67;flow.getRange('B1:F33').format.columnWidth=14;
wb.recalculate();
const observed=res.getRange('D24:G31').values;
for(let i=0;i<countrows.length;i++){const e=countrows[i];if(observed[i][0]!==e.NED||observed[i][1]!==e.Non_NED||observed[i][2]!==e.Total)throw new Error('Workbook count mismatch');}
for(const [row,e]of [[7,main[0]],[9,main[2]],[15,high[0]],[17,high[2]]]){
 const v=res.getRange(`D${row}:F${row}`).values[0];if(Math.abs(v[0]-e.OR)>1e-10||Math.abs(v[1]-e.CI_lower)>1e-10||Math.abs(v[2]-e.CI_upper)>1e-10)throw new Error('Workbook OR or interval mismatch');
}
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:1800});console.log(errors.ndjson);
const verify=path.join(root,'verification');await fs.mkdir(verify,{recursive:true});
for(const[name,range,file]of [['Results','A1:G49','results_preview.png'],['Cohort flow','A1:F32','flow_preview.png'],['Protein records','A1:F12','records_preview.png'],['Protein records','V1:Y10','sources_preview.png']]){
 const b=await wb.render({sheetName:name,range,scale:1.3,format:'png'});await fs.writeFile(path.join(verify,file),new Uint8Array(await b.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);await output.save(path.join(root,'Mouse_NCLE_NED_association_results.xlsx'));
console.log('Workbook count/formula checks passed; XLSX exported.');
