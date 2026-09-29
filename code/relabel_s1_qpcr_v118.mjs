import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {FileBlob, SpreadsheetFile} from '@oai/artifact-tool';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const file=path.join(root,'supplementary/S1_laboratory_source_data.xlsx');
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(file));
const previewOnly=process.argv.includes('--inspect');
if (!previewOnly) {
  const readme=wb.worksheets.getItem('README');
  const updates={
    B2:'TARDBP RT-qPCR: four biological replicates per group. Raw Ct, dCt, ddCt and 2^-ddCt are retained.',
    B3:'Four target genes: raw target and GAPDH Ct values from four biological replicates per group. The target-gene assays share an RNA set; TARDBP uses a separate set.',
    B4:'Relative expression values plotted in Figure 1B. Four biological replicates per target and group.',
    B8:'Group means and SEM across biological replicates for RT-qPCR. Two-sided Welch tests on dCt; Holm correction across four targets and separately across two TARDBP contrasts. Fura-2 and WST-1 summaries describe well-to-well variation only.'
  };
  for (const [cell,value] of Object.entries(updates)) readme.getRange(cell).values=[[value]];
  const st=wb.worksheets.getItem('Summary_stats');
  const csv=await fs.readFile(path.join(root,'source_data/qpcr_biological_replicate_tests.csv'),'utf8');
  const lines=csv.trim().split('\n').map(x=>x.split(',')); const keys=lines.shift();
  const tests=lines.map(v=>Object.fromEntries(keys.map((k,i)=>[k,k==='holm_adjusted_p'?Number(v[i]):v[i]])));
  const values=st.getRange('A2:H13').values;
  for(let i=0;i<values.length;i++){
    let label='four biological replicates per group; mean and SEM',p=null;
    if(values[i][1]==='TARDBP silencing'){
      label='Welch dCt; max Holm-adjusted p of two TARDBP contrasts';
      p=Math.max(...tests.filter(x=>x.gene==='TARDBP').map(x=>x.holm_adjusted_p));
    } else if(values[i][0]==='1B' && values[i][2]==='shTDP-43'){
      label='two-sided Welch dCt; Holm-adjusted across four targets';
      p=tests.find(x=>x.gene===values[i][1].split(' ')[0]).holm_adjusted_p;
    }
    st.getRange(`G${i+2}:H${i+2}`).values=[[label,p]];
  }
  st.getRange('H2:H13').setNumberFormat('0.00E+00');
  for(const [sheet,cell] of [['TARDBP_qPCR','B1'],['Target_qPCR_Ct','C1'],['Target_qPCR_rel','C1']]){
    const h=wb.worksheets.getItem(sheet).getRange(cell);
    h.values=[['biological_replicate']]; h.format.columnWidth=22;
  }
  readme.getRange('B2:B10').format.wrapText=true;
  readme.getRange('A1:A10').format.columnWidth=26;
  readme.getRange('B2:B10').format.columnWidth=100;
  readme.getRange('B2:B10').format.rowHeight=65;
  readme.getRange('B5').format.rowHeight=180;
  for(const [col,width] of Object.entries({A:9,B:35,C:35,D:8,E:10,F:10,G:58,H:12}))
    st.getRange(`${col}1:${col}21`).format.columnWidth=width;
  st.getRange('G1:G21').format.wrapText=true;
  st.getRange('A2:H21').format.rowHeight=38;
  wb.recalculate();
  console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!',options:{useRegex:true,maxResults:20},maxChars:1000})).ndjson);
  const out=await SpreadsheetFile.exportXlsx(wb); await out.save(file);
}
for(const [sheet,range] of [['README','A1:B10'],['Summary_stats','A1:H13']]){
  const blob=await wb.render({sheetName:sheet,range,scale:1,format:'png'});
  await fs.writeFile(`/tmp/S1_${sheet}_${previewOnly?'before':'after'}.png`,new Uint8Array(await blob.arrayBuffer()));
}
console.log((await wb.inspect({kind:'table',range:'Summary_stats!G1:H13',include:'values',tableMaxRows:13,tableMaxCols:2,maxChars:1800})).ndjson);
