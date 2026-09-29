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
    B2:'TARDBP RT-qPCR: four technical measurements of the same biological sample per group. Raw Ct, dCt, ddCt and 2^-ddCt are retained. Recorded assay dates do not identify independent cultures.',
    B3:'Four target genes: raw target and GAPDH Ct values. Four measurements per group are technical repeats of the same biological sample, not independent RNA isolations.',
    B4:'Relative expression values plotted in Figure 1B. Four technical measurements per target and group; descriptive only.',
    B8:'Group means and technical-repeat or well-to-well SEM. All laboratory comparisons are descriptive; no inferential p values are reported.'
  };
  for (const [cell,value] of Object.entries(updates)) readme.getRange(cell).values=[[value]];
  const st=wb.worksheets.getItem('Summary_stats');
  const values=st.getRange('A2:H21').values;
  for(let i=0;i<values.length;i++){
    if(['1A','1B'].includes(values[i][0])){
      st.getRange(`G${i+2}`).values=[['descriptive only; four technical RT-qPCR repeats per group']];
      st.getRange(`H${i+2}`).values=[[null]];
    }
  }
  for(const [sheet,cell] of [['TARDBP_qPCR','B1'],['Target_qPCR_Ct','C1'],['Target_qPCR_rel','C1']])
    wb.worksheets.getItem(sheet).getRange(cell).values=[['technical_repeat']];
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
