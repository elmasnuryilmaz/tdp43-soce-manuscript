import { FileBlob, SpreadsheetFile } from '@oai/artifact-tool';

const path = '/Users/elmas/Desktop/MAKALE/09_YAYIN_PAKETI/supplementary/S1_laboratory_source_data.xlsx';
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(path));

const readme = workbook.worksheets.getItem('README');
const pairs = workbook.worksheets.getItem('Fura2_per_culture');
const stats = workbook.worksheets.getItem('Summary_stats');

const b5 = readme.getRange('B5').values[0][0];
if (!b5.includes('n = 3 independent cultures per group')) throw new Error('Unexpected S1 README B5');
readme.getRange('B5').values = [[b5.replace(
  'n = 3 independent cultures per group, and replicate i of both phases is the same recording.',
  'three wells per group from one culture plate, and well i of both phases is the same recording. This is one biological experiment; the well comparison is descriptive.')]];
readme.getRange('B6').values = [['Both Fura-2 phases per well with the readdition-to-release ratio of that well; one culture plate only. The sheet name is historical.']];
readme.getRange('B8').values = [['Group means and well-to-well SEM behind Figures 1 and 2; the Fura-2 and WST-1 comparisons are descriptive only.']];
pairs.getRange('B1').values = [['well']];
stats.getRange('G17').values = [['descriptive only; three wells per group on one culture plate']];
stats.getRange('H17').values = [[null]];
stats.getRange('G19').values = [['descriptive only; three wells per group on one culture plate']];
stats.getRange('H19').values = [[null]];
stats.getRange('B20:B21').values = [['readdition/ER release ratio per well'], ['readdition/ER release ratio per well']];

await workbook.recalculate();
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(path);
console.log('Relabelled S1 Fura-2 experimental unit');
