# RT-qPCR experimental-unit correction for v1.0.7

The author clarified that RT-qPCR comprises four biological replicates per group.
This supersedes the qPCR experimental-unit description in release v1.0.6 and its
September 29 correction notes. Fura-2 remains three wells per group on one plate;
WST-1 remains four wells from the available experiment.

The abstract, methods, results, limitations, captions, supplementary inventory,
highlights, figures and S1 workbook now describe biological qPCR replicates.
Raw Ct and relative-expression measurements are unchanged.

RT-qPCR comparisons were recomputed using two-sided Welch t-tests on Delta Ct.
Holm adjustment is applied across the four targets and separately across the two
TARDBP knockdown-versus-control comparisons. Adjusted p values for TRPC1, STIM1,
ORAI1 and ATP2A3 are 0.003962, 0.003962, 0.002609 and 0.00006766. Both TARDBP
contrasts also remain significant. Means and SEM in the plots use the preserved
relative-expression values; inference uses the raw Delta Ct values.

The calculation is reproducible with code/qpcr_biological_replicates_v121.py and
published in source_data/qpcr_biological_replicate_tests.csv. S1 includes the adjusted
comparisons. Fura-2 and WST-1 remain descriptive.
