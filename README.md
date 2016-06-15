# coloc FOR INTEGRATIVE ANALYSIS #

* This folder contains scripts to loop thorugh eQTL and biomarker datasets

### Steps ###

* For eQTL and biomarker dataset: require data.frame in R with
*   eQTL data (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N ProbeID # optional: Gene.name/ensemblID
*   biom quant data (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N
*   biom case control (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N Ncases
*   MAF in at least one of the datasets
* Soon to add scripts for epigenetic annotations
