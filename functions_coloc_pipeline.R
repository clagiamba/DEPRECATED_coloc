#####################################
# This function takes in:
#  eqtl.df -  eQTL data frame
#  biom.df -  biomarker data frame
  # eQTL:
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N ProbeID # optional: Gene.name/ensemblID
  # if output from matrixEQTL, get SE from: eQTL.data$se.beta  = eQTL.data$beta / eQTL.data$t.stat
  # biom quant
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N
  # biom case control
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N Ncases
  # MAF in at least one of the datasets
  #!! Please make sure that the betas refer to the same allele!
#  p12        -  probability of trait 1 and trait 2 
#  useBETA  - should BETA and SE be used for coloc, or only PVAL?
#  outfile   -  name of output file

## Example use: res = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=1e-6, useBETA=TRUE, outfile="summary.csv")
coloc.eqtl.biom <- function(eqtl.df, biom.df, p12=1e-6, useBETA=TRUE, outfile) {

suppressPackageStartupMessages({
  require(coloc);
});

###
  # Set Variables 
  maf_filter = 0.001 # 0.05  #MAF filter applied to datasets
  rsq_filter = 0.3 #Imputation quality filter applied to datasets

###
if ("Ncases" %in% names(biom.df)) cc=TRUE
maf.eqtl = ifelse("MAF" %in% names(eqtl.df), TRUE, FALSE)
maf.biom = ifelse("MAF" %in% names(biom.df), TRUE, FALSE)
if (!maf.eqtl & !maf.biom) stop("There is no MAF information in neither datasets")

## check all columns exist
if (useBETA) cols.eqtl = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "ProbeID", "N") else cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N")
if (!all(  cols.eqtl %in% names(eqtl.df))) stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
if (useBETA) cols.biom = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "N") else cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N")
if (cc) cols.biom = c(cols.biom, "Ncases")
if (!all(  cols.biom %in% names(biom.df))) stop("These columns are missing from the biomarker data: ", cols.biom[!cols.biom %in% names(biom.df)])

# Filter by imputation quality if column exists
info.columns <- grep( names(biom.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    biom.df = subset(biom.df, biom.df[,info.columns] > rsq_filter)
    }
info.columns <- grep( names(eqtl.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    eqtl.df = subset(eqtl.df, eqtl.df[,info.columns] > rsq_filter)
    }

# use only one of the MAFs from the two datasets
# Filter by MAF
if (maf.eqtl) {
   cols.eqtl = c(cols.eqtl, "MAF")
   eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
   }
if (!maf.eqtl & maf.biom) {
   cols.biom = c(cols.biom, "MAF")
   biom.df = subset(biom.df, biom.df$MAF > maf_filter)
   }

eqtl.df = eqtl.df[,cols.eqtl]
biom.df = biom.df[,cols.biom]
# Remove missing data
eqtl.df = eqtl.df[complete.cases(eqtl.df),]
biom.df = biom.df[complete.cases(biom.df),]

# if there is a "chr" in front of CHR column
hasChr=ifelse(any(grep("chr",eqtl.df$CHR))>0, TRUE,FALSE)
  if (!hasChr) (eqtl.df$CHR=paste("chr", eqtl.df$CHR, sep=""))
hasChr=ifelse(any(grep("chr",biom.df$CHR))>0, TRUE,FALSE)
  if (!hasChr) (biom.df$CHR=paste("chr", biom.df$CHR, sep=""))


   res.all <- data.frame()

  ############################### now start the loop
  list.probes <- unique(eqtl.df$ProbeID)

  # There can be more than one probe per gene 
  for (i in 1:length(list.probes)) {  ### find each gene with a cis-eQTL

    ProbeID = as.character(list.probes[i]) ##the character bit is important for probe names that are numbers
    region.eqtl <- subset(eqtl.df, ProbeID == as.character(list.probes[i]))
    pos.start <- min(region.eqtl$POS)
    pos.end   <- max(region.eqtl$POS)
    my.chr = unique(region.eqtl$CHR)
    
      matches <- which(biom.df$CHR==my.chr & biom.df$POS > pos.start & biom.df$POS < pos.end )
      region.biom <- biom.df[matches, ]
      # matches <- which(my_split_list[[as.character(my.chr)]]
      # region.biom <- subset(my_split_list[[as.character(my.chr)]], biom.df[matches, ])

      # Loop over each biomarker 
      # message(ProbeID, ": ", length(matches), " snps in biomarkers. From: ", pos.start, " To: ", pos.end)

      if (cc) {
          type= "cc"
          #  s = proportion of individuals that are cases (cases / N)
         region.biom$s1 = region.biom$Ncases/region.biom$N
         } else {  
         type = "quant"
         region.biom$s1=rep(0.5, length(region.biom$N)) ## This will be ignored since the type is "quant"
         } 


         # merged.data <- merge(region.biom[, c("SNPID", colname.pval, colname.N, "s1", colname.beta, colname.se)], region.eqtl, by = "SNPID")
         merged.data <- merge(region.biom, region.eqtl, by = "SNPID",  suffixes=c(".biom", ".eqtl"))
         
                  
         nsnps = nrow(merged.data)

         message(ProbeID, ": ", nsnps, " snps in both biomarker and eQTL data. From: ", pos.start, " To: ", pos.end)

         if (nsnps <= 2 ) ("There are not enough common snps in the region")
         if (nsnps > 2 ) {

         # For now run with p-values (better for cc data)
         # dataset.biom = list(snp = merged.data$SNPID, beta = merged.data[, colname.beta],varbeta = merged.data[,colname.se]^2,
         if (!useBETA) {
         dataset.biom = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.biom,
                         N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF)
         dataset.eqtl = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.eqtl,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF)
         } else {
         # ‘beta’ and ‘varbeta’
         dataset.biom = list(snp = merged.data$SNPID, beta = merged.data$BETA.biom, varbeta= (merged.data$SE.biom)^2,
                         N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF)
         dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$BETA.eqtl, varbeta= (merged.data$SE.eqtl)^2,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF)
         }
         suppressMessages(capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12)))
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])
         snp.biom <- merged.data[which.min(merged.data$PVAL.biom), "SNPID"]
         snp.eqtl <- merged.data[which.min(merged.data$PVAL.eqtl), "SNPID"]
         min.pval.biom <- min(merged.data$PVAL.biom)
         min.pval.eqtl <- min(merged.data$PVAL.eqtl)
         best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])

         res.temp = data.frame(ProbeID = ProbeID, Chr = my.chr, pos.start=pos.start, pos.end=pos.end, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, pp3, pp4)

         res.all <- rbind(res.all, res.temp)

     }
   }

   res.all <- data.frame(res.all)
   res.all$ProbeID <- as.character(res.all$ProbeID)
   res.all$snp.eqtl <- as.character(res.all$snp.eqtl)
   res.all$best.causal <- as.character(res.all$best.causal)

   res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   #outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfile, row.names = FALSE, quote = FALSE, sep = '\t')

   return(res.all)
}


