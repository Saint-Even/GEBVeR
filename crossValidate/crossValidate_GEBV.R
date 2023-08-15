
##           (\.---./)
##            /.-.-.\
##           /| 0_0 |\
##          |_`-(v)-'_|
##          \`-._._.-'/      .-.
##    -~-(((.`-\_|_/-'.)))-~' <_
##           `.     .'
##             `._.'
##
##    -----~--~---~~~----~-`.-;~
##            GEBVeR



#### setup ####
if (!require("pacman")) install.packages("pacman")
# official version works on linux
p_load_current_gh("eurominister/BWGS")
#for vcf file manipulation
p_load("vcfR")
# for spreadsheet import
p_load("openxlsx")

#Requires a input dir containing uncompressed .vcf files, as processed by mergeImpute
# in dir have xlsx of BLUE/BLUP
#Requires an environment with up to date R version

#### Clean up for new run ####

setwd("/mnt/QNAP/holdens/PROC/GEBVeR/crossValidate/")
home <- getwd()
dirs <- c("data", "output")
for (d in dirs) {
  unlink(d, recursive = TRUE)
  sys::exec_wait("mkdir", d)
  setwd(d)
  sys::exec_wait("touch", ".gitkeep")
  setwd(home)
}
rm(d, dirs)

#### load test data  ####
if (FALSE) {
  # load test data for cross validation
  data(inra)
  geno_cv <- TRAIN47K
  pheno_cv <- YieldBLUE
  FIXED_cv <- "NULL"

  #!!! TEST ONLY: set NA in geno to -1 (homozygous reference)
  geno_cv[is.na(geno_cv)] <- -1
}

#### load data for cross validation: geno, pheno, FIXED ####

#geno: Matrix (n x m) of genotypes for the training population: n lines with m markers. Genotypes should be coded -1, 0, 1. Missing data are allowed and coded as NA.
#pheno: Vector (n x 1) of "phenotypes", i.e. observations or pre-processed, corrected values. This vector should have no missing values, otherwise missing values (NA) will be omitted in both pheno and geno. In a first step, bwgs.cv checks whether rownames(geno) match with names(pheno). If not the case, the common elements (intersect) are selected in both geno and pheno for further analyses. If a MAP file is provided, the selected set of markers are also sorted out in MAP.
#FIXED: A matrix of fixed effect, to be used with some methods such as those included in BGLR, MUST have same rownames as geno and coded(-1 0 1)

#interactive file selection
setwd("input")
# USER: select a .vcf genotype file
geno_cv_file <- file.choose(new = FALSE)
# USER: select a .xlsx phenotype file
pheno_cv_file <- file.choose(new = FALSE)

FIXED_cv <- "NULL"
pheno_cv <- read.xlsx(pheno_cv_file, colNames = TRUE, rowNames = FALSE)
ALLcols <- colnames(pheno_cv)
for (i in 1:length(ALLcols)) {
  print(paste(i, ":", ALLcols[i]))
}
# USER: enter a single number to indicate the column containing genotype names
# eg: genoColumn <- 1
genoColumn <- 1
# USER: enter a single number to indicate the column containing trait of interest
# eg: genoColumn <- 1
phenoColumn <- 8

#### process data ####
setwd(home)
setwd("input")
geno_cv <- vcfR::read.vcfR(geno_cv_file, verbose = FALSE)
ids <- getID(geno_cv)
u <- length(unique(ids))
l <- length(ids)
print(paste("The number of duplicated IDs:", l-u))
geno_cv <- geno_cv[!duplicated(ids, incomparables = NA)]
geno_cv <- vcfR::extract.gt(geno_cv, element = "GT", as.numeric = TRUE)
geno_cv <- as.matrix(geno_cv)
geno_cv <- t(geno_cv)

genotype <- pheno_cv[, genoColumn]
trait <- pheno_cv[, phenoColumn]
names(trait) <- genotype
pheno_cv <- trait

#test...
listGENO=rownames(geno_cv)
listPHENO <- names(pheno_cv)
inter=intersect(listGENO, listPHENO)
length(inter)
sum(is.na(geno_cv))

rm(ids, geno_cv_file, l, pheno_cv_file, genotype, trait, u, inter, listGENO, listPHENO)
rm(ALLcols, genoColumn, phenoColumn, i)
#### cross validate parameters ####
pval <- 0.5 #strict example: 0.001
nFolds <- 5
nTimes <- 5
#                1        2         3     4        5     6      7
ALLmethods <- c("GBLUP", "EGBLUP", "RR", "LASSO", "EN", "BRR", "BL",
                "BA", "BB", "BC", "RKHS", "RF", "SVM", "BRNN")
#                8     9     10    11      12    13     14
# can accept NA in snp matrix: 1, 8, 11
# can run with FIXED = NULL: 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14?

# BRNN fails could be internally reparamterized
#   expect long run times, good results to use the model free ML  method.
# SVM LASSO EN perform poorly in publication and on our data.

methods <- ALLmethods[c(1,2,3,4,5,6,7,8,9,10,11,12)]
print(paste("Cross validation method:", methods))

#### calculate cross validation ####

for (method in methods){
  print(paste("Validating:", method))
  result <- bwgs.cv(
    geno = geno_cv,
    pheno = pheno_cv,
    FIXED = FIXED_cv,
    geno.impute.method = "NULL",
    geno.reduct.method = "ANO",
    pval = pval,
    predict.method = method,
    nFolds = nFolds,
    nTimes = nTimes
  )
  name <- paste0("result_", method)
  assign(name, result)
  rm(method)
  rm(result)
  rm(name)

}

# save workspace after run
setwd(home)
setwd("data")
save.image(file="statePostCV")
# reload if needed
if (FALSE) {
  load("statePostCV")
}

#cleanup all runs
rm(methods)
rm(method)
rm(ALLmethods)
rm(nFolds)
rm(nTimes)
rm(geno_cv)
rm(FIXED_cv)
rm(pheno_cv)
rm(pval)

#### visualize methods comparisons ####

# rerun from here if dropping data >>>
results <- grep("result_", ls(), value = TRUE)
compare <- c() # model$cv
deviations <- c() # model std devs
names <- c() #clean names

# result_<method> has summary, cv, sd, MSEP, SDMSEP, bv_table
for (result in results){
  e <- eval(as.name(result))
  s <- e$cv
  compare <- cbind(compare, s)
  s <- e$SD_cv
  deviations <- append(deviations, s)
  s <- sub("result_", "", result)
  names <- append(names, s)
  rm(e)
  rm(s)
  rm(result)
}

pdf("comparison.pdf")
colnames(compare) <- names
boxplot(compare,
  xlab = "Prediction method",
  ylab = "predictive ability",
  main = "Predictive ability of methods"
  )
dev.off()

barplot(deviations,
        xlab = "Prediction method",
        ylab = "Std Dev",
        main = "Std Deviation of methods",
        names.arg = names)

print(names)

# !!!drop undesirable models from memory
if (FALSE) {
  rm(result_SVM)
  rm(result_LASSO)
  rm(result_EN)
  rm(result_RF)
  rm(result_BL)
  rm(result_RKHS)
  rm(result_EGBLUP)
  rm(result_)
}

#cleanup vis
rm(results)
rm(compare)
rm(deviations)
rm(names)

#### report on best method ####

                                        # result_<method> has summary, cv, sd, MSEP, SDMSEP, bv_table
sink("Summary.txt")
print("Summary")
result_GBLUP
sink()
