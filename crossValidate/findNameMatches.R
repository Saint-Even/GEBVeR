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
p_load("languageserver") # emacs repl etc
p_load("vcfR") # vcf file manipulation

## Requires an environment with up to date R version
## Requires a input dir, data dir, output dir
## input dir containing uncompressed .vcf files
## in dir have csv of BLUE/BLUP trait scores

#### Clean up for new run ####
#USER: ensure working directory is correct for this project
getwd()
setwd("")
home <- getwd()

# interactive file selection
setwd(home)
setwd("input")
#USER: select a .vcf genotype file
dir()
geno_file <- file.choose(new = FALSE)
#USER: select a .csv phenotype file
pheno_file <- file.choose(new = FALSE)
pheno <- read.csv(pheno_file, header=TRUE)
#pheno <- read.table(pheno_file, sep = "\t", header = TRUE)
ALLcols <- colnames(pheno)
for (i in 1:length(ALLcols)){
  print(paste(i, ":", ALLcols[i]))
}
# USER: enter a single number to indicate the column containing genotype names
# eg: genoColumn <- 1
genoColumn <- 1

#### process data ####
geno <- vcfR::read.vcfR(geno_file, verbose = FALSE)
ids <- getID(geno)
l <- length(ids)
u <- length(unique(ids))
geno <- geno[!duplicated(ids, incomparables = NA)]
geno <- vcfR::extract.gt(geno, element = "GT", as.numeric = TRUE)
geno <- as.matrix(geno)
geno <- t(geno)

g_genotypes <- rownames(geno)
p_genotypes <- pheno[,genoColumn]

shared_genotypes <- intersect(g_genotypes, p_genotypes)
g_unique_genotypes <- setdiff(g_genotypes, p_genotypes)
p_unique_genotypes <- setdiff(p_genotypes, g_genotypes)

#### report ####
print("")
print("++++++++++++++++++")
print(paste("Pheno File:", pheno_file))
print(paste("Geno File:", geno_file))
print("")
print(paste("The number of (removed) duplicated IDs in geno id col:", l-u))

# missing test
print("")
print(paste("Total NA in geno:", sum(is.na(geno))))

#mimic internal BWGS tests
print("")
print("Shared Genotypes:")
print(shared_genotypes)
print(paste("Lines shared in pheno and geno:", length(shared_genotypes)))

# print unique in each as a name check for any lines which should match but are misnamed
print("")
print(paste("Unmatched lines in Pheno:", length(p_unique_genotypes),"(consider editing where appropriate)"))
print(p_unique_genotypes)
print("")
print(paste("Unmatched lines in Geno:", length(g_unique_genotypes),"(compare to pheno above)"))
print(g_unique_genotypes)
print("++++++++++++++++++")

#### auto process data ####
setwd(home)
regex <- glob2rx("*.vcf")
location <- "input/predict"

flist <- list.files(location, regex, full.names = TRUE)

for (i in 1:length(flist)) {
  geno_file <- flist[i]

  geno <- vcfR::read.vcfR(geno_file, verbose = FALSE)
  ids <- getID(geno)
  l <- length(ids)
  u <- length(unique(ids))
  geno <- geno[!duplicated(ids, incomparables = NA)]
  geno <- vcfR::extract.gt(geno, element = "GT", as.numeric = TRUE)
  geno <- as.matrix(geno)
  geno <- t(geno)

  g_genotypes <- rownames(geno)
  p_genotypes <- pheno[,genoColumn]

  shared_genotypes <- intersect(g_genotypes, p_genotypes)
  g_unique_genotypes <- setdiff(g_genotypes, p_genotypes)
  p_unique_genotypes <- setdiff(p_genotypes, g_genotypes)

  #### report ####
  print("")
  print("++++++++++++++++++")
  print(paste("Pheno File:", pheno_file))
  print(paste("Geno File:", geno_file))
  print("")
  print(paste("The number of duplicated IDs in geno id col:", l-u))

                                        # missing test
  print("")
  print(paste("Total NA in geno:", sum(is.na(geno))))

                                        #mimic internal BWGS tests
  print("")
  print(paste("Lines shared in pheno and geno:", length(shared_genotypes)))
  print(shared_genotypes)

                                        # print unique in each as a name check for any lines which should match but are misnamed
  print("")
  print("Unmatched lines in Pheno, consider editing spelling where appropriate.")
  print(p_unique_genotypes)
  print("")
  print("Unmatched lines in Geno, compare to pheno above.")
  print(g_unique_genotypes)
  print("++++++++++++++++++")
}

#### cleanup auto ####
mem <- ls()
for (i in 1:length(mem)){
  print(paste0(i, ": ", mem[i]))
}
                                        #USER: add any items to keep, eg keep <- c("home", "gnome")
keep <- c("home")

                                        #calc non in other set
del <- setdiff(mem, keep)
nomatch <- setdiff(keep, mem)

                                        #prevent deletions due to name changes or mispellings
if (length(nomatch) > 0){
  for (i in 1:length(nomatch)){
    print(paste("An item has no match:", nomatch[i]))
  }
  stop("Check keep is correct")
}

                                        #report and execute
print(paste("Keeping:", keep))
print(paste("Deleting:", del))
rm(list=del)
rm(mem, keep, del, nomatch)
gc() #garbage collect to immediately free RAM
