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
#p_load_current_gh("eurominister/BWGS")
p_load("BWGS")
#for vcf file manipulation
p_load("vcfR")
# for spreadsheet import
p_load("openxlsx")
# for file name manipulation
p_load(tools)
# for plotting
p_load("ggplot2")
p_load("viridis")
# for overlay plot
p_load("reshape")
p_load("plyr")


#Requires an input dir containing uncompressed .vcf files, as processed by merge tools.
# in dir have xlsx of BLUE estimates of traits
# in dir input a dir predict with .vcf files can be loaded automatically and processed
# Requires an environment with up to date R version

#### Clean up for new run ####
getwd()
setwd("")
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

#### load test data ####
if (FALSE) {
  # load test data for prediction
  data(inra)
  geno_train <- TRAIN47K
  pheno_train <- YieldBLUE
  FIXED_train <- "NULL"

  geno_target <- TARGET47K
  FIXED_target <- "NULL"

  #!!! TEST ONLY: set NA in geno to -1 (homozygous reference)
  geno_...[is.na(geno_cv)] <- -1
}

#### load data for prediction: geno_train, pheno_train, FIXED_train, geno_target, FIXED_target ####

#geno_train: Matrix (n x m) of genotypes for the training population: n lines with m markers. Genotypes should be coded as -1, 0, 1, NA. Missing data are allowed and coded as NA.
#pheno_train: Vector (n x 1) of phenotype for the training phenotypes. Named vector, where each genotyope name corresponds to a value. This vector should have no missing values. Otherwise, missing values (NA) will be omitted in both pheno_train and geno_train.
#FIXED_train: A matrix of fixed effect for training, to be used with some methods such as those included in BGLR, MUST have same rownames as geno and coded(-1 0 1)
#geno_target: Matrix (z x m) of genotypes for the target population: z lines with the same m markers as in geno_train. Genotypes should be coded as -1, 0, 1, NA. Missing data are allowed and coded as NA. Other arguments are identical to those of bwgs.cv, except pop_reduct_method, nTimes and nFolds, since the prediction is run only once, using the whole training population for model estimation, then applied to the target population.
#FIXED_target: A matrix of fixed effect for targeting, to be used with some methods such as those included in BGLR, MUST have same rownames as geno and coded(-1 0 1)

#### interactive selection ####
#USER: move through this section to interactively select training data
setwd(home)
setwd("input")
dir()
geno_train_file <- file.choose(new = FALSE)
pheno_train_file <- file.choose(new = FALSE)
FIXED_train <- "NULL"
FIXED_target <- "NULL"

#### get all *.vcf target files from input/predict ####
setwd(home)
regex <- glob2rx("*.vcf")
flist <- list.files("input/predict", regex, full.names = TRUE)

#### interactive single target selection ####
if (FALSE) {
  # overwrites flist with one selected file
  setwd(home)
  flist <- file.choose(new = FALSE)
}

#### process data ####

# set geno_train from file
setwd(home)
setwd("input")
geno_train <- vcfR::read.vcfR(geno_train_file, verbose = FALSE)
ls()

setwd(home)
setwd("data")
ids <- vcfR::getID(geno_train)
u <- length(unique(ids))
l <- length(ids)
print(paste("The number of duplicated IDs:", l-u))
geno_train <- geno_train[!duplicated(ids, incomparables = NA)]
geno_train <- vcfR::extract.gt(geno_train, element = "GT", as.numeric = TRUE)
geno_train <- as.matrix(geno_train)
geno_train <- t(geno_train)
ls()

                                        # set pheno_train vector
setwd(home)
setwd("input")
## pheno_train <- read.xlsx(pheno_train_file)
&&& change to read csv
pheno_train <- read.delim(pheno_train_file, sep =",")

setwd(home)
setwd("data")
print("Use the column numbers to identify the columns you want to process.")
print("Enter the numbers for genoColumn, and traitColumns")
ALLcols <- colnames(pheno_train)
for (i in 1:length(ALLcols)){
  print(paste(i, ":", ALLcols[i]))
}

                                        # USER: enter a single number to indicate the column containing genotype names
                                        # eg: genoColumn <- 1
genoColumn <- 1
                                        # USER: enter comma separated column numbers to process for prediction
                                        # eg: traitColumns <- c(3,5,9)
traitColumns <- c(2,3)

#### prediction loop of traits and files ####
                                        #get user set of traits and loop on set
for (traitCol in traitColumns){
  traitName <- ALLcols[traitCol]

                                        #reset pheno-train vector
  setwd(home)
  setwd("input")
  #pheno_train <- read.xlsx(pheno_train_file)
  pheno_train <- read.delim(pheno_train_file, sep =",")

  setwd(home)
  setwd("data")
  ALLcols <- colnames(pheno_train)
  genotype <- pheno_train[,genoColumn]

                                        #extract one col
  trait <- pheno_train[,traitCol]
  names(trait) <- genotype
  pheno_train <- trait

  # optional name match testing, use fineNameMatches.R
  test <- TRUE
  stop <- FALSE

  if (test) {

    #test exactly like internal BWGS tests
    listGENO <- rownames(geno_train)
    listPHENO <- names(pheno_train)
    inter <- intersect(listGENO, listPHENO)
    print(paste("Lines shared in pheno and geno:", length(inter)))

    # missing test
    print(paste("Total NA in geno:", sum(is.na(geno_train))))

    # print unique in each as a name check for any lines which should match but are misnamed
    gU <- setdiff(listGENO, listPHENO)
    pU <- setdiff(listPHENO, listGENO)
    print("Unmatched lines in Pheno, consider editing spelling where appropriate.")
    print (pU)
    print("Unmatched lines in Geno, compare to pheno above.")
    print(gU)

    if (stop) {
      print("Stopping early, set stop to false to continue.")
      exit(0)
    }
  }

  #### predict parameters ####
  # pval should be set to same as in cross validation runs!
  pval <- 0.45 #strict example: 0.001
  #                1        2         3     4        5     6      7
  ALLmethods <- c("GBLUP", "EGBLUP", "RR", "LASSO", "EN", "BRR", "BL",
                  "BA", "BB", "BC", "RKHS", "RF", "SVM", "BRNN")
  #                8     9     10    11      12    13     14
  # can accept NA in snp matrix: 1, 8, 11
  # can run with FIXED = NULL: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13?, 14
  # ...requires FIXED matrix

  method <- ALLmethods[1] # select one best methods from cross validation.
  print(paste("Predicting with method:", method))

  #### loop through target files ####

  for (i in 1:length(flist)) {
    geno_target_file <- flist[i]
    #... add trait to reporting
    print(paste("Predicting for:", traitName, ":", geno_target_file))
    # set geno_target from file
    setwd(home) # because file path starts at home
    geno_target <- vcfR::read.vcfR(geno_target_file, verbose = FALSE)
    setwd("data") # now drop down
    ids <- vcfR::getID(geno_target)
    u <- length(unique(ids))
    l <- length(ids)
    print(paste("The number of duplicated target IDs:", l-u))
    geno_target <- geno_target[!duplicated(ids, incomparables = NA)]
    geno_target <- vcfR::extract.gt(geno_target, element = "GT", as.numeric = TRUE)
    geno_target <- as.matrix(geno_target)
    geno_target <- t(geno_target)
    predicted <- bwgs.predict(
      geno_train = geno_train,
      pheno_train = pheno_train,
      FIXED_train = FIXED_train,
      geno_target = geno_target,
      FIXED_target = FIXED_target,
      MAP = "NULL",
      geno.impute.method = "NULL",
      geno.reduct.method = "ANO",
      pval = pval,
      predict.method = method
    )

    #predicted cols: gpred, gpredSD, CD, lines
    predicted <- as.data.frame(predicted)
    predicted <- predicted[!(row.names(predicted)) %in% "empty",]
    predicted["lines"] <- rownames(predicted)
    predicted <- predicted[order(predicted$gpred),]

    setwd(home)
    setwd("output")
    n <- file_path_sans_ext(basename(geno_target_file))
    name <- paste0("prediction_", traitName, "_", n)

    #sorted pred
    ggplot(predicted, aes(x=reorder(lines,gpred), y=gpred, fill=CD)) +
      geom_point(size=3, colour="black", shape=23) +
      scale_fill_viridis(option="turbo", begin=0.9, end=0.3, direction=1) +
      ggtitle("Predictions, Sorted:", subtitle=name) +
      theme(plot.subtitle=element_text(size=8))

    name_spec <- paste0(name, "_plotPred.pdf")
    ggsave(name_spec)

    #sorted pred with sd
    ggplot(predicted, aes(x=reorder(lines,gpred), y=gpred, fill=CD)) +
      geom_pointrange(aes(ymin=gpred-gpredSD, ymax=gpred+gpredSD),
                      fatten=7, colour="black", shape=23) +
      scale_fill_viridis(option="turbo", begin=0.9, end=0.3, direction=1) +
      ggtitle("Predictions and Std.Dev, Sorted:", subtitle=name) +
      theme(plot.subtitle=element_text(size=8))

    name_spec <- paste0(name, "_plotPredSD.pdf")
    ggsave(name_spec)

    #Hist of predictions
    ggplot(predicted, aes(x=gpred, y=after_stat(density))) +
      geom_histogram(binwidth=10, colour="lightblue", fill="white") +
      geom_vline(aes(xintercept=mean(gpred)), colour="blue", linetype="dashed") +
      geom_density(alpha=.5, fill="blue") +
      ggtitle("Distribution of Predicted Values", subtitle=name) +
      theme(plot.subtitle=element_text(size=8))

    name_spec <- paste0(name, "_plotDistOfPred.pdf")
    ggsave(name_spec)

    #Hist of training data
    train <- pheno_train
    names(train) <- NULL
    train <- data.frame(train)

    ggplot(train, aes(x=train, y=after_stat(density))) +
      geom_histogram(binwidth=10, alpha=0.2,
                     colour="palegreen", fill="white") +
      geom_vline(aes(xintercept=mean(train)), colour="springgreen4",
                 linetype="dashed") +
      geom_density(alpha=.4, fill="springgreen4") +
      ggtitle("Distribution of Training Values", subtitle=name) +
      theme(plot.subtitle=element_text(size=8))

    name_spec <- paste0(name, "_plotDistOfTrain.pdf")
    ggsave(name_spec)

    #Hist overlay

    pred <- predicted$gpred
    train <- pheno_train
    names(train) <- NULL
    w <- cbind(pred, train)
    both <- melt(w)
    both$X1 <- NULL
    names(both)[1] <- "group"

    means <- ddply(both, "group", summarise, groupMean=mean(value))
    groupColours <- c(train = "springgreen4", pred = "blue")

    ggplot(data=both, aes(x=value, y=after_stat(density),
                          colour=group, fill=group)) +
      geom_vline(data=means, aes(xintercept=groupMean, colour=group),
                 linetype="dashed") +
      geom_density(alpha=0.5) +
      ggtitle("Distribution and Mean for Train and Predict Data", subtitle=name) +
      theme(plot.subtitle=element_text(size=8)) +
      scale_colour_manual(values=groupColours) +
      scale_fill_manual(values=groupColours)

    name_spec <- paste0(name, "_plotDistOverlay.pdf")
    ggsave(name_spec)

    #write report
    name_spec <- paste0(name, "_report.csv")
    write.csv(predicted, name_spec)

  } #end for each .vcf in input/predict

} # end for each trait in traits
#### END prediction loop of traits and files ####


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
if (length(nomatch) > 0) {
  for (i in 1:length(nomatch)) {
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

exit(0)
