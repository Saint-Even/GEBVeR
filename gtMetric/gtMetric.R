#### requirements inputs and assumptions

## requires an env with current version of R

## Expected Directory Structure

## - input
##   - ground-truth
##     - 1x .csv
##   - phenotype
##     - 1x .csv
##   - predictions
##     - nx .csv
## - data
## - output

## Input File Specifications

## Ground Truth file (.csv):
##    - Contains validated trait values for genotypes
##    - Used as the baseline for comparison with predictions
## &&& format cols and rows

## Phenotype file (.csv):
##    - same as used for BWGS training, determining the set of output predictions
##    - Contains observed trait data for genotypes
##    - rows genotypes
##    - col 1 Genotype
##    - n cols traits
##    - trait col names are here used to identify prediction files

## Prediction files (.csv in input/predictions/):
##    - Prediction files must correspond to specific traits in phenotype columns
##    - the file names must be formated like 'fileName_<trait>_theRest_ofTheFileName.csv'
##    - From the BWGS process, it is assumed that a col named 'gpred' contains the predicted trait value
##    - Contains predicted trait values (GEBVs) for genotypes
##    - the  line names in the genotype col must match some of the ground truth line names

## Key relationships:
## - Genotype col names should be consistent across Phenotype and Ground Truth files, but may be checked and corrected
## - Trait columns in phenotype and ground truth files should match, but may be checked and corrected
## - Analysis focuses on intersecting genotypes and traits across these files, before testing
## - The script aligns and compares these files to evaluate prediction accuracy for each trait.

#### setup ####
                                        # import
if (!require("pacman")) install.packages("pacman")
#p_load("languageserver") #use if language server disconnected
p_load("openxlsx") # for spreadsheet import
p_load("sys") # for file system interaction
p_load("glue") # for convent printing
p_load("tidyr")
p_load("recometrics") #for ndcg
p_load("Hmisc") # stats tools
p_load("corrplot") # corrplot
p_load("ggplot2") #plotting functions
p_load("gridExtra") #grid print multiple plots

                                        # set home
getwd()
# setwd("")
home <- getwd()

                                        # clean up
dirs <- c("data", "output")
for (d in dirs) {
  unlink(d, recursive = TRUE)
  sys::exec_wait("mkdir", d)
  setwd(d)
  sys::exec_wait("touch", ".gitkeep")
  setwd(home)
}
rm(d, dirs)

                                        # file selection
setwd(home)
setwd("input/phenotype")
flist <- list.files(pattern = "*.csv", full.names = TRUE)
if (length(flist) != 1 ){
  stop("input/phenotype must contain 1 .csv file")
}
phen_file <- flist[1]
phen <- read.csv(phen_file, header = TRUE, row.names = NULL)

setwd(home)
setwd("input/ground-truth")
flist <- list.files(pattern = "*.csv", full.names = TRUE)
if (length(flist) != 1) {
  stop("input/ground-truth must contain 1 .csv file")
}
grtr_file <- flist[1]
grtr <- read.csv(grtr_file, header=TRUE, row.names=NULL)

setwd(home)
setwd("input/predictions")
flist <- list.files(pattern = "*.csv", full.names = TRUE)
if (length(flist) <= 0) {
  stop("input/predictions must contain 1 or more .csv files")
}
## USER: optional select a single .csv gebv predictions file, override flist
if (FALSE) {
  setwd(home)
  setwd("input/predictions")
  dir()
  flist <- file.choose(new = FALSE)
}

                                        # rename grtr cols
## names(grtr)[names(grtr) == "<change-me>"] <- "<to-this>"
names(grtr)[names(grtr) == "name"] <- "Genotype"

                                        # get cols
phen_cols <- colnames(phen)
grtr_cols <- colnames(grtr)
                                        #intersect cols for expected matches
phen_only <- setdiff(phen_cols, grtr_cols)
grtr_only <- setdiff(grtr_cols, phen_cols)
shared <- intersect(grtr_cols, phen_cols)

                                        # report for confirm
print("We expect the prediction files in input/predictions to match some of the phenotype columns.")
print("Here the unmatched columns are presented to check if expected matches are being missed due to mispelling. If that is the case, rename the columns in the ground truth file, then reload the file.")
print("NOTE: the name of the column for genotypes must match!")
print("")
print("Unmatched cols in groundtruth file (rename these):")
print(grtr_only) # ground truth cols should be edited if a match is expected
print("")
print("Unmatched cols in phenotype file (to match these):")
print(phen_only) # because we know that prediction file names will match some of these
print("")
print("Matched columns:")
print(shared)

                                        #present matches
for (i in 1:length(shared)) {
  print(paste(i, ":", shared[i]))
}

## selection for validation
## head"
## ht"
## mat"
## yield"
## agr"
## twt"
## kwt"
## plump"
## thin"
## gpro"
## aa"
## dp"
## xf"
## spro"
## s_t"
## bglu3"

## USER: enter the numbers to indicate the columns containing traits. eg: traits <- c(2,3,4)
traitCols <- c(2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)
## USER: enter the number of the column with genotype names eg: gtypeCol <- 1
genoCol <- 1
genoColStr <- shared[genoCol]

intersect_selectedTraits_predFiles <- function(traitCols, shared, flist) {
  matched_files <- character()
  matched <- character()
  non_matched <- character()
# find
  for (i in traitCols) {
    trait <- shared[i]
    matching_files <- flist[grepl(trait, flist, fixed = TRUE)]
    if (length(matching_files) > 0) {
      matched_files <- c(matched_files, matching_files)
      matched <- c(matched, trait)
    } else {
      non_matched <- c(non_matched, trait)
    }
  }
# report
  print("Traits with matching prediction file")
  print(matched)
  print("Traits with no matching prediction file")
  print(non_matched)
# test
  if (length(matched) == 0) {
    stop("No prediction files were found to match selected traits")
  }
  return(matched_files)
}

matchedPredictions <- intersect_selectedTraits_predFiles(traitCols, shared, flist)

##########################
#### metric functions ####
dataFrameValues <- function(vec1, vec2, trait){
  df <- data.frame(vec1, vec2)
  tgt <- paste0(trait, "_Ground_Truth")
  tp <- paste0(trait, "_Predicted")
  colnames(df)[colnames(df)=="vec1"] <- tgt
  colnames(df)[colnames(df)=="vec2"] <- tp
  return(df)
}
calc_correlation <- function(vec1, vec2){
  spearman <- cor(vec1, vec2, method="spearman")
  pearson <- cor(vec1, vec2, method="pearson")
  kendall <- cor(vec1, vec2, method="kendall")
  v <- c(Spearman = spearman, Pearson=pearson, Kendall=kendall)
  return(v)
}
calc_basicStats <- function(vec1, vec2){
  summaryTrue = summary(vec1)
  summaryPredicted = summary(vec2)
  v <- c(Ground_Truth=summaryTrue, Predicted=summaryPredicted)
  return(v)
}
calc_performanceMetrics <- function(vec1, vec2){
  mse <- mean((vec1 - vec2)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(vec1 - vec2))
  r_sq <- cor(vec1, vec2)^2
  df <- data.frame(MSE=mse, RMSE=rmse, MAE=mae, R_SQUARED=r_sq)
  return(df)
}
calc_correlationTest <- function(vec1, vec2){
  corrTest <- rcorr(vec1, vec2)
  v <- c(r=corrTest$r[2], P=corrTest$P[2])
  return(v)
}
                                        # Scatter plot
plot_scatter <- function(vec1, vec2, trait) {
  range <- range(vec1, vec2)
  ggplot(data.frame(True = vec1, Predicted = vec2), aes(x = True, y = Predicted)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    coord_cartesian(xlim = range, ylim = range) +
    labs(title = glue("True vs Predicted Values {trait}"),
         x = "True Values",
         y = "Predicted GEBV") +
    theme_minimal()
}
                                        # Bland-Altman plot
plot_blandAltman <- function(vec1, vec2, trait) {
  mean_values <- (vec1 + vec2) / 2
  differences <- vec2 - vec1
  ggplot(data.frame(Mean = mean_values, Diff = differences), aes(x = Mean, y = Diff)) +
    geom_point() +
    geom_hline(yintercept = mean(differences), color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = mean(differences) + 1.96 * sd(differences), color = "blue", linetype = "dashed") +
    geom_hline(yintercept = mean(differences) - 1.96 * sd(differences), color = "blue", linetype = "dashed") +
    labs(title = glue( "Bland-Altman Plot {trait}"),
         x = "Mean of True and Predicted Values",
         y = "Difference (Predicted - True)") +
    theme_minimal()
}
                                        # Histogram of differences
plot_diffHisto <- function(vec1, vec2, trait) {
  differences <- vec2 - vec1
  mindif <- min(differences)
  maxdif <- max(differences)
  bin <- ((maxdif - mindif) / 10)
  ggplot(data.frame(Diff = differences), aes(x = Diff)) +
    geom_histogram(binwidth = bin, fill = "lightblue", color = "black") +
    labs(title = glue( "Histogram of Differences {trait} (Predicted - True)"),
         x = "Difference",
         y = "Frequency") +
    theme_minimal()
}
                                        #agreement matrix
plot_agreeMatrix <- function(vec1, vec2, trait) {
  common_names <- intersect(names(vec1), names(vec2))
  vec1 <- vec1[common_names]
  vec2 <- vec2[common_names]
  df <- data.frame(
    Name=common_names,
    True=vec1,
    Predicted=vec2,
    Difference= vec2 - vec1 #overestimates are positive values
  )
  df <- df[order(df$Difference),] #sort
  df$Name <- factor(df$Name, levels=df$Name) #convert name to factor to preserve order
  ggplot(df, aes(x= Name))+
    geom_col(aes(y = True, fill = Difference), alpha=0.5)+
    geom_col(aes(y=Predicted, fill=Difference), alpha=0.5)+
    geom_point(aes(y=True, x=Name))+
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0)+
    labs(
      title=glue("Agreement: {trait}"),
      x="Genotypes",
      y="Values",
      fill="Difference (Predicted - True)"
    )+
    theme_minimal()+
    theme(
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
      axis.text.y = element_text(size=6)
    )+
    coord_flip() #rotate for horizontal bars
}
#### metric functions ####
##########################

##########################
#### utility functions ####
printer <- function(data){
  sprintf("%s: %s", names(data), data)
}
calculate_metrics <- function(vec1, vec2, trait) {
  ## vec1 is True values
  ## vec2 is Predicted values
  ## assumed: both are named, names match, length matches, no NANs
                                        # list of calls
  ##data
  dfval <- dataFrameValues(vec1, vec2, trait)
  ##calcs
  cbas <- calc_basicStats(vec1, vec2)
  cper <- calc_performanceMetrics(vec1, vec2)
  ccor <- calc_correlation(vec1, vec2)
  ccot <- calc_correlationTest(vec1, vec2)
  ##plots
  psca <- plot_scatter(vec1, vec2, trait)
  pbla <- plot_blandAltman(vec1, vec2, trait)
  pdif <- plot_diffHisto(vec1, vec2, trait)
  pmat <- plot_agreeMatrix(vec1, vec2, trait)
                                        # write values used to csv file
  setwd(home)
  setwd("data")
  fname=glue("values_{trait}.csv")
  write.csv(dfval, fname, row.names=TRUE)
                                        # write text stats results to file
  data <- c(
    "Trait:", trait,
    "",
    "Basic Trait Stats:", printer( cbas ),
    "",
    "Model Performance Statistics:", printer( cper ),
    "",
    "Correlation Statistics:", printer( ccor ),
    "",
    "Correlation Test:", printer( ccot )
  )
  fname <- glue("stats_{trait}.txt")
  writeLines(as.character(data), fname)
                                        # write plots to pdf files
  fname <- glue("plots_{trait}.pdf")
  pdf(fname)
  grid.arrange(
    psca,
    pbla,
    pdif,
    pmat,
    ncol=2
  )
  dev.off()
}
#### utility functions ####
##########################

for (f in matchedPredictions) {

  # get predicton file
  setwd(home)
  setwd("input/predictions")
  pred <- read.csv(f, header = TRUE)
  # get trait from filename
  splits <- strsplit(f, "_")[[1]]
  trait <- splits[2]

  setwd(home)
  setwd("data")
  ## extract named vector from grtr and pred
  namedPredictions <- setNames(pred$gpred, pred$lines)
  namedGrtr <- setNames(grtr[[trait]], grtr[[genoColStr]])
  ## drop  any nans, report count of remaining lines
  nansPred <- sum(is.na(namedPredictions))
  nansGrtr <- sum(is.na(namedGrtr))
  anyNans <- nansPred + nansGrtr
  if (anyNans > 0) {
    print("Dropping NaN from data.")
    print(glue("count of NAN in Predictions: {nansPred}"))
    print(glue("count of NAN in GroundTruth: {nansGrtr}"))
    namedPredictions <- na.omit(namedPredictions)
    namedGrtr <- na.omit(namedGrtr)
  }

  ## intersect by names
  namedMatches <- intersect(names(namedPredictions), names(namedGrtr))
  ## report count of matching lines
  print("Count of lines in data set. Matching lines are final data set.")
  print(glue("Ground Truth Lines: {length(namedGrtr)}"))
  print(glue("Predictions Lines: {length(namedPredictions)}"))
  print(glue("Matching Lines: {length(namedMatches)}"))
  print(namedMatches)
  ## drop non match lines
  namedPredictions <- namedPredictions[namedMatches]
  namedGrtr <- namedGrtr[namedMatches]
  ## run comparison metrics args, two named vectors, trait string
  setwd(home)
  setwd("data")
  results <- calculate_metrics(namedGrtr, namedPredictions, trait)
  # move files to output
  setwd(home)
  setwd("data")
  pattern <- paste0(".*_", trait, "\\.[^.]+$")
  matching_files <- list.files(path = ".", pattern = pattern, full.names = TRUE)
  target_dir <- file.path("..", "output", trait)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  for (file in matching_files) {
    path <- file.path(target_dir, basename(file))
    file.copy(file, path)
    file.remove(file)
  }
  ## cleanup anything left in data
  setwd(home)
  setwd("data")
  pattern <- paste0(".*")
  matching_files <- list.files(path = ".", pattern = pattern, full.names = TRUE)
  for (file in matching_files) {
    file.remove(file)
  }
setwd(home)
}

exit()
#########################################
