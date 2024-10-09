#### requirements inputs and asumptions

## requires an env with current version of R

## Input File Specification

## Training Phenotype file (.csv):
##    - Contains observed trait data for genotypes for the BWGS training, determines its output
##    - its trait cols are used to indentify prediction files
## phen0 is trainig data


## Ground Truth Phenotype file (.csv):
##    - Contains validated trait values for genotypes
##    - Used as the baseline for comparison with predictions
## phen1 is validation data

this requires named overlap between the two data sets, this is untested.

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
setwd("")
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

                                        # interactive file selection
setwd(home)
setwd("input")
dir()

##USER: select a .csv phenotype file, this should be the same as used for predictions in BWGS
phen0_file <- file.choose(new = FALSE)
##USER: select a .csv groundTruth file
phen1_file <- file.choose(new = FALSE)

phen0 <- read.csv(phen0_file, header=TRUE, row.names = NULL)
phen1 <- read.csv(phen1_file, header=TRUE, row.names=NULL)


                                        # rename phen1 cols

if(FALSE){
  names(phen1)[names(phen1) == "name"] <- "Genotype"
}

                                        # get cols
phen0_cols <- colnames(phen0)
phen1_cols <- colnames(phen1)

                                        #intersect cols for expected matches
phen0_only <- setdiff(phen0_cols, phen1_cols)
phen1_only <- setdiff(phen1_cols, phen0_cols)
shared <- intersect(phen1_cols, phen0_cols)

                                        # report for confirm
print("")
print("Unmatched cols in phen1 file (rename these):")
print(phen1_only) #cols should be edited if a match is expected

print("")
print("Unmatched cols in phen0 file (to match these):")
print(phen0_only) # because we know that prediction file names will match some of these

print("")
print("Matched columns:")
print(shared)
                                        #present matches

for (i in 1:length(shared)){
  print(paste(i, ":", shared[i]))
}

## selection for correlation check
## [1] "2 : head"
## [1] "3 : ht"
## [1] "4 : mat"
## [1] "5 : yldp"
## [1] "6 : yield"
## [1] "7 : agr"
## [1] "8 : plump"
## [1] "9 : thin"
## [1] "10 : twt"
## [1] "11 : kwt"
## [1] "12 : gpro"
## [1] "13 : aa"
## [1] "14 : dp"
## [1] "15 : xf"
## [1] "16 : spro"
## [1] "17 : s_t"
## [1] "18 : bglu3"

## USER: enter the numbers to indicate the columns containing traits. eg: traits <- c(2,3,4)
traitCols <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)
## USER: enter the number of the column with genotype names eg: gtypeCol <- 1
genoCol <- 1
genoColStr <- shared[genoCol]

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


#### loop on traits
for ( i in traitCols){
                                        #get trait from shared traits
  trait <- shared[i]

      setwd(home)
      setwd("data")

      ## extract named vectors
      namedphen0 <- setNames(phen0[[trait]], phen0[[genoColStr]])
      namedphen1 <- setNames(phen1[[trait]], phen1[[genoColStr]])
      ## drop  any nans, report count of remaining lines
      nansP0 <- sum(is.na(namedphen0))
      nansP1 <- sum(is.na(namedphen1))
      anyNans <- nansP0+nansP1
      if(anyNans > 0) {
        print("Dropping NaN from data.")
        print(glue("count of NAN in Phen0: {nansP0}"))
        print(glue("count of NAN in Phen1: {nansP1}"))
        namedphen0 <- na.omit(namedphen0)
        namedphen1 <- na.omit(namedphen1)
      }

      ## intersect by names
      namedMatches <- intersect(names(namedphen0), names(namedphen1))
      ## report count of matching lines
      print("Count of lines in data set. Matching lines are final data set.")
      print(glue("phen0 Lines: {length(namedphen0)}"))
      print(glue("phen1 Lines: {length(namedphen1)}"))
      print(glue("Matching Lines: {length(namedMatches)}"))
      print(namedMatches)
      ## drop non match lines
      namedphen0 <- namedphen0[namedMatches]
      namedphen1 <- namedphen1[namedMatches]
      ## run comparison metrics args, two named vectors, trait string
      setwd(home)
      setwd ("data")
      results <- calculate_metrics(namedphen1, namedphen0, trait)
      # move files to output
      setwd(home)
      setwd ("data")
      pattern <- paste0(".*_", trait, "\\.[^.]+$")
      matching_files <- list.files(path=".", pattern=pattern, full.names=TRUE)
      target_dir <- file.path("..", "output", trait)
      if(!dir.exists(target_dir)) {
        dir.create(target_dir, recursive=TRUE)
      }
      for(file in matching_files){
        path <- file.path(target_dir, basename(file))
        file.copy(file, path)
        file.remove(file)
      }
      ## cleanup anything left in data
      setwd(home)
      setwd ("data")
      pattern <- paste0(".*")
      matching_files <- list.files(path=".", pattern=pattern, full.names=TRUE)
      for(file in matching_files){
        file.remove(file)
      }
}#end trait
