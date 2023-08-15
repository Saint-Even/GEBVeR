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
#### import and install ####
if (!require("pacman")) install.packages("pacman")
# for system cleanup etc
pacman::p_load(sys)
#for parallel processing
pacman::p_load(parallel, foreach, doParallel)
# UpsetR for radiator
pacman::p_load(UpSetR)
# radiator for filter
pacman::p_load_gh("thierrygosselin/radiator")
# package install correction
radiator_pkg_install()

install.packages('ggplot2', dependencies = TRUE)
install.packages('ggtern', dependencies = TRUE)
#for radiator plots
pacman::p_load("ggplot2")
pacman::p_load("ggtern")


#Requires a input dir containing uncompressed .vcf files, as processed by commonMarker or intersectMarker
#Rquires an environment with up to date R version


#### Clean up for new run ####
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

####  Manually load vcfs ####
if (FALSE) {

  # move through this section to interactively select and filter files
  setwd("input")
  vcf <- file.choose(new = FALSE)
  setwd(home)
  # interactive filter
  setwd("data")
  f <- radiator::filter_rad(vcf)
  # USER: set the default filter setings below

#### auto filter, loads all *_common.vcf in input ####
  #from cassetteGBS -> commonMarker -> intersectMarker
  #command line call flows here
} else {
  # command line filter settings
  #...
  #missingness blacklist
  y outliers
  #heterozygosity blacklist
  y outlier
  #minor allele statistic
  maf 0.10
  #short distance LD threshold
  mac
  #long distance LD filter
  n

  # ...filter returns a single vcf
  setwd(home)
  setwd("data")
  f <- radiator::filter_rad(vcf)
  f <- radiator::filter_rad(vcf, interactive.filter = FALSE)
  radiator::write_vcf(f$gds, filename = "name")


  f$output
  f$gds
  rm(f)

  # copy inputs to data
  setwd(home)
  regex <- glob2rx("*_common.vcf")
  flist <- list.files("input", regex, full.names = TRUE)

  cl <- makeCluster(detectCores(logical = TRUE) - 2)
  registerDoParallel(cl)
  foreach(i = 1:length(flist)) %dopar% {
    vcf <- flist[i]
    file.copy(vcf, "data")
  }
  stopCluster(cl)
  registerDoSEQ()

  # loop filter through all
  setwd(home)
  regex <- glob2rx("*_common.vcf")
  flist <- list.files("data", regex, full.names = TRUE)

  cl <- makeCluster(detectCores(logical = TRUE) - 2)
  registerDoParallel(cl)
  foreach(i = 1:length(flist)) %dopar% {
    vcf <- flist[i]
    #...filter
    #...file handle to output
  }
  stopCluster(cl)
  registerDoSEQ()

  print("Complete")
}
##############################
# visualize
# https://knausb.github.io/vcfR_documentation/
# https://speciationgenomics.github.io/filtering_vcfs/

# save to disk
