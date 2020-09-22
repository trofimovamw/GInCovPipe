#!/usr/local/bin Rscript --vanilla --slave

# clear the environment variables
rm(list = ls())


# install needed packages
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")"))))
    return(TRUE)

  install.packages(package)
  return (eval(parse(text=paste("require(",package,")"))))
}

#"ggformula"
for(p in c("ggplot2", "mgcv","grid")) {
  dynamic_require(p)
}


#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  cat("\nCall the script with at least 2 arguments: inputFile outputFile\n
The input file contains a tab separated table with headers and has 3 columns: t value variance\n
The output is written to the outputFile with the given output directory,
      which is created if it does not exist yet.\n
      An optional third argument trueN can be given, in order to plot the estimate versus the ground truth. For this,
      the table needs an additional column true_N\n\n")

  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])
input.table = read.table(inputFile, header=T, sep = "\t")

trueN<-FALSE
if(length(args) == 3) {
  if(args[3]=="trueN" & "trueN" %in% colnames(input.table))
    trueN=TRUE
}

outputFile<-file.path(args[2])
print(outputFile)
fileName<-basename(outputFile)
outputDir<-dirname(outputFile)
if(!dir.exists(outputDir))
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
outputFile<-paste0(normalizePath(outputDir),"/",fileName)

# set working directory to call other R Scripts
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

this.dir <- getSrcDirectory(function(x) {x})

if (rstudioapi::isAvailable()) {
  this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}else {
  this.dir <- getScriptPath()
}
setwd(this.dir)
#getwd()

# load other r routines
source("splineRoutines.r")
source("plotRoutines.r")
source("dateFormatRoutines.r")
#source("dateFormatRoutines.r")

input.table$t <- as.days_since_d0(input.table$meanBinDate)
spline.table <- computeSplineTable(input.table)
#outputFile = paste0(outputpath, "/estimatedTheta_splineFit.pdf")
plotSpline(input.table, spline.table, outputFile)


if(trueN) {
  # for test purposes: give true N in table and plot both
  spline.table <-addSplineValuesForTrueN(input.table, spline.table)
  ratio <-computeRatio(spline.table$value, spline.table$value_trueN)
  outputFile_ratio <- paste0(normalizePath(outputDir),"/ratio_",fileName)
  outputFile_sampsize <- paste0(normalizePath(outputDir),"/sample_size_",fileName)

  plotRatio(ratio, outputFile_ratio)

  outputFile_trueN <- paste0(normalizePath(outputDir),"/reportedNewCases_vs_",fileName)
  #outputFile = paste0(outputDir, "/esti_vs_real_PopulationSize_splineFit.pdf")
  plotSplineWithNewCases(input.table, spline.table, outputFile_trueN)
  plotCummulativeSampleSize(input.table,outputFile_sampsize)
}
