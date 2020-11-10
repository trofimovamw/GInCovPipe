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
for(p in c("ggplot2", "mgcv")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  cat("\nCall the script with at least 2 arguments: inputFile outputFile\n
The input file contains a tab separated table with headers and has 3 columns: t value variance\n
The output is written to the outputFile with the given output directory,
      which is created if it does not exist yet.\n")
  #terminate without saving workspace
  quit("no")
}


cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])
input.table = read.table(inputFile, header=T, sep = "\t")

outputFile<-file.path(args[2])
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

# Take log of the estimate values
pseudoNumber=10^-10
input.table$value <- log(input.table$value + pseudoNumber)
# Make GAM
gam.deriv <- computeSpline(input.table)
spline.table <- computeSplineTable(input.table)
t1 = spline.table$t[1]
t2 = spline.table$t[length(spline.table$t)]
t = seq(t1, t2, len = t2-t1+1)
# Infectiousness period
infPer = 5

# Compute the R0 table
spline.deriv.table <- computeSplineDerivativeTable(t,gam.deriv,infPer)
write.csv(spline.deriv.table,paste0(normalizePath(outputDir),"/spline_table_deriv.csv"), row.names = T)
# Plot the R0 values
minDate <- min(as.Date(input.table$meanBinDate))
plotSplineDerivative(spline.deriv.table,minDate,normalizePath(outputDir),infPer)
