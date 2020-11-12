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
for(p in c("ggplot2", "mgcv","grid","gridExtra")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  cat("\nCall the script with 4 arguments: inputFile1 inputFile2 outputFile\n
The first input file contains a tab separated table with headers and has at least 3 columns: t value variance.\n
The second input file is a comma separated file and contains a table with reported cases on each date.\n
The third argument is the grouping variable, e.g. country or city.\n
The output is written to the outputFile with the given output directory,
      which is created if it does not exist yet.\n")
  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])

input.table = read.table(inputFile, header=T, sep = "\t")

casesFile<-normalizePath(args[2])

cases.table = read.table(casesFile, header=T, sep=",")

outputFile<-file.path(args[4])

country = toString(args[3])

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

# Define output file and output directory
fileName<-basename(outputFile)
outputDir<-dirname(outputFile)

if(!dir.exists(outputDir))
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
outputFile<-paste0(normalizePath(outputDir),"/",fileName)

# Compute the splines and dot sizes
input.table$t <- as.days_since_d0(input.table$meanBinDate)
pointSize <- c()
for (i in (1:nrow(input.table))) {
  # rand*(UB1-LB1) + LB1
  point <- ((10-1)*(input.table$sampleSize[i]-
            min(input.table$sampleSize)))/(max(input.table$sampleSize)-
            min(input.table$sampleSize)) + 1
  pointSize <- c(pointSize, point)
}
input.table$pointSize <- pointSize

spline.table <- computeSplineTable(input.table)

plotSpline(input.table, spline.table, outputFile)

# Replace negative number of new cases with 0 - happens if calculated from cumulative confirmed
# cases count
if (!"new_cases" %in% colnames(cases.table)) {
  cat("\n Reported cases table has no column named 'new_cases'. Rename existing columns
              or provide a new table\n")
  #terminate without saving workspace
  quit("no")}

cases.table$new_cases[cases.table$new_cases<0] <- 0

# Write spline table
write.csv(spline.table,paste0(normalizePath(outputDir),"/spline_",country,".csv"), row.names = T)
# Write reported cases table
write.csv(cases.table,paste0(normalizePath(outputDir),"/cases_",country,".csv"), row.names = T)
# Write reported cases table
write.csv(input.table,paste0(normalizePath(outputDir),"/raw_",country,".csv"), row.names = T)
