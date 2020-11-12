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
if(length(args) < 2) {
  cat("\nCall the script with at least 2 arguments: inputFile outputFile\n
The input file contains a tab separated table with headers and has 3 columns: t value variance\n
The output is written to the outputFile with the given output directory,
      which is created if it does not exist yet.\n
      An optional third argument trueN can be given, in order to plot the estimate versus the ground truth. For this,
      the table needs an additional column true_N\n
      An optional fourth argument is added to add a cumulative count of sequences used to calculate population
      size estimates\n
      And an optional fifth argument is added together with fourth to put a vertical marker line on some date in the plot (lockdown,
      peak reported cases etc.). Date should be in format %Y-m-%d.\n\n")
  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])

input.table = read.table(inputFile, header=T, sep = "\t")

outputFile<-file.path(args[2])

trueN<-FALSE
trueN.from.table<-FALSE
if(length(args) > 2) {
  if(args[3]=="trueN" & "trueN" %in% colnames(input.table))
    trueN=TRUE
}
if(length(args) == 6) {
  trueN.from.table<-TRUE
  
}
meta=FALSE
if(length(args) >= 5) {
  date_m<-args[5]
  var1 <- !is.na(as.Date(date_m, "%Y-%m-%d"))
  #var2
  if(var1){
    meta=TRUE
    metaFile<-normalizePath(args[4])
    date_m<-args[5]
  } else {
    cat("\n Please check if the date is in the right format and in the right place.
            Proceeding without metadata.\n\n")
  }
}

if(length(args) == 4) {
  cat("\n Please supply both meta file path and the marker line date!")
}


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

fileName<-basename(outputFile)
outputDir<-dirname(outputFile)
if(!dir.exists(outputDir))
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
outputFile<-paste0(normalizePath(outputDir),"/",fileName)



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

#TODO kann weg! auch 
#if(trueN) {
if(FALSE) {
  # if available: give true N in table and plot both
  spline.table <-addSplineValuesForTrueN(input.table, spline.table)
  outputFile_trueN <- paste0(normalizePath(outputDir),"/reportedNewCases_vs_",fileName)
  outputFile_trueNSE <- paste0(normalizePath(outputDir),"/reportedNewCases_vs_SE_",fileName)
  plotSplineWithNewCases(input.table, spline.table, outputFile_trueN)
  plotSplineWithNewCasesSE(input.table, spline.table, outputFile_trueNSE)
  if(meta) {
    meta.table = read.table(metaFile, header=T, sep = "\t")
    meta.table.freq <- as.data.frame(table(meta.table$Collection_date))
    meta.table.freq$t <- as.days_since_d0(meta.table.freq$Var1)
    plotSplineWithNewCasesSeqData(input.table, spline.table, meta.table.freq, date_m, outputFile_trueN)
    plotSplineWithNewCasesSESeqData(input.table, spline.table, meta.table.freq, date_m, outputFile_trueNSE)
  }

}

if(trueN.from.table) {
  cases.table<-read.table(casesFile, header=T, sep=',')
  cases.table$t <- as.days_since_d0(cases.table$date)
  
  # TODO Table contains negative values ??? How? for now check and substitute with 0
  #cases.table$new_confirmed[which(cases.table$new_confirmed<0)] <- 0
  
  outputFile_trueN <- paste0(normalizePath(outputDir),"/reportedNewCasesFromTable_vs_",fileName)
  outputFile_trueNSE <- paste0(normalizePath(outputDir),"/reportedNewCasesFromTable_vs_SE_",fileName)
  outputFile_trueNSESim <- paste0(normalizePath(outputDir),"/reportedNewCasesFromTable_vs_SESim_",fileName)
  
  outputFile_trueN_spline <- paste0(normalizePath(outputDir),"/reportedNewCasesFromTableSplined_",fileName)
  outputFile_trueN_splineSE <- paste0(normalizePath(outputDir),"/reportedNewCasesFromTableSplined_vs_SE_",fileName)
  
  if("new_cases" %in% colnames(cases.table)) {
    cases.table["value"] = cases.table["new_cases"]
  } else if("new_cases_smoothed" %in% colnames(cases.table)) {
    cases.table["value"] = cases.table["new_cases_smoothed"]
  }
  
  if("value" %in% colnames(cases.table)) {
    spline.cases.table <- computeSplineNewCasesTable(cases.table)
    if(meta) {
      meta.table = read.table(metaFile, header=T, sep = "\t")
      meta.table.freq <- as.data.frame(table(meta.table$Collection_date))
      meta.table.freq$t <- as.days_since_d0(meta.table.freq$Var1)
      #plotSplineWithNewCasesSeqDataRepData(input.table, cases.table, spline.table, spline.cases.table, meta.table.freq, date_m, outputFile_trueN)
      plotSplineWithNewCasesSeqDataRepDataSplined(input.table, cases.table, spline.table, spline.cases.table, meta.table.freq, date_m, outputFile_trueN_spline, F)
      plotSplineWithNewCasesSeqDataRepDataSplined(input.table, cases.table, spline.table, spline.cases.table, meta.table.freq, date_m, outputFile_trueN_splineSE, T)
  
      #plotSplineWithNewCasesSESeqDataRepData(input.table, spline.table, spline.cases.table, meta.table.freq, date_m, outputFile_trueNSE)
    }
  } else {"\n Reported cases table has no column named 'new_cases'. Rename existing columns
              or provide a new table"}
}



# Derivatives tables for intervals from BEAST pasrt of the analysis - depends
# on how many there will be

# Italy:
# 1) from start until February 9th,
# 2) from February 10th until March 9th,
# 3) from March 10th until May 10th
#
# China:
# 1) from start until January 23rd,
# 2) from January 24th until February 24th,
# 3) from February 25th until March 29th
#
# South Korea:
# 1) from start until January 24th,
# 2) from January 25th until March 6th,
# 3) from March 7th until April 17th,
# 4) from April 18th until May 29th,
# 5) from May 30th until July 10th,
# 6) from July 11th until August 11th
print(min(as.Date(input.table$meanBinDate)))
#Italy
tt = c(toString(min(as.Date(input.table$meanBinDate))),"2020-03-09","2020-03-10","2020-05-10")
#China
#tt = c(toString(min(as.Date(input.table$meanBinDate))),"2020-01-23","2020-01-24","2020-02-24","2020-02-25","2020-03-29")
# South Korea
#tt = c(toString(min(as.Date(input.table$meanBinDate))),"2020-03-06","2020-03-07","2020-04-17","2020-04-18","2020-05-29","2020-05-30","2020-07-10","2020-07-11","2020-08-11")

pseudoNumber=10^-10
input.table$value <- log(input.table$value + pseudoNumber)

spline.deriv.tables <- data.frame(t=c(),value=c(),interval=c())
for(i in (1:length(tt))) {
  if((i-1)%%2==0){

    #log transform, to avoid negative values

    gam.deriv <- computeSpline(input.table)
    print(gam.deriv)
    t1 = tt[i]
    t2 = tt[i+1]
    t1.since.d0 <- day.as.days_since_d0(min(as.Date(input.table$meanBinDate)),t1)
    t2.since.d0 <- day.as.days_since_d0(min(as.Date(input.table$meanBinDate)),t2)

    t = seq(t1.since.d0, t2.since.d0, len = t2.since.d0-t1.since.d0+1)
    spline.deriv.table <- computeSplineDerivativeTable(t,gam.deriv)
    spline.deriv.table$interval <- rep(t1,nrow(spline.deriv.table))
    spline.deriv.tables<-rbind(spline.deriv.tables,spline.deriv.table)
  } else {print("Skip")}
}
print("Spline derivatives table")
print(spline.deriv.tables)
outputFile_box <- paste0(normalizePath(outputDir),"/boxderiv_",fileName)

plotDerivativesBoxPLot(spline.deriv.tables,outputFile_box)

#write.csv(spline.deriv.tables,"/Users/mariatrofimova/Documents/GitHub/nCovPopDyn/italy.csv", row.names = T)
