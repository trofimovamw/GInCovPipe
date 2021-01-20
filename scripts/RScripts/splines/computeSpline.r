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
if(length(args) < 8) {
  cat("\nCall the script with 4 arguments: estimatesFile reportedCasesFile label outputFile\n\n
The first input file contains a tab separated table with headers and has at least 3 columns: \n
t value variance.\n\n
The second input file is a comma separated file and contains a table with reported cases on each date. Columns must be named: \n
date new_cases \n\n
The third argument is a label, e.g. country or city.\n\n
The fourth argument is the output path. The output tables are written to the given output directory,
      which is created if it does not exist yet.\n")
  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])

#cases.list<- args[2]
table_name = args[2]
table_delim = args[3]
table_date_col = args[4]
table_active_col = args[5]
table_date_format = args[6]
country = toString(args[7])
outputFile<-file.path(args[8])

cases.table = read.table(normalizePath(table_name), header=T, sep=table_delim)
input.table = read.table(inputFile, header=T, sep = "\t")

# Define output file and output directory
fileName<-basename(outputFile)
outputDir<-normalizePath(dirname(outputFile))

if(!dir.exists(outputDir)) 
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

outputFile<-paste0(normalizePath(outputDir),"/",fileName)

### All absolute paths are set, change working directory
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

# load other r routines
source("splineRoutines.r")
source("plotRoutines.r")
source("dateFormatRoutines.r")



# change column names as they can be chosen flexibly, but we require date and new_cases
names(cases.table)[names(cases.table) == table_date_col] <- "date"
names(cases.table)[names(cases.table) == table_active_col] <- "new_cases"

#change date format to yyy-mm-dd
cases.table$date <- as.Date(cases.table$date, table_date_format)

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
cases.table$t <- as.days_since_d0(cases.table$date)

cases.spline.table <- computeSplineNewCasesTable(cases.table)

# Write spline table
write.csv(spline.table,paste0(normalizePath(outputDir),"/spline_",country,".csv"), row.names = T)
# Write reported cases table
write.csv(cases.table,paste0(normalizePath(outputDir),"/cases_",country,".csv"), row.names = T)
# Write reported cases spline table
write.csv(cases.spline.table,paste0(normalizePath(outputDir),"/cases_spline_",country,".csv"), row.names = T)
# Write estimated theta table
write.csv(input.table,paste0(normalizePath(outputDir),"/theta_",country,".csv"), row.names = T)

