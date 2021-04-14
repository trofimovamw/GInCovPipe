#!/usr/local/bin Rscript --vanilla --slave

# clear the environment variables
rm(list = ls())


# install needed packages
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")"))))
    return(TRUE)
  install.packages(package,repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

#"ggformula"
for(p in c("ggplot2", "mgcv","grid","gridExtra","MASS","R0","scales")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  cat("\nCall the script with 4 base arguments and 5 optional arguments: estimatesFile label outputFile R0 (reportedCasesFile delim dateColumName newCasesColumnName dateFormat )\n\n
1. The estimates file contains a tab separated table with headers and has at least 3 columns: \n
t value variance.\n
2. A label for the data set, e.g. a country or city.\n
3. The output path. The output tables are written to the given directory,
      which is created if it does not exist yet.\n
4. Compute R0 - 'y' or 'n'\n\n
Optionally: \n
5. The reported cases file contains a table with reported cases on each date. \n
The separator and the column names can be chosen arbitrarily and are defined with the following parameters.\n
6. Delim gives the delimiter in the reported cases table. \n
7. The column name for the the date in the reported cases table. \n
8. The column name for the the number of cases in the reported cases table. \n
9. The format of the date, e.g. %Y-%m-%d. \n\n")
  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])
print(inputFile)
country = toString(args[2])
print(country)
outputFile<-file.path(args[3])
print(outputFile)
r0 = args[4]
table_name = args[5]
table_delim = args[6]
table_date_col = args[7]
table_active_col = args[8]
table_date_format = args[9]

# A pre-set date format for metric output
input_date_format = "%Y-%m-%d"

# Read input
input.table = read.table(inputFile, header=T, sep = "\t")
# Define output file and output directory
fileName<-basename(outputFile)
outputDir<- dirname(outputFile)
if(!dir.exists(outputDir))
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
outputDir <- normalizePath(outputDir)
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

cases.table.full = data.frame(matrix(ncol=0,nrow=0))
cases.table = data.frame(matrix(ncol=0,nrow=0))
# Some fake date
minDate1 <- integer(0)
class(minDate1) <- "Date"
if (!is.na(table_name)) {
  cases.table.full = read.table(normalizePath(table_name), header=T, sep=table_delim)
  # change column names as they can be chosen flexibly, but we require date and new_cases
  names(cases.table.full)[names(cases.table.full) == table_date_col] <- "date"
  names(cases.table.full)[names(cases.table.full) == table_active_col] <- "new_cases"
  minDate1 = min(as.Date(cases.table.full$date, table_date_format))
  #minDates <- c(minDates,minDate1)
  cases.table.full$date <- as.Date(cases.table.full$date, table_date_format)
  # Make cases table with t, new_cases, dates
  cases.table <- data.frame(new_cases=cases.table.full$new_cases,date=cases.table.full$date)
  if (!"new_cases" %in% colnames(cases.table)) {
    cat("\n Reported cases table has no column named 'new_cases'. Rename existing columns
                or provide a new table\n")
    #terminate without saving workspace
    quit("no")}
  cases.table$new_cases[cases.table$new_cases<0] <- 0
  cases.table$new_cases[is.na(cases.table$new_cases)] <- 0
  cases.table$new_cases_avrg <- filter(cases.table$new_cases, rep(1/7,7))
  cases.table <- na.omit(cases.table)
  cases.table$country <- rep(country,nrow(cases.table))

}

# minDate to use in all plotted tables
minDate2 = min(as.Date(input.table$meanBinDate,input_date_format))
minDates <- c(minDate1,minDate2)
minDate = min(minDates)
print("minDate")
print(minDate)
# Compute the splines and dot sizes
cat("--- Compute spline and interpolation ---\n\n")
input.table$t <- as.days_since_global_d0(input.table$meanBinDate,minDate)
print("h")
if (!is.na(table_name)) {cases.table$t <- as.days_since_global_d0(cases.table$date,minDate)}
pointSize <- c()
for (i in (1:nrow(input.table))) {
  # rand*(UB1-LB1) + LB1
  point <- ((10-1)*(input.table$sampleSize[i]-
            min(input.table$sampleSize)))/(max(input.table$sampleSize)-
            min(input.table$sampleSize)) + 1
  pointSize <- c(pointSize, point)
}
input.table$pointSize <- pointSize

# Replace negative number of new cases with 0 - happens if calculated from cumulative confirmed
# cases count
print(seq(min(input.table$t), max(input.table$t)))
interp.table <-  computeInterpolation(input.table, seq(min(input.table$t), max(input.table$t)), input.table$sampleSize)
# Plot interpolated curve with 95% CI starting on global minDate
interp.table["date"] = days.as.Date(interp.table$t, minDate)
print(interp.table["date"])
#interp.table <- interp.table[!is.na(interp.table$median) ,]
interp.table[is.na(interp.table)] = 0
# Remove rows with zeros in smooth median
interp.table <- interp.table[interp.table$smoothMedian != 0,]

# Write tables and plot
write.csv(interp.table,paste0(outputDir,"/interpolation_",country,".csv"), row.names = F, col.names = T)
outputFileInter<-paste0(normalizePath(outputDir),"/",fileName)
outputFileInterDots<-paste0(normalizePath(outputDir),"/","wdots_",fileName)
plotInterpolationWithNewCases(cases.table, interp.table, input.table, minDate, outputFileInter, outputFileInterDots,country)

# Write estimated theta table
write.csv(input.table,paste0(outputDir,"/theta_",country,".csv"), row.names = F)

# R0 package
# Wallinga and Teunis (2004)
# Generation intervals distribution - gamma with mean = 5, sd = 1.9
if (r0=='y') {
  GT <- generation.time(type = "gamma",
                        val = c(5,1.9), truncate = NULL, step = 1, first.half = TRUE,
                        p0 = TRUE)
  td <- est.R0.TD(as.numeric(unlist(round(interp.table$smoothMedian))),GT=GT,t=days.as.Date(interp.table$t, minDate))
  td.table <- data.frame(t=interp.table[1:length(td$R),]$t,value=as.vector(td$R),lower=as.vector(td$conf.int$lower),upper=as.vector(td$conf.int$upper))
  outputFileWT <- paste0(normalizePath(outputDir),"/","wt04_",fileName)
  plotR0Package(td.table,"WT04",outputFileWT)
}
