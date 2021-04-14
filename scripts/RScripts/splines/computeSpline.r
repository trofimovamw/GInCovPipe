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
for(p in c("ggplot2", "mgcv","grid","gridExtra","MASS","R0","shadowtext","scales","plyr","lubridate")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 9) {
  cat("\nCall the script with 4 arguments: estimatesFile reportedCasesFile delim dateColumName newCasesColumnName dateFormat label outputFile\n\n
1. The estimates file contains a tab separated table with headers and has at least 3 columns: \n
t value variance.\n\n
2. The reported cases file contains a table with reported cases on each date. \n
The separator and the column names can be chosen arbitrarily and are defined with the following parameters.\n\n
3. Delim gives the delimiter in the reported cases table. \n\n
4. The column name for the the date in the reported cases table. \n\n
5. The column name for the the number of cases in the reported cases table. \n\n
6. The format of the date, e.g. %Y-%m-%d. \n\n
7. A label for the data set, e.g. a country or city.\n\n
8. The output path. The output tables are written to the given directory,
      which is created if it does not exist yet.\n
9. BEAST table for R0 values comparison.\n")
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

cases.table.full = read.table(normalizePath(table_name), header=T, sep=table_delim)
input.table = read.table(inputFile, header=T, sep = "\t")
# Just for the plots, REMOVE
meta.table = read.table("/Users/mariatrofimova/Documents/GitHub/nCovPopDyn/results/meta/meta_dates.tsv", header=T, sep = "\t")

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

remotes::install_github("wilkelab/ggtext")
library(ggtext)

# load other r routines
source("splineRoutines.r")
source("plotRoutines.r")
source("dateFormatRoutines.r")

# change column names as they can be chosen flexibly, but we require date and new_cases
names(cases.table.full)[names(cases.table.full) == table_date_col] <- "date"
names(cases.table.full)[names(cases.table.full) == table_active_col] <- "new_cases"

#change date format to yyy-mm-dd
# minDate to use in all plotted tables
minDate1 = min(as.Date(cases.table.full$date, table_date_format))
minDate2 = min(as.Date(input.table$meanBinDate, table_date_format))
minDate = min(c(minDate1,minDate2))
cases.table.full$date <- as.Date(cases.table.full$date, table_date_format)
meta.table$Collection_date <- as.Date(meta.table$Collection_date, table_date_format)

# Make cases table with t, new_cases, dates
cases.table <- data.frame(new_cases=cases.table.full$new_cases,
                          date=cases.table.full$date)

# Compute the splines and dot sizes
cat("--- Compute interpolation ---\n\n")
input.table$t <- as.days_since_global_d0(input.table$meanBinDate,minDate)
meta.table$t <- as.days_since_global_d0(meta.table$Collection_date,minDate)

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

# Replace negative number of new cases with 0 - happens if calculated from cumulative confirmed
# cases count
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

cases.table$t <- as.days_since_global_d0(cases.table$date,minDate)

cases.spline.table <- computeSplineNewCasesTable(cases.table)

interp.table <-  computeInterpolation(input.table, seq(min(input.table$t), max(input.table$t)), input.table$sampleSize)
# Plot interpolated curve with 95% CI starting on global minDate
interp.table["date"] = days.as.Date(interp.table$t, minDate)
#interp.table <- interp.table[!is.na(interp.table$median) ,]
interp.table[is.na(interp.table)] = 0
# Remove rows with zeros in smooth median
interp.table <- interp.table[interp.table$smoothMedian != 0,]
# Measures table
measure.table = read.table("/Users/mariatrofimova/Desktop/restrictions_new.csv", header =T, sep=";", as.is=T, quote = "\"",allowEscapes=TRUE)

# Write tables and plot
write.csv(interp.table,paste0(outputDir,"/interpolation_",country,".csv"), row.names = F, col.names = T)
outputFileInter<-paste0(normalizePath(outputDir),"/","rep_cases_interp_",fileName)
outputFileInterDots<-paste0(normalizePath(outputDir),"/","rep_cases_interp_wdots_",fileName)
plotInterpolationWithNewCases(cases.table, interp.table, input.table, meta.table, minDate, outputFileInter, outputFileInterDots,measure.table[which(measure.table$country==country),],country)

# Plot spline with daily new cases data - no negative values
# With global minDate
outputFileRC<-paste0(normalizePath(outputDir),"/","rep_cases_",fileName)
plotSplineWithNewCases(cases.table, input.table, spline.table, outputFileRC, minDate)

#Write tables and plot
# outputFileInter<-paste0(normalizePath(outputDir),"/rep_cases_interp_smoothed_",fileName)
# outputFileInterDots<-paste0(normalizePath(outputDir),"/","rep_cases_interp_smoothed_wdots_",fileName)
# plotInterpolationWithNewCases(cases.table, interp.table2, input.table, meta.table, minDate, outputFileInter, outputFileInterDots,measure.table[which(measure.table$country==country),],country)
# outputFileInter<-paste0(normalizePath(outputDir),"/100K_rep_cases_interp_smoothed_",fileName)
# outputFileInterDots<-paste0(normalizePath(outputDir),"/","100K_rep_cases_interp_smoothed_wdots_",fileName)
# plotInterpolationWithNewCasesPer100K(cases.table, interp.table2, input.table, meta.table, minDate, outputFileInter, outputFileInterDots,measure.table[which(measure.table$country==country),],country)

# Plot spline with daily new cases data - no negative values
outputFileRC<-paste0(normalizePath(outputDir),"/","rep_cases_",fileName)
plotSplineWithNewCases(cases.table, input.table, spline.table, outputFileRC, minDate)

# Simple output spline plot
plotSpline(input.table, spline.table, outputFile)

# Write spline table
write.csv(spline.table,paste0(outputDir,"/spline_",country,".csv"), row.names = F)
# Write reported cases table
write.csv(cases.table,paste0(outputDir,"/cases_",country,".csv"), row.names = F)
# Write reported cases spline table
write.csv(cases.spline.table,paste0(outputDir,"/cases_spline_",country,".csv"), row.names = F)
# Write estimated theta table
write.csv(input.table,paste0(outputDir,"/theta_",country,".csv"), row.names = F)

# Cases with pseudocount
cases.table$new_cases_pseudo <- cases.table$new_cases_avrg+1

# Generation time distribution
GT <- generation.time(type = "gamma",
                      val = c(5,1.9), truncate = NULL, step = 1, first.half = TRUE,
                      p0 = TRUE)
tdc <- est.R0.TD(as.numeric(unlist(round(cases.table$new_cases_pseudo))),GT=GT,t=days.as.Date(cases.table$t, minDate))
tdc.table <- data.frame(t=cases.table[1:length(tdc$R),]$t,value=as.vector(tdc$R),lower=as.vector(tdc$conf.int$lower),upper=as.vector(tdc$conf.int$upper))
outputFileWT <- paste0(normalizePath(outputDir),"/","cases_wt04_smooth_",fileName)
plotR0Package(tdc.table,"WT04",outputFileWT)

###############################################################################################################
#BEAST comparison
# R0 package
# Wallinga and Teunis (2004)
# Generation intervals distribution - lognormal with mean = 5, sd = 1.9
interp.table$smoothMedian <- interp.table$smoothMedian+1
td <- est.R0.TD(as.numeric(unlist(round(interp.table$smoothMedian))),GT=GT,t=days.as.Date(interp.table$t, minDate))
tdc.minDate <- td$begin
print(interp.table)
print(tdc.minDate)
print(minDate)
print(minDate1)
print(minDate2)
td.table <- data.frame(t=interp.table[1:length(td$R),]$t,value=as.vector(td$R),lower=as.vector(td$conf.int$lower),upper=as.vector(td$conf.int$upper))
outputFileWT <- paste0(normalizePath(outputDir),"/","wt04_",fileName)
plotR0Package(td.table,"WT04",outputFileWT)

# Bin the values - R0s
values.vec <- c()
#Switzerland
#border.dates <- as.Date(c(toString(minDate),"2020-03-16","2020-05-11","2020-06-15","2020-08-20","2020-10-19","2021-03-30"))
#Victoria
border.dates <- as.Date(c(toString(minDate),"2020-03-16","2020-05-11","2020-07-19","2020-09-13"))#,"2020-10-26","2021-03-30"))
#Denmark
#border.dates <- as.Date(c("2020-03-02","2020-03-18","2020-04-28","2020-06-15","2020-08-22","2020-10-26","2020-12-16","2021-02-09"))
# Scotland
#border.dates <- as.Date(c("2020-03-04","2020-03-24","2020-05-29","2020-07-15","2020-09-23","2020-11-17","2021-01-05","2021-02-02"))
td.table$date <- days.as.Date(td.table$t, minDate)
group.means <- c()
# 16.03-11.05-15.06-20.08-19.10
for (i in (1:nrow(td.table))){
  for (j in (2:length(border.dates))){
    if ((td.table$date[i]>=border.dates[j-1]) & (td.table$date[i]<border.dates[j])){
      values.vec <- c(values.vec,(j-1))
      group.means <- c(group.means,(border.dates[j]-border.dates[j-1])/2)
    }
  }
}
print("Here")
td.table$group <- values.vec
print("Here")
#mean.tdc.table <- aggregate(tdc.table$value, list(tdc.table$group), FUN = 'quantile', probs = 0.05)
td.table.95 <- ddply(td.table, "group", summarise, WQ95 = quantile(value, .95))
print("Here")
td.table.5 <- ddply(td.table, "group", summarise, WQ50 = quantile(value, .5))
print("Here")
td.table.05 <- ddply(td.table, "group", summarise, WQ05 = quantile(value, .05))
#
print("Here")
beast.table <- read.table(args[9], header =T, sep=" ", as.is=T, quote = "\"",allowEscapes=TRUE)
beast.table.95 <- ddply(beast.table, "interval", summarise, WQ95 = quantile(R, .95))
beast.table.5 <- ddply(beast.table, "interval", summarise, WQ50 = quantile(R, .5))
beast.table.05 <- ddply(beast.table, "interval", summarise, WQ05 = quantile(R, .05))

rzero.table <- data.frame(period=td.table.5$group,est.median=td.table.5$WQ50,
                          est.lower=td.table.05$WQ05, est.upper=td.table.95$WQ95,
                          beast.median=beast.table.5$WQ50, date=border.dates[1:(length(border.dates)-1)],
                          beast.lower=beast.table.05$WQ05, beast.upper=beast.table.95$WQ95)
rzero.table <- rbind(rzero.table, rzero.table[nrow(rzero.table),])
rzero.table$date[nrow(rzero.table)] <- max(td.table$date)
print(rzero.table)
outputFileR0 = paste0(normalizePath(outputDir),"/","rzero_beast_",fileName)
plotR0BEAST(rzero.table,td.table,outputFileR0,minDate,country)
