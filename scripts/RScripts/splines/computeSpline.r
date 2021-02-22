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
for(p in c("ggplot2", "mgcv","grid","gridExtra","MASS","R0")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 8) {
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

# load other r routines
source("splineRoutines.r")
source("plotRoutines.r")
source("dateFormatRoutines.r")



# change column names as they can be chosen flexibly, but we require date and new_cases
names(cases.table)[names(cases.table) == table_date_col] <- "date"
names(cases.table)[names(cases.table) == table_active_col] <- "new_cases"

#change date format to yyy-mm-dd
# minDate to use in all plotted tables
minDate1 = min(as.Date(cases.table$date, table_date_format))
minDate2 = min(as.Date(input.table$meanBinDate, table_date_format))
minDate = min(c(minDate1,minDate2))
cases.table$date <- as.Date(cases.table$date, table_date_format)
meta.table$Collection_date <- as.Date(meta.table$Collection_date, table_date_format)

# Compute the splines and dot sizes
cat("--- Compute spline and interpolation ---\n\n")
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

# cat("--- Compute Rnaught ---\n\n")
# # Calculate spline derivatives for Rnaught plot
# infPer = 5
# gam.deriv <- computeSpline(input.table)
# t = input.table$t
# spline.deriv.table <- computeSplineDerivativeTableGratia(t,gam.deriv,infPer)
#
# # Plot R0
# outputFileR0<-paste0(outputDir,"/","rnaught_",fileName)
# plotRzero(spline.deriv.table, infPer, input.table, outputFileR0)

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


cases.table$t <- as.days_since_global_d0(cases.table$date,minDate)

cases.spline.table <- computeSplineNewCasesTable(cases.table)

### TODO: decide for which/how many time points the smoothing should be calculated and exchange
### or add to the spline table (depending if you want to keep both for the evaluation)

interp.table <-  computeInterpolation(input.table, seq(min(input.table$t), max(input.table$t)), input.table$sampleSize)
# Plot interpolated curve with 95% CI starting on global minDate
interp.table["date"] = days.as.Date(interp.table$t, minDate)
#interp.table <- interp.table[!is.na(interp.table$median) ,]
interp.table[is.na(interp.table)] = 0
# Remove rows with zeros in smooth median
interp.table <- interp.table[interp.table$smoothMedian != 0,]

# R0 package
# Wallinga and Teunis (2004)
# Generation intervals distribution - lognormal with mean = 5, sd = 1
GT <- generation.time(type = "lognormal",
                      val = c(5,1), truncate = NULL, step = 1, first.half = TRUE,
                      p0 = TRUE)
GT_obj <- R0::generation.time("empirical", val=GT$GT)
td <- est.R0.TD(as.numeric(unlist(round(interp.table$smoothMedian))),GT=GT_obj,t=days.as.Date(interp.table$t, minDate))
td.table <- data.frame(t=interp.table[1:length(td$R),]$t,value=as.vector(td$R),lower=as.vector(td$conf.int$lower),upper=as.vector(td$conf.int$upper))
outputFileWT <- paste0(normalizePath(outputDir),"/","wt04_",fileName)
plotR0Package(td.table,"WT04",outputFileWT)

# Bayesian
te <- est.R0.SB(as.numeric(unlist(round(interp.table$smoothMedian))),GT=GT_obj,t=days.as.Date(interp.table$t, minDate))
te.table <- data.frame(t=interp.table[1:length(te$R),]$t,value=as.vector(te$R),lower=as.vector(te$conf.int$CI.lower),upper=as.vector(te$conf.int$CI.upper))
outputFileE <- paste0(normalizePath(outputDir),"/","bayes_",fileName)
plotR0Package(te.table,"Bayes.",outputFileE)

# Compute q(t) = y(t)/y(t-1)
qtm <- c()
qt5 <- c()
qt95 <- c()

qtm <- c(qtm,0)
qt5 <- c(qt5,0)
qt95 <- c(qt95,0)

for (i in 2:nrow(interp.table)) {
  qtm <- c(qtm,interp.table$smoothMedian[i]/interp.table$smoothMedian[i-1])
  qt5 <- c(qt5,interp.table$smooth5[i]/interp.table$smooth5[i-1])
  qt95 <- c(qt95,interp.table$smooth95[i]/interp.table$smooth95[i-1])
}
interp.table$qtMedian <- qtm
interp.table$qt5 <- qt5
interp.table$qt95 <- qt95

# Write tables and plot
write.csv(interp.table,paste0(outputDir,"/interpolation_",country,".csv"), row.names = F, col.names = T)
outputFileInter<-paste0(normalizePath(outputDir),"/","rep_cases_interp_",fileName)
outputFileInterDots<-paste0(normalizePath(outputDir),"/","rep_cases_interp_wdots_",fileName)
plotInterpolationWithNewCases(cases.table, interp.table, input.table, meta.table, minDate, outputFileInter, outputFileInterDots)

# Qt plot
outputFileQt <- paste0(normalizePath(outputDir),"/","qt_interp_",fileName)
plotQT(interp.table,minDate,outputFileQt)
# Plot spline with daily new cases data - no negative values
# With global minDate
outputFileRC<-paste0(normalizePath(outputDir),"/","rep_cases_",fileName)
plotSplineWithNewCases(cases.table, input.table, spline.table, outputFileRC, minDate)

# Smoothed inerpolation
interp.table2 <-  computeSmoothedInterpolation(input.table, seq(min(input.table$t), max(input.table$t)), input.table$sampleSize)
interp.table2["date"] = days.as.Date(interp.table2$t, minDate)
interp.table2[is.na(interp.table2)] = 0

# Generation intervals distribution
td2 <- est.R0.TD(as.numeric(unlist(round(interp.table2$smoothMedian))),GT=GT_obj,t=days.as.Date(interp.table2$t, minDate))
td2.table <- data.frame(t=interp.table2[1:length(td2$R),]$t,value=as.vector(td2$R),lower=as.vector(td2$conf.int$lower),upper=as.vector(td2$conf.int$upper))
outputFileWT <- paste0(normalizePath(outputDir),"/","wt04_smooth_",fileName)
plotR0Package(td2.table,"WT04",outputFileWT)

# Bayesian
te2 <- est.R0.SB(as.numeric(unlist(round(interp.table2$smoothMedian))),GT=GT_obj,t=days.as.Date(interp.table2$t, minDate))
te2.table <- data.frame(t=interp.table2[1:length(te2$R),]$t,value=as.vector(te2$R),lower=as.vector(te2$conf.int$CI.lower),upper=as.vector(te2$conf.int$CI.upper))
outputFileE <- paste0(normalizePath(outputDir),"/","bayes_smooth_",fileName)
plotR0Package(te2.table,"Bayes.",outputFileE)


# Compute q(t) = y(t)/y(t-1)
qtm <- c()
qt5 <- c()
qt95 <- c()

qtm <- c(qtm,0)
qt5 <- c(qt5,0)
qt95 <- c(qt95,0)

for (i in 2:nrow(interp.table2)) {
  qtm <- c(qtm,interp.table2$smoothMedian[i]/interp.table2$smoothMedian[i-1])
  qt5 <- c(qt5,interp.table2$smooth5[i]/interp.table2$smooth5[i-1])
  qt95 <- c(qt95,interp.table2$smooth95[i]/interp.table2$smooth95[i-1])
}
interp.table2$qtMedian <- qtm
interp.table2$qt5 <- qt5
interp.table2$qt95 <- qt95

#Write tables and plot
outputFileInter<-paste0(normalizePath(outputDir),"/","rep_cases_interp_smoothed_",fileName)
outputFileInterDots<-paste0(normalizePath(outputDir),"/","rep_cases_interp_smoothed_wdots_",fileName)
write.csv(interp.table2,paste0(outputDir,"/interpolation_smooth_",country,".csv"), row.names = F, col.names = T)
plotInterpolationWithNewCases(cases.table, interp.table2, input.table, meta.table, minDate, outputFileInter, outputFileInterDots)

# Qt plot
outputFileQt <- paste0(normalizePath(outputDir),"/","qt_interp_smoothed_",fileName)
plotQT(interp.table2,minDate,outputFileQt)

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


############################# Rnaught #############################
#interp.table and interp.table2 contain already t, derivativeMedian, derivative5, derivative95 and date
#deriv.table <- data.frame(t=head(ts,-1),deriv=derivs.vec)
#With Method I
# width = 21
# interp.table$rzeroMedian <- filter(interp.table$rzeroMedian, rep(1/width,width))
# interp.table$rzero5 <- filter(interp.table$rzero5, rep(1/width,width))
# interp.table$rzero95 <- filter(interp.table$rzero95, rep(1/width,width))
#
# interp.table2$rzeroMedian <- filter(interp.table2$rzeroMedian, rep(1/width,width))
# interp.table2$rzero5 <- filter(interp.table2$rzero5, rep(1/width,width))
# interp.table2$rzero95 <- filter(interp.table2$rzero95, rep(1/width,width))

outputFileR0<-paste0(normalizePath(outputDir),"/","rzero_5days_",fileName)
plotInterpolationR0(interp.table,minDate,outputFileR0)
# With Method II
outputFileR02<-paste0(normalizePath(outputDir),"/","rzero_5days_smooth_",fileName)
plotInterpolationR0(interp.table2,minDate,outputFileR02)

outputFileDeriv<-paste0(normalizePath(outputDir),"/","deriv_5days_",fileName)
plotInterpolationDerivative(interp.table,minDate,outputFileDeriv)
# Witth Method II
outputFileDeriv2<-paste0(normalizePath(outputDir),"/","deriv_5days_smooth_",fileName)
plotInterpolationDerivative(interp.table2,minDate,outputFileDeriv2)


########################## Estimate R0 real cases ######################
# Generation intervals distribution
tdc <- est.R0.TD(as.numeric(unlist(round(cases.table$new_cases_avrg))),GT=GT_obj,t=days.as.Date(cases.table$t, minDate))
tdc.table <- data.frame(t=cases.table[1:length(tdc$R),]$t,value=as.vector(tdc$R),lower=as.vector(tdc$conf.int$lower),upper=as.vector(tdc$conf.int$upper))
outputFileWT <- paste0(normalizePath(outputDir),"/","cases_wt04_smooth_",fileName)
plotR0Package(tdc.table,"WT04",outputFileWT)
# Bayesian
tec <- est.R0.SB(as.numeric(unlist(round(cases.table$new_cases_avrg))),GT=GT_obj,t=days.as.Date(cases.table$t, minDate))
tec.table <- data.frame(t=cases.table[1:length(tec$R),]$t,value=as.vector(tec$R),lower=as.vector(tec$conf.int$CI.lower),upper=as.vector(tec$conf.int$CI.upper))
outputFileE <- paste0(normalizePath(outputDir),"/","cases_bayes_smooth_",fileName)
plotR0Package(te.table,"Bayes.",outputFileE)

# Identity DFs - method I
inp = merge(td.table, tdc.table, by.x="t", by.y="t", sort = TRUE)
outputFile <- paste0(normalizePath(outputDir),"/","cases_wt04_ident_",fileName)
plotRzeroIdentity(inp,outputFile)

# Identity DFs - method II
inp2 = merge(td2.table, tdc.table, by.x="t", by.y="t", sort = TRUE)
outputFile <- paste0(normalizePath(outputDir),"/","cases_wt04_ident_smooth_",fileName)
plotRzeroIdentity(inp2,outputFile)

# Plot both curves on one image
outputFileD <- paste0(normalizePath(outputDir),"/","cases_wt04_double_",fileName)
plotRzeroDouble(inp,minDate,outputFileD)

outputFileD2 <- paste0(normalizePath(outputDir),"/","cases_wt04_double_smooth",fileName)
plotRzeroDouble(inp2,minDate,outputFileD2)
