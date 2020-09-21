#%j 	Day of the year
#%u Weekday as a decimal number (1–7, Monday is 1).
#%U Week of the year as decimal number (00–53) using Sunday as the first day 1 of the week (and typically with the first Sunday of the year as day 1 of week 1). The US convention.
#%W Week 00-53 with Monday as first day of the week

# return the calender week of the date (ISO 8601)
as.calWeek <- function(date) {
  return(as.numeric(format(date, "%W")))
}


# return the day of the year of the date
as.dayOfYear <- function(date) {

  return(as.numeric(format(date, "%Y%j")))
}

doy.as.Date <- function(doy) {
  return(as.Date(x = as.character(doy), format = "%Y%j"))
}

as.days_since_d0 <- function(dates) {
  dates <- as.Date(dates)
  minDate <- min(dates)
  return(as.numeric(dates - minDate))
}

days.as.Date <- function(days, minDate) {
  return(as.Date(as.Date(minDate) + days))
}

# return a representative of the calender week. For now: the last day
as.lastDayOfWeek <- function(week, year) {
  return(as.Date(x = paste0(year, "-",week,"-7"), format = "%Y-%W-%u"))
}
