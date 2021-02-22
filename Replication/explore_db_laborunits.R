
# computing entry costs (ne*tauE) in units of labor, using labor share from PWT

rm(list = ls())  # clean workspace

ne = 31.65 # calibrated value of entry cost in model's undistorted stationary equilibrium
library(foreign)
library(dplyr)
setwd('~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Replication')
db <- read.csv('Entry_Costs_sector_analysis.csv')
pwt <- read.csv('pwt_labshare_2014.csv')
# selecting year 2014 only from DB
db$year <- 2014
pwt$country <- as.character(pwt$country)
pwt$country[pwt$country=='El Salvador'] <- 'Salvador'

dblab <- merge(db,pwt, by=c('country','year'))
dblab$labsh[is.na(dblab$labsh)] <- 0.6          # labor share equals 0.6 when missing in PWT
dblab <- dblab %>% mutate(
  start_L = start_gdp*(1/labsh),
  electr_L = electr_gdp*(1/labsh),
  total_L  = totalec_gdp*(1/labsh)
)
#dblab_sample <- merge(dblab, entry, by=c('country'))

dblab_sample <- subset(dblab, select=c('year','country','case','start_L','electr_L','total_L'))
write.table(dblab_sample,'db_laborunits_sector_analysys.txt')

dblabsort <- dblab_sample[order(dblab_sample$case), ]

dblabunits <- dblabsort$total_L/100*(1/ne)

write.table(dblabunits, 'db_indicator_laborunits.txt',col.names = F,row.names = F)

