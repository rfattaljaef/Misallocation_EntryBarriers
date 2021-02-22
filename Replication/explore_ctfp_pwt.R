
rm(list = ls())  # clean workspace

# prior to executing the code, please refer to https://cran.r-project.org/web/packages/pwt9/index.html for downloading the Penn World Tables Database 9.0.
# In this code, this file is being labeled pwt90.csv

library(foreign)
setwd('~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Replication')
pwt <- read.csv('pwt90.csv')

pwt_sh <- subset(pwt, year==2014 ,select=c('country','countrycode','year','ctfp'))
pwt_sh$country <- as.character(pwt_sh$country)
pwt_sh$country[pwt_sh$country=='El Salvador'] <- 'Salvador'  # relabeling to match country name elsewhere
write.table(pwt_sh,'ctfp_pwt.dat',col.names=T, row.names=F, sep=',')
