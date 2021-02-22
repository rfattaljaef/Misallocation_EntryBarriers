rm(list = ls()) 
setwd('~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Data_Misallocation')
library(haven)
usa <- read.csv("~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Data_Misallocation/County_Business_Patterns/SBA_2007.csv")
usa$total <- as.numeric(as.character(usa$total))           #converting from factor type to numeric type
usa$naics2002 <- as.numeric(as.character(usa$naics2002))           #converting from factor type to numeric type


conver <- read.csv("~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Data_Misallocation/County_Business_Patterns/2002 NAICS_to_ISIC_rev3.csv")
conver$isic <- as.numeric(as.character(conver$isic))  #converting from factor type to numeric type
conver$isic2 <- as.integer(conver$isic/100) # creating variable with 2 digit isic code associated with the 4 digit isic one

emp <- subset(usa, select=c('naics2002','total'), data=='empl.' ) #& naics2002>=310000 & naics2002<=340000
names(emp)[names(emp) == 'total'] <- 'emp'
firms <- subset(usa, select=c('naics2002','total'), data=='firms')
names(firms)[names(firms) == 'total'] <- 'firms'
estab <- subset(usa, select=c('naics2002','total'), data=='estab.')
names(estab)[names(estab) == 'total'] <- 'estab'  
merged <- merge(emp, conver, by='naics2002', all=TRUE)
merged <- merge(merged,firms,by.y='naics2002')
merged <- merge(merged,estab,by.y='naics2002',all=TRUE)
merged$naics2002 <- as.numeric(as.character(merged$naics2002))  #converting from factor type to numeric type
manuf <- subset(merged, merged$naics2002>310000 & merged$naics2002<=340000)
library(modeest)
library(dplyr)

#grouping naics into 2 digit isic
col <- manuf %>% group_by(naics2002) %>% mutate(
  isic2mode_int = mlv(isic2, method = "mfv",na.rm=TRUE)[1],
  isic2mode = lapply(isic2mode_int,min)                       #when multiple mode, assing the naics industry to the isic industry with lowest 2 digit
)
col$isic2mode <- as.numeric(as.character(col$isic2mode)) #from list to numeric

col <- col %>% group_by(naics2002) %>%
  summarize_all(mean,na.rm=TRUE)

#library(dplyr)
collapsed <- subset(col, select=c('isic2','naics2002','isic2mode','emp','firms','estab'))
collapsed <- collapsed %>% group_by(isic2mode) %>%
  summarize_at(c('emp','firms','estab'),sum, na.rm=TRUE)

manufcol <- subset(collapsed,isic2mode>=15 & isic2mode<=37)

write.csv(manufcol, "~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Data_Misallocation/County_Business_Patterns/Manufacturing_Data_2digit_ISIC31.csv")



