rm(list = ls()) 
library(haven)

manuf <- read.csv("Manufacturing_Data_2digit_ISIC31.csv")
conver <- read.csv("ISIC31_to_ISIC4_2digits.csv")
names(conver)[names(conver) == 'isic31_2'] <- 'isic2mode'  #renaming to facilitate merge
manuf <- merge(manuf, conver, by='isic2mode')

#restricting to manufacturing codes according to rev 4 (some rev 3 manuf industries classified out of manuf in rev 4, hence #firms, Agg Empl, etc, will not match with original data)
manuf <- subset(manuf, isic4_2mode>=10 & isic4_2mode<=33)

manuf$firm_share <- manuf$firms/sum(manuf$firms)
manuf$estab_share <- manuf$estab/sum(manuf$estab)
manuf$emp_share <- manuf$emp/sum(manuf$emp)
manuf$avsize_firms <- manuf$emp/manuf$firms
manuf$avsize_estab <- manuf$emp/manuf$estab

manuf <- subset(manuf, select=c('isic2mode','isic4_2mode','emp','firms','estab','firm_share','emp_share','estab_share','avsize_firms','avsize_estab'))
write.csv(manuf, "Manuf_Distribution_2digit_Isic4.csv")
