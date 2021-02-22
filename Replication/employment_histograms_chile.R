
# CONSTRUCTING CUMULATIVE DISTRIBUTION FUNCTION FOR COMPUTATION OF SHARE OF FIRMS WITH 250 WORKERS OR MORE
# the firm-level data is not provided in the repository. Please readme file for instructions on how to download it from Chile's statistical agency
# Below, I shall assume that the firm-level data is stored in the file named Enia_2014.dta




library(haven)
library(foreign)


rm(list = ls())  # clean workspace

data <- read_dta("Enia_2014.dta")
data=data.frame(data)
Nco = 1 # number of countries for which the c.d.f is constructed ( equal to 1 here, since just focusing on Chile for replication purpuses)

Mshare250vec <- matrix(0,nrow=Nco,ncol=1)

for (i in 1:Nco){

  names(data)[names(data)=='emptot'] <- 'L'
  data <- subset(data, ciiu4>=1000 & ciiu4<=3400 & L>=10)}

emp <- data$L
emp <- data.frame(emp)

Lbin <- c(10,20,50,100, 250, 500, 1000, 2500,5000,10000)
Lbinnames <- c('10-20','20-50','50-100','100-250','250-500','500-1000','1000-2500','2500-5000','5000-10000','10000+')


# Histogram
#---------------------

Mhist <- matrix(0,10,1)

for (j in 1:10) {
  if(j==10) {
    Mhist[j,1] <- sum(emp>=Lbin[j]) }
  else {
    Mhist[j,1] <- sum(emp>=Lbin[j] & emp<Lbin[j+1]) 
  }
}

Mhist[,1] <- Mhist[,1]/sum(emp>=10)  #turning into shares of agg employment

cdfM <- cumsum(Mhist)


Mshare250vec[i] <- sum(emp>=250)/sum(emp>0)
}
write.table(cbind(case,Mshare250vec),'Mshare250censuses.txt',col.names = F,row.names = F)

