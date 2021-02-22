
# loading results and data for replication of figures in the paper 

rm(list = ls())  # clean workspace

library(foreign)
library(ggplot2)
library(ggrepel)
library(gridExtra)

#*******************************************************************************************************************************************
# Constructing Support of Productivity Distribution and Probability Distribution Funciton for Productivity Distribution of Entrants
#*******************************************************************************************************************************************

nty=180              # number of productivity types
jump=0.25            # jump in binomial process for idiosyncratic productivity, calibrated to match StDev of the log of employment growth rates in the US
rho = 4              # elasticity of substitution
a=1-1./(rho-1.) 

ZZ <- matrix(0,ncol=1, nrow=nty)   # vector of idiosyncratic productivity
ZZ[1]=0.1 
for (i in 2:nty) {
  ZZ[i]=ZZ[i-1]*(exp((1.-a)*jump))
}

EZZ=ZZ^(1./(1.-a)) 
z=log(EZZ)


#*******************************************************************************************************************************************
# LOADING
#*******************************************************************************************************************************************

setwd('~/Dropbox/My Research@Dropbox/Misallocation_EntryCosts/Replication')  #set to working directory where files replications files are downloaded
data <- read.table("ResMasterAvs10indWeightedRegTr5.dat") # results from various stationary equilibriums (by country, and by distortion combination)
avs <- read.csv( "Av_Size_bySector_All_short_WeightedReg_tr5.csv")  # average size data
gdp <- read.csv( "GDP_percapita_sector_analysis.csv", header=T,sep=",") # loading data from gdp per capita
entry <- read.csv( "Entry_Costs_sector_analysis.csv", header=T,sep=",") # loading data from DB entry costs
entryL <- read.table("db_laborunits_sector_analysys.txt") # Doing Business in units of labor 
regcoeffs <- read.csv( "RegCoeffs_by_Country_WeightedReg_tr5.csv", header=T,sep=",") # loading data reg.coefficients distortions
pwt <- read.csv('pwt90_gdp_kstock.csv')   # penn world table
ctfp <- read.csv('ctfp_pwt.dat')          # TFP relative to US, from penn world tables
res_db <- read.table('res_db_weighted_reg_tr5.txt')                   # results from stationary equilibriums imputing entry barriers from Doing Business Indicators
res2_db <- read.table('resgains_db_weighted_reg_tr5.txt')             # results from stationary equilibriums imputing entry barriers from DBI: counterfactuals that remove one distortion, setting the other one to zero
LDage <-  read.table('LdAWeightedRegTr5.pc')                          # life-cycle dynamics distorted countries
agefless <- read.table('LdAFlessAvs10ind.pc')                         # life-cycle dynamics undistorted economy (US)
Mshare250census <- read.table("Mshare250censuses.txt")                # share of firms with 250 workers or more, from Census-based databases
Mshare250amadeus <-read.table("Mshare250amadeus.txt")                 # share of firms with 250 workers or more, from Amadeus-based databases

levels(avs$country)[levels(avs$country)=="ethiopia"] <- "Ethiopia"   # consistency upper case and lower case in country names
levels(avs$country)[levels(avs$country)=="ghana"] <- "Ghana"         # consistency upper case and lower case in country names
levels(avs$country)[levels(avs$country)=="salvador"] <- "Salvador"   # consistency upper case and lower case in country names
merged <- merge(avs, gdp, by.y='country') #,all=T)
merged <- merge(merged, regcoeffs, by.y='country')#,all=T)
merged <- merge(merged, entry, by.y=c('country','case'))#,all=T)
merged$regcoeff <- merged$regcoeff - 0.15  #controlling for USA misallocation
dist_db <- read.table('dist_weighted_reg_tr5_DB.txt')               # loading results from varios SS where TauE is fed with Doing Business Indicator's
merged_db <- merge(avs,dist_db, by=c('country','year'))
avsize_usa <- 118    # average size US manufacturing L>=10

#Parameters
#-----------
ne <- 31.65   #entry cost

# Labelling
#-------------
namesdata <- c('case', 'eta', 'TauE', 'Me', 'Y', 'TFP','TFPpm','AvZa', 'Lp'
               , 'Wage','EntryR','ExitR', 'M','TauEWage_Y','AvSize','AvSize10','slopeL','LI','AvRDint'
               , 'AvLe_LP','Lshare90','LshareTop','MshareTop','MshareBot')
colnames(data) <- namesdata

# Constructing Tables with Data to be Plotted in various figures
#----------------------------------------------------------------
TauE<- subset(data, TauEWage_Y!=0 & eta>0)       #varying entry tax, and varying misallocation
dist <- merge(gdp,TauE,by='case')
dist <- merge(dist,entryL,by=c('case','country'))
dist <- dist[order(dist$case),]   #sorting so that scatter plots dont' get messed up


#Frictionless Vector
Yfless <- data$Y[data$eta==0 & data$TauEWage_Y==0]
AvSize10fless <- data$AvSize10[data$eta==0 & data$TauEWage_Y==0]


#Subsets of data
fless <- subset(data, TauEWage_Y==0 & eta==0)     #frictionless
TauE<- subset(data, TauEWage_Y!=0 & eta>0)       #varying entry tax, and varying misallocation
TauEeta0 <- subset(data, TauEWage_Y!=0 & eta==0) #varying entry tax, misallocation=0
EtaTauE0 <- subset(data,TauEWage_Y==0 & eta>0)  #varying misallocation, TauE=0
TauE$Ynorm <- TauE$Y/Yfless
TauEeta0$Ynorm <- TauEeta0$Y/Yfless
TauE$AvSize10norm <- TauE$AvSize10/AvSize10fless
TauEeta0$AvSize10norm <- TauEeta0$AvSize10/AvSize10fless

logAvSize10normEta0 <- log(TauEeta0$AvSize10norm)

EtaTauE0$Ynorm <- EtaTauE0$Y/Yfless
EtaTauE0$AvSize10norm <- EtaTauE0$AvSize10/AvSize10fless
logAvSize10normTauE0 <- log(EtaTauE0$AvSize10norm)

Ygain_all <- (TauE$Ynorm^-1)
Mgain_all <- (TauE$M/fless$M)^-1

Ygain_eta <- TauEeta0$Y/TauE$Y       #gain from removing misallocatoin only, relative to distorted eq.
Ygain_tauE <- EtaTauE0$Y/TauE$Y      #gain from removing entry barrier only, relative to distorted eq.
Ygain_eta_tauE0 <- Yfless/EtaTauE0$Y  #gain from removing misallocation starting from eq. with misallocation only
Ygain_tauE_eta0 <- Yfless/TauEeta0$Y  #gain from removing entry barrier starting from entry barrier only

Mgain_eta <- TauEeta0$M/TauE$M
Mgain_tauE <- EtaTauE0$M/TauE$M
Mgain_eta_fless <- fless$M/EtaTauE0$M
Mgain_tauE_fless <- fless$M/TauEeta0$M

Res <- matrix(c(Ygain_all, Ygain_eta,Ygain_tauE,TauE$case),ncol=4)
colnames(Res) <- c("totalgain","gainmisalloc","gaineb","case")

Res <- merge(Res,gdp, by='case')
#************************************************************************************************
# FIGURE AVERAGE SIZE (controlling for US shares) and GDPpc
#  Appendix: COMPARING AV.SIZE CONTROLLING FOR PRODUCTION STRUCTURES, VS WITHOUT
#************************************************************************************************

# Figure 1: AvSize10 vs Log_GDPpc
#*******************************************************************
corr_size_dev <- cor(log(merged$gdppc),log(merged$AvSize10_Usa_weight)) # log-log correlation
corr_size_dev <- round(corr_size_dev,digits=2)


f1 <- ggplot(merged, aes(log(gdppc),AvSize10_Usa_weight)) +
  geom_point(color='darkgray',size=5)+
  xlab('log GDP per capita')+
  ylab('Average Size')+
  xlim(5.5,11)+
  ylim(30,200)+
  theme_light()+
  theme(axis.text.x=element_text(size=12),axis.text.y = element_text(size=12))+
  theme(axis.title=element_text(size=13))+
  theme(panel.grid.major=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('AvSize_GDPpc_10plus_USweighted.pdf',  width = 7,
       height = 5,
       units = c("in"))

fapp1 <- ggplot(avs, aes(AvSize10,AvSize10_Usa_weight)) +
  geom_point(color='darkgray',size=5)+
  xlab('Av Size country-specific firm shares')+
  ylab('Av Size U.S. firm shares')+
  xlim(10,250)+
  ylim(10,250)+
  theme_light()+
  theme(axis.text.x=element_text(size=12),axis.text.y = element_text(size=12))+
  theme(axis.title=element_text(size=13))+
  theme(panel.grid.major=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('AvSize_CountryShares_vs_USshares.pdf',  plot=fapp1, width = 7,
       height = 5,
       units = c("in"))

#************************************************************************************************
# FIGURES CHARACTERIZING SLOPE TFPR-TFPQ
#************************************************************************************************

f2a <- ggplot(merged, aes(log(gdppc),regcoeff)) +
  geom_point(color='darkgray',size=5)+
  xlab('log GDP per capita')+
  ylab('Slope TFPR-TFPQ')+
  xlim(5.5,11)+
  ylim(-0.1,0.5)+
  theme_light()+
  theme(axis.text.x=element_text(size=12),axis.text.y = element_text(size=12))+
  theme(axis.title=element_text(size=13))+
  theme(panel.grid.major=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)

f2b <- ggplot(merged, aes(AvSize10_Usa_weight,regcoeff)) +
  geom_point(color='darkgray',size=5)+
  xlab('Average Size')+
  ylab('Slope TFPR-TFPQ')+
  xlim(30,200)+
  ylim(-0.1,0.5)+
  theme_light()+
  theme(axis.text.x=element_text(size=12),axis.text.y = element_text(size=12))+
  theme(axis.title=element_text(size=13))+
  theme(panel.grid.major=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)

f2 <- grid.arrange(f2a, f2b, nrow=1)
ggsave('RegCoeff_vs_AvSize_vs_GDPpc.pdf', plot = f2,  width = 7,
       height = 5,
       units = c("in"))

#************************************************************************************************
# FIGURE: Average Size with Idiosyncratic Distortions Only
#************************************************************************************************

avsizemisaloc <- cbind(EtaTauE0$AvSize10norm, EtaTauE0$case)
colnames(avsizemisaloc) <- c('AvSize10norm', 'case')
merged <- merge(merged, avsizemisaloc, by='case')
merged$AvSize10_Usa_weight_norm <- merged$AvSize10_Usa_weight/avsize_usa

f3 <- ggplot(merged, aes(AvSize10_Usa_weight_norm,AvSize10norm)) +
  geom_point(color='darkgray',size=5)+
  xlab('Data/USA')+
  ylab('Misallocation/Undistorted')+
  xlim(0,1.5)+
  ylim(0,1.5)+
  theme_light()+
  theme(axis.text.x=element_text(size=12),axis.text.y = element_text(size=12))+
  theme(axis.title=element_text(size=13))+
  theme(panel.grid.major=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()

ggsave('AvSize_Misaloc_Only_weighted_reg_tr5.pdf', plot = f3,  width = 7,
       height = 5,
       units = c("in"))


#************************************************************************************************
# FIGURE: Model Based Entry barriers and GDP per capita
#************************************************************************************************

f4 <- ggplot(dist, aes(log(gdppc) ,log(1+TauE))) +
  geom_point(color='darkgray',size=5)+
  #xlim(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1)+
  #ylim(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1)+
  scale_x_continuous(name='log GDP per capita',limits=c(5,12), breaks=seq(5,12,by=2))+
  scale_y_continuous(name='log(1+TauE)',limits=c(min(log(1+dist$TauE))-0.5,max(log(1+dist$TauE))+0.5), breaks=seq(-1,2.5,by=0.5))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('TauEModel_vs_GDPpc.pdf', plot = f4,  width = 7,
       height = 5,
       units = c("in"))

#************************************************************************************************
# FIGURE: MODEL BASED VS DOING BUSINESS INDICATOR'S ENTRY BARRIERS:
# Figure 3: TauEmodel vs TauE data
# Figure 4: Average size ratio feeding TauE model vs TauE data: Equilibrium vs Data
#************************************************************************************************

f5 <- ggplot(dist, aes(log(1+TauE) ,log(1+dist$total_L/100*(1/ne)))) +
  geom_point(color='darkgray',size=5)+
  #xlim(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1)+
  #ylim(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1)+
  scale_x_continuous(name='log(1+TauE Model)',limits=c(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1), breaks=seq(-0.5,2.5,by=0.5))+
  scale_y_continuous(name='log(1+Doing Business)',limits=c(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1), breaks=seq(-0.5,2.5,by=0.5))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('TauEModel_vs_DB_logs.pdf', plot = f5,  width = 7,
       height = 5,
       units = c("in"))

f6<- ggplot(merged_db, aes(AvSize10_Usa_weight/avsize_usa ,AvSize10norm)) +
  geom_point(color='darkgray',size=5)+
  #xlim(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1)+
  #ylim(min(log(1+dist$TauE))-0.2,max(log(1+dist$TauE))+0.1)+
  scale_x_continuous(name='Data_country/Data USA',limits=c(0.2,1.4), breaks=seq(0.2,1.4,by=0.2))+
  scale_y_continuous(name='Distorted/Undistorted',limits=c(0.2,1.4), breaks=seq(0.2,1.4,by=0.2))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('AvSizeRatio_DB_Data.pdf', plot = f6,  width = 7,
       height = 5,
       units = c("in"))

#************************************************************************************************
# FIGURE 7: Life-Cycles
#************************************************************************************************
time <- 1:40
pdf('Life_Cycles_cases_vs_USA.pdf', width=7,height=5)
par(mfrow=c(1,1))
plot(time, agefless[time,5]/agefless[1,5],ylim=c(0.9,11),ylab='Employment relative to Age=1',xlab='Age',pch=21,bg="black",log='y')
points(time, LDage[time,6]/LDage[1,6],pch=21,bg="darkgray")
points(time, LDage[time,7]/LDage[1,7],pch=8,bg="black")
points(time, LDage[time,18]/LDage[1,18],pch=5,bg="darkgray")
legend(1,6,legend=c('USA   ','France   ','Ghana   ','India   '),pch=c(21,21,8,5),pt.bg=c('black','darkgray','black','darkgray'))
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
dev.off()

#************************************************************************************************
# FIGURE 8: Share of Firms with 250 workers or more
#************************************************************************************************
Mshare250model <- subset(dist, select=c('case','country','MshareTop'))   # selecting information on firm shares with L>=250 from the model's output file
Mshare250amacensus <- rbind(Mshare250amadeus, Mshare250census)           # binding together the amadeus and census based data of firm shares with L>=250
colnames(Mshare250amacensus) <- c('case','Mshare250')                    # labeling
Topgraph <- merge(Mshare250model, Mshare250amacensus,by='case')          # merging model and data
f8 <- ggplot(Topgraph, aes(MshareTop,Mshare250)) +
  ggtitle("Share of firms with more than 250 workers")+ 
  xlab("Model") + ylab("Data")+
  geom_point(color='darkgray',size=5)+
  xlim(0,0.3)+
  ylim(0,0.3)+
  theme_light()+
  theme(axis.text.x=element_text(size=12),axis.text.y = element_text(size=12))+
  theme(axis.title=element_text(size=13))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('MshareTop_Model_vs_Data_weighted_tr5.pdf', plot = f8,  width = 7,
       height = 5,
       units = c("in"))

#************************************************************************************************
# FIGURE: COUNTERFACTUAL TFP GAINS, VARIOUS EXERCISES
# Figure 9: TFP gains vs GDPpc
# Figure 4: Average size ratio feeding TauE model vs TauE data: Equilibrium vs Data
#************************************************************************************************

#TOTAL GDP GAINS vs CURRENT GDPpc
#-------------------------------------

f9<- ggplot(Res, aes(log(gdppc),totalgain)) +
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5.5,12), breaks=seq(5.5,12,by=1))+
  scale_y_continuous(name='TFP-efficient / TFP-distorted',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()
ggsave('TFPgains_vs_GDPpc.pdf', plot = f9,  width = 7,
       height = 5,
       units = c("in"))


# SHARE OF ACTUAL TFP GAINS ACCOUNTED FOR BY ALL DISTORTIONS
#---------------------------------------------------------------

res_tfp <- merge(Res,ctfp,by='country')

# Computing TFP gaps for which PWT has no data on ctfp variable
gha <- subset(pwt, country=='Ghana' & year==2014, select=c('country','year','rgdpna','rkna','emp'))
eth <- subset(pwt, country=='Ethiopia' & year==2014, select=c('country','year','rgdpna','rkna','emp'))
slv <- subset(pwt, country=='El Salvador'& year==2014, select=c('country','year','rgdpna','rkna','emp'))
usa <- subset(pwt, country=='United States'& year==2014, select=c('country','year','rgdpna','rkna','emp'))
bgdsh <- subset(pwt, country=='Bangladesh' & year==2014, select=c('country','year','rgdpna','rkna','emp'))
pktn <- subset(pwt, country=='Pakistan' & year==2014, select=c('country','year','rgdpna','rkna','emp')) 

gha$tfp <- gha$rgdpna/((gha$rkna^(1/3))*(gha$emp^(2/3)))
eth$tfp <- eth$rgdpna/((eth$rkna^(1/3))*(eth$emp^(2/3)))
slv$tfp <- slv$rgdpna/((slv$rkna^(1/3))*(slv$emp^(2/3)))
bgdsh$tfp <- bgdsh$rgdpna/((bgdsh$rkna^(1/3))*(bgdsh$emp^(2/3)))
pktn$tfp <- pktn$rgdpna/((pktn$rkna^(1/3))*(pktn$emp^(2/3)))
usa$tfp <- usa$rgdpna/((usa$rkna^(1/3))*(usa$emp^(2/3)))

gha$tfpgap <- gha$tfp/usa$tfp
eth$tfpgap <- eth$tfp/usa$tfp
slv$tfpgap <- slv$tfp/usa$tfp
bgdsh$tfpgap <- bgdsh$tfp/usa$tfp
pktn$tfpgap <- pktn$tfp/usa$tfp

res_tfp$ctfp[res_tfp$case==15] <- slv$tfpgap  #allocating computation of tfp gap to missing ctfp
res_tfp$ctfp[res_tfp$case==4] <- eth$tfpgap  #allocating computation of tfp gap to missing ctfp
res_tfp$ctfp[res_tfp$case==7] <- gha$tfpgap  #allocating computation of tfp gap to missing ctfp
res_tfp$ctfp[res_tfp$case==20] <- bgdsh$tfpgap  #allocating computation of tfp gap to missing ctfp
res_tfp$ctfp[res_tfp$case==21] <- pktn$tfpgap  #allocating computation of tfp gap to missing ctfp


f10<- ggplot(res_tfp, aes(ctfp^-1,totalgain/(ctfp^-1))) +
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='TFP(USA) / TFP(i)',limits=c(1,12), breaks=seq(1,12,by=1))+
  scale_y_continuous(name='Fraction of TFP gap closed',limits=c(0.05,1.1), breaks=seq(0.05,1.1,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
  
ggsave('Fraction_TFPgap_USA.pdf', plot = f10,  width = 7,
       height = 5,
       units = c("in"))

# TFP GAIN DECOMPOSITION: CONTRIBUTION OF EACH DISTORTION
#---------------------------------------------------------

f11a <- ggplot(Res, aes(log(gdppc),totalgain)) +
  ggtitle('Full Reform')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Relative to Initial SS',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f11a <- f11a+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f11b <- ggplot(Res, aes(log(gdppc),gainmisalloc)) +
  ggtitle('Misalloc. Only')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Relative to Initial SS',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f11b <- f11b+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f11c <- ggplot(Res, aes(log(gdppc),gaineb)) +
  ggtitle('Entry Barriers Only')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Relative to Initial SS',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f11c <- f11c+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f11<- grid.arrange(f11a, f11b, f11c, nrow=1)
ggsave('TFPgain_Decomp.pdf', plot = f11,  width = 7,
       height = 5,
       units = c("in"))

# INTERACTION: TOTAL GAINS VS SUM OF PARTIALS
#-----------------------------------------------
Res$gain_sum <- (Res$gaineb-1)+(Res$gainmisalloc-1)  #sum of partial gains

f12 <- ggplot(Res, aes(gain_sum,(totalgain-1))) +
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='Sum of Partials',limits=c(-0.05,0.6), breaks=seq(0,.6,by=.1))+
  scale_y_continuous(name='Total Gain',limits=c(-0.05,0.6), breaks=seq(0,.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)+
  geom_abline()

ggsave('TotalGains_vs_SumOfPartials.pdf', plot = f12,  width = 7,
       height = 5,
       units = c("in"))

# COMPARING TFP GAINS OF REMOVING ENTRY BARRIERS: MODEL VS DB INDICATOR'S
#----------------------------------------------------------------------------

dframe <- data.frame(Res,res_db)
f13a <- ggplot(dframe, aes(log(gdppc),gaineb-gaineb.1)) +
  ggtitle('Remove Entry Barrier with Misallocation')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Differencial TFP gain',limits=c(-0.05,0.3), breaks=seq(-0.05,.3,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f13a <- f13a+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

dframe$gaineb_eff <- Ygain_tauE_eta0-1
dframe$gaineb_eff_DB <- res2_db$gaineb_eff

f13b <- ggplot(dframe, aes(log(gdppc),gaineb_eff-gaineb_eff_DB)) +
  ggtitle('Remove Entry Barrier with Misallocation')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Differencial TFP gain',limits=c(-0.05,0.3), breaks=seq(-0.05,.3,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f13b <- f13b+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f13<- grid.arrange(f13a, f13b,nrow=1)
ggsave('TFPgainDifferential_ModelBarrier_vs_DB.pdf', plot = f13,  width = 9,
       height = 5,
       units = c("in"))

# FIGURE TFP GAINS DECOMPOSITION FOR DOING BUSINESS' ENTRY BARRIERS
#----------------------------------------------------------------------

f14a <- ggplot(res_db, aes(log(gdppc),totalgain)) +
  ggtitle('Full Reform')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Relative to Initial SS',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f14a <- f14a+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f14b <- ggplot(res_db, aes(log(gdppc),gainmisalloc)) +
  ggtitle('Misalloc. Only')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Relative to Initial SS',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f14b <- f14b+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f14c <- ggplot(res_db, aes(log(gdppc),gaineb)) +
  ggtitle('Entry Barriers Only')+
  geom_point(color='darkgray',size=5)+
  scale_x_continuous(name='log GDPpc',limits=c(5,12), breaks=seq(5,12,by=1))+
  scale_y_continuous(name='Relative to Initial SS',limits=c(0.9,1.6), breaks=seq(0.9,1.6,by=0.1))+
  theme_light()+
  theme(axis.text.x=element_text(size=11),axis.text.y = element_text(size=11))+
  theme(axis.title=element_text(size=12))+
  theme(panel.grid.major=element_line(linetype='dotted'), panel.grid.minor=element_line(linetype='dotted'))+
  geom_text_repel(aes(label=country),size=3.5)
f14c <- f14c+theme(plot.title=element_text(size=12, face="bold",hjust=0.5))

f14<- grid.arrange(f14a, f14b, f14c, nrow=1)
ggsave('TFPgain_Decomp_DB.pdf', plot = f14,  width = 7,
       height = 5,
       units = c("in"))
