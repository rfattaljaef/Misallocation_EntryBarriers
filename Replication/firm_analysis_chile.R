# Computes US weighted average firm-size in Chile, and computes regression coefficient of logTFPR and logTFPQ following Hsieh and Klenow's methodology
# the firm-level data is not provided in the repository. Please readme file for instructions on how to download it from Chile's statistical agency
# Below, I shall assume that the firm-level data is stored in the file named Enia_2014.dta

rm(list = ls())  # clean workspace
library(dplyr)
library(haven)

data <- read_dta("Enia_2014.dta")
data=data.frame(data)

names(data)[names(data) == 'anio'] <- 'year'
names(data)[names(data) == 'va'] <- 'PY'
names(data)[names(data) == 'emptot'] <- 'L'
names(data)[names(data) == 'thhano'] <- 'LH'  #hours worked
names(data)[names(data) == 'rsc'] <- 'WLs'    #wagebill without contract
names(data)[names(data) == 'rcc'] <- 'WLc'    #wagebill with contract
names(data)[names(data) == 'salter'] <- 'Kter'  #land
names(data)[names(data) == 'saledi'] <- 'Kb' #end of period buildings
names(data)[names(data) == 'salmaq'] <- 'Kmq' #end of period machinery
names(data)[names(data) == 'salveh'] <- 'Kv'   #end of period vehicles
names(data)[names(data) == 'salmue'] <- 'Kf'   #end of period furniture
names(data)[names(data) == 'saloat'] <- 'Koa'   #end of period other tangible assets

data$WL <-  data$rempag  # wage bill
data$K <- data$Kb + data$Kmq + data$Kv + data$Kf + data$Koa + data$Kter   # capital stock

# ISIC rev 4, 4 digit
data$ciiu4 <- as.numeric(as.character(data$ciiu4))
data <- subset(data, ciiu4>=1000 & ciiu4<=3400)
data$rev4 <- data$ciiu4

# Labor Share in value added by 4 digit sector, ISIC rev 4, USA
data_share <- read_dta("US_ISIC4_rev4_stata12.dta")
data <- merge(data, data_share, by.y="rev4", all=TRUE)  #may be merging rows that correspond to rev4's in the US data outside manufacturing. Don't worry will be dropped later when we impose L>=0
data$alpha_s <- data$alpha_US_rev4    #relabeling capital share



#************************
# COMPUTING AVERAGE SIZE
#************************

data$isic2 <- substr(data$rev4,1,2) #turning 4 digit rev 4 ISIC code into a 2 digit one

digit2 <- subset(data, select=c('rev4','isic2','year','L') )  #constraining data
digit2 <- subset(digit2, L>=0) # dropping observations with missing employment data
digit2 <- digit2 %>% group_by(isic2) %>% mutate(
  firms=n(),
  firms10=sum(L>=10),  #number of firms in isic2 category conditional on 10+
  L10=sum(L[L>=10])  #employment in isic2 category conditional on 10+
)
digit2 <- digit2 %>% group_by(isic2) %>%
  summarize_at(c('L','L10','year','firms','firms10'),mean,na.rm=T)
digit2 <- digit2 %>% group_by(isic2) %>% mutate(
  avsize2=L, # Unlike L10, L is the firm level employment, not the unconditional sum of employment within a 2 digit industry
  avsize2_10=L10/firms10
)

usa2digit <- read.csv('Manufacturing_Data_2digit_ISICrev4.csv') # loading US average size by 2-digit industry according to ISIC rev 4 classification
merged <- merge(digit2, usa2digit, by.y.=isic2, all=T)
merged[is.na(merged)] <- 0.   #if there's no production in a given sector, replace NA with zero

merged <- transform( merged,
                     TotFirmsChi =sum(firms),
                     TotFirmsChi10=sum(firms10),
                     TotFirmsUsa=sum(firmsIsic2),
                     ShFirmsChi=firms/sum(firms),
                     ShFirmsChi10=firms10/sum(firms10),
                     ShFirmsUsa=firmsIsic2/sum(firmsIsic2)
)
merged$AvSize10=sum(merged$L10)/merged$TotFirmsChi10                 #average size conditional on L>=10 aggregating according to Chile's share of firms in each 2 digit industry
merged$AvSize10_alt=sum(merged$avsize2_10*merged$ShFirmsChi)         #average size conditional on L>=10 but aggregating according to Share of firms in each 2 digit industry for entire size distribution (comparability with US)
merged$AvSize10_Usa_weight <- sum(merged$avsize2_10*merged$ShFirmsUsa) # average size conditional on L>=10 aggregating accoring to US' share of firms across 2 digit industries.

short <- merged %>%  summarise(AvSize = mean(AvSize), AvSize10 = mean(AvSize10),AvSize10_alt = mean(AvSize10_alt), AvSizeUsa = mean(AvSizeUsa), AvSize_Usa_weight = mean(AvSize_Usa_weight), AvSize10_Usa_weight = mean(AvSize10_Usa_weight))
short$country='Chile'
short$year=2014


#****************************************************************************************************
# COMPUTATION OF ELASTICITY OF IDIOSYNCRATIC DISTORTIONS WITH RESPECT TO FIRM-LEVEL PRODUCTIVITY
#****************************************************************************************************

#HOUSEKEEPING
tail <- 1
w <- 1         # wage
R <- .1       # interest rate
sigma <- 3     # elasticity of substitution
data$alpha <- data$alpha_s    # share of payments to capital in value added
data$k <- data$K              # relabelling capital stock
data$z <- data$WL #data$WL    # relabeling wage bill
data <- subset(data, PY>0 & K>0 & z >0  & L>=10)  # restricting to non missing values of value added, capital stock, and wage bill; and conditioning on L>=10



data$isic_4 <- data$rev4 #relabelling


data <- data %>%
  group_by(isic_4) %>%
  mutate(avsize_bysec=mean(L)) %>%   #average by isic_4 sector
  group_by(isic_4) %>%
  mutate(va_4d=sum(PY))
  
data <- data %>% group_by(year) %>% mutate(va_ag=sum(PY))
data$vashare_4d <- data$va_4d/data$va_ag
data$weights_PY <- data$PY/data$va_4d
theta <- data$vashare_4d


#----------------------------------------------------
#  Analysis of Misallocation: Hsieh and Klenow 2009
#---------------------------------------------------

data <- data %>% mutate(
  tau_k1 = 1- (alpha/(1-alpha))*z /(R*k),
  tau_y1= 1+ (sigma/(sigma-1))*z/((1-alpha)*PY) ,     #tau_y only
  tau_k = (alpha/(1-alpha))*z /(R*k) ,               # 1+tau_k
  tau_y = (sigma/(sigma-1))*z /((1-alpha)*PY)  ,       #1-tau_y
  kappa = 1,
  
  TFPQ = kappa*PY^(sigma/(sigma-1))/(k^(alpha) * z ^(1-alpha)),
  TFPR = (sigma/(sigma-1) )* (1/(1-alpha))^(1-alpha)*tau_k^(alpha)*(R/alpha)^(alpha) / tau_y # footnote 10
)

data <- data %>% group_by(isic_4) %>% mutate(
  Abar = sum(TFPQ^(sigma-1))^(1/(sigma-1)),
  temp1 =  sum(tau_y* weights_PY/tau_k) ,
  temp2 = sum( weights_PY*tau_y),
  TFPRbar = (sigma/(sigma-1) )* (R/(alpha*temp1))^(alpha)*(1/((1-alpha)*temp2))^(1-alpha), # footnote 11
  M=n()
)

data <- data %>% mutate(
  log_TFPQ = log(TFPQ*M^(1/(sigma-1))/Abar) ,
  log_TFPR=  log(TFPR/TFPRbar)
)

# Trimming Outliers at 5%
#--------------------------
data$p1_TFPR=quantile(data$log_TFPR, probs=0.05, na.rm=TRUE)
data$p99_TFPR=quantile(data$log_TFPR, probs=0.95, na.rm=TRUE)


data$p1_TFPQ=quantile(data$log_TFPQ, probs=0.05, na.rm=TRUE)
data$p99_TFPQ=quantile(data$log_TFPQ, probs=0.95, na.rm=TRUE)
data <- subset(data, log_TFPQ>+p1_TFPQ & log_TFPQ<=p99_TFPQ & log_TFPR>=p1_TFPR & log_TFPR<=p99_TFPR)

data <- data %>% group_by(year) %>% mutate(Nobs=n())


# Re-calculation Post-Trimming
#-------------------------------------
data <- data %>%
  group_by(isic_4) %>%
  mutate(avsize_bysec=mean(L)) %>%   #average by isic_4 sector
  group_by(isic_4) %>%
  mutate(va_4d=sum(PY))

data <- data %>% group_by(year) %>% mutate(va_ag=sum(PY),AggL = sum(L))
data$vashare_4d <- data$va_4d/data$va_ag
data$weights_PY <- data$PY/data$va_4d
data$weights_regr <- data$L/data$AggL
theta <- data$vashare_4d

data <- data %>% mutate(
  tau_k1 = 1- (alpha/(1-alpha))*z /(R*k),
  tau_y1= 1+ (sigma/(sigma-1))*z/((1-alpha)*PY) ,     #tau_y only
  tau_k = (alpha/(1-alpha))*z /(R*k) ,               # 1+tau_k
  tau_y = (sigma/(sigma-1))*z /((1-alpha)*PY)  ,       #1-tau_y
  kappa = 1,
  
  TFPQ = kappa*PY^(sigma/(sigma-1))/(k^(alpha) * z ^(1-alpha)),
  TFPR = (sigma/(sigma-1) )* (1/(1-alpha))^(1-alpha)*tau_k^(alpha)*(R/alpha)^(alpha) / tau_y # footnote 10
)

data <- data %>% group_by(year, isic_4) %>% mutate(
  Abar = sum(TFPQ^(sigma-1))^(1/(sigma-1)),
  temp1 =  sum(tau_y* weights_PY/tau_k) ,
  temp2 = sum( weights_PY*tau_y),
  TFPRbar = (sigma/(sigma-1) )* (R/(alpha*temp1))^(alpha)*(1/((1-alpha)*temp2))^(1-alpha), # footnote 11
  M=n()
)

data <- data %>% mutate(
  log_TFPQ = log(TFPQ*M^(1/(sigma-1))/Abar) ,
  log_TFPR=  log(TFPR/TFPRbar)
)

data <- data %>% group_by(isic_4) %>% mutate(
sd_logTFPR=sd(log_TFPR),
sd_logTFPQ=sd(log_TFPQ),
Nind=n()
)


# Weighted Least Squares Regression of Log(TFPR)  vs log(TFPQ)
#-----------------------------------------------------------------

regcoeff <- 0
year <- max(data$year)
reg <- lm(log_TFPR ~log_TFPQ, data=data,weights=weights_regr) 
regcoeff <- coefficients(reg)[2]                              #selecting slope coefficient
data_regcoef <- data.frame(regcoeff,year)
names(data_regcoef)[names(data_regcoef) == 'yyear'] <- 'year'  #relabelling

