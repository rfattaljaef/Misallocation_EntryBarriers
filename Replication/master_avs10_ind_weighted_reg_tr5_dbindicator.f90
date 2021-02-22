

!STEADY STATE CODE, ENDOGENOUS INNOVATION
! Master file solving cases for all countries at once

!avs10: here we calibrate entry barriers and measure TFP losses using the calibration to the US's average size L>=10.

!ind: indicates we are using the average size, conditional on 10+, fixing the distribution of firms across 2-digit ISIC-rev4 industries in the US

!weighted_reg stands for computing TFPR-TFPQ regression coefficient using weighted least squares with employment shares at the firm level as weights

!tr5 stands for estimating regcoeff TFPR-TFPQ trimming 5% of tails of TFPQ, TFPR

!dbindicator stands for feeding the doing business indocator's estimate of entry barrier. We use the value defined in units of labor, computed as : (DBI/100)*(1/labshare)*(1/ne)  where labshare is the labor share from PWT, and we divide by the entry cost to isolate the TauE component

program ssEndo
implicit none
real(8), parameter :: beta=1./1.05, delta=0.0250, s=0.25, nf=3.165, L=10.0000
real(8), parameter :: rho=4., ne=31.650000 !Labor supply, elast of subst, entry cost, transition length
integer, parameter :: N=180, nty=180, Nco=21, NC=4, NR=22, NuE=1, NIN=6 !length of vector of Productivities, number of cases
real(8), parameter :: AvSizeFless=116.64, AvSizeUSA10=118.77 ! average size in frictionless model(will vary depending we condition on L>10 or not)
real(8) :: AvSizeData(Nco),RegCoeff(Nco), AvSizeTarget  ! each country's av size, to be loaded from file
real(8) :: Measure(N), MeasureH(N), Mass(N), MeasureGuess(N), MassCases(N,NC*Nco) ! Productivity distributions
real(8) :: Tau(N), Corr(N), TauFunc(N), TauEvec(Nco), eta, TauE, TauEH, TauEL, TauMax, TauMin !Distortions
real(8) :: z(N) !vector of productivities
real(8) :: Pr(N,N),PrH(N,N) !Transition Matrix
real(8) :: ValueT(N) !Value Function at New Steady State
integer :: IndexZbar
real(8) :: Zbar !Cutoff    
real(8) :: ProfitDom, ProfitDomH, ProfitDomL
real(8) :: Wage, Y, Lp, Me, r !ProfitDomregates
real(8) :: Wpm, LPpm, YWr, LFixedCostProd, Za, ZaP
real(8) :: ValueF(N), ValueG(N), ValueZ(N), EV(N), ProfitD(N),ValCases(N,NC*Nco)
real(8) :: CumM(N), CumL(N), cdfM(N), cdfL(N), LpFirm(N), LpInd(N), EmpTot(N), EmpTotTop(N)
real(8) :: pCases(N,NC*Nco), RDintCases(N,NC*Nco), LpFirmCases(N,NC*Nco),EmpTotCases(N,NC*Nco), TransInputCases(NC,NIN)
real(8) :: slopeM, slopeL
real(8) :: ProfitDomZa, ProfitDomZaP  !Average Productivities (NOT per unit of entrant)
real(8) :: test, EntryTest, test3, test4, testM, testTE, test6, Tests(NC*Nco,5)
integer :: i,ii,j,jj,t,it, itm, iy, ic, itR, itMed, itE,ivf, Nfp, counter2 !indexes in loops
integer :: IndexZ0, IndexZmedian, IndexM, IndexL, IndexZ0med !39. This number depends on the support I am using
real(8) :: Lp10, LI10, M10, AvSize10, EntryR, ExitR, Lfc10
real(8) :: ProfitDomIntan, ProfitDomIntanSh, IntanMeSh, IntanPrSh
real(8) :: TFP, TFPpm, AvSize, M, ExitRInv, Results(NC*Nco,NR)
integer :: ResultsInt(3)
real(8) :: TFPcf, ZaCF, ZmedG, Zmed, AvTFP, ShareTax
real(8) :: VcompareL(2), VcompareH(2)

!Calibration Targets
real(8) :: AvLE, AvLP,Lshare90, LshareTop
integer :: Index90

!for interpolation of median
integer, parameter :: sr=100
real(8) :: ZI(sr), CDFint(sr), counter
integer :: IndexMint

!for interpolation of value function at cutoff
integer :: IzbarG
real(8) :: VVint(sr), ZG(sr), ED  ! ED stands for fraction of marginal type that stays in operation or exits (interpolation

real(8), parameter :: ga=4.5, ZZmin=1. , a=1-1./(rho-1.)	!parameters of pareto a is the DRS equivalente of rho=1. just to be able to use code from DRS model
real(8) :: ZZmax, cdfPar(nty), ZZ(nty), EZZ(nty), PP(nty) 
real(8), parameter :: npercmin=0.1, npercmax=0.9999, jump=0.25

!Innovation parameters and variable declaration
real(8), parameter :: phi=15., h=0.00056 ! phi=5., h=0.02102
real(8) :: p(N), c(N), LIpm, LI !innovation
real(8) :: Sales(N), RDint(N), AvRDint
!----------------------------------------------------------

open (1,file="avsize10_usa_weight.txt")   
read (1,*) AvSizeData(1:Nco)
close (1)

open (2,file="regcoeff_weighted_reg_tr5.txt")    !list of regcoeffs for countries , trimming tails at 5%
read (2,*) RegCoeff(1:Nco)
close (2)

open (2,file="db_indicator_laborunits.txt")    !doing business indicator's estiamte of TauE
read (2,*) TauEvec(1:Nco)
close (2)

! Adjusting or USA's regression coefficient. no country can by less misallocated than US!
do ic =1,Nco-1
if (regcoeff(ic) < 0.15) then
   regcoeff(ic) = 0.001
else
   regcoeff(ic) = regcoeff(ic)-0.15
endif
enddo



print *, RegCoeff

!***********************************************************!

! CONSTRUCTING SUPPORT (comparably with OC and DRS model)
!***********************************************************!


!ZZ is the support in the DRS, OC model
ZZ(1)=0.1 !(1./(1.-npercmin))**(1./ga)*ZZmin;

do i=2,nty
  ZZ(i)=ZZ(i-1)*(exp((1.-a)*jump)); ! constructed so that log[Z(i)^(1/(1-al-th))]-log[Z(i-1)^(1/(1-al-th))]=jump
enddo

do i=1,N
cdfPar(i)=max(1.-(ZZmin/ZZ(i))**ga,0.);
enddo

cdfPar=cdfPar/cdfPar(nty);
PP(1)=cdfPar(1);
do i=2,nty;
    PP(i)=cdfPar(i)-cdfPar(i-1);    ! given that DRS was calibrated to CES, then ZZ is equal to exp(Z)^(1/(rho-1))
                                    ! so PP still defines de distribution of entrants in the  CES model
enddo

EZZ=ZZ**(1./(1.-a))   ! ZZ^(1/(1-a)) corresponds to exp(Z) in CES model where exp(Z)^(1/(rho-1)) is TFP
Z=log(EZZ)

IndexZ0med=minloc(abs(cdfPar-0.5),1)

TauMax=0.8
TauMin=-0.8

!**************************************************
! STARTING EXPERIMENTS
!************************************************
counter2=0

do iy=1,Nco   !loop on country
AvSizeTarget=(AvSizeData(iy)/AvSizeUSA10)*AvSizeFless  !target for identification of entry tax
print *, 'AV SIZE TARGET', AvSizeTarget
do ic=1,NC    ! loop on cases of distortions

counter2=counter2+1
ValueG=exp(z)

if (abs(ic-1)< 0.001) then   !frictionless
   eta=0.
   Nfp=1
   TauE=0.
elseif (abs(ic-2)<0.001) then  ! benchmark with 2 distortions
  itE=1
  TauE=TauEvec(iy)
  eta=RegCoeff(iy)
  Nfp=40
elseif (abs(ic-3)<0.001) then
   TauE=TauEvec(iy)  !! keep the TauE we just estimated
   eta=0.     !! no misallocation
   Nfp=1
elseif (abs(ic-4)<0.001) then !!no entry costs
    TauE=0.
   eta=RegCoeff(iy)
   Nfp=1
endif
   
print *, 'CASE', TauE, eta, Nfp, ic

   IndexZmedian=IndexZ0med !Initial Guess for Median
   ZmedG=z(IndexZmedian)
   TauFunc=(exp(z))**((-eta)/(rho-1.)) 
   tau=1-TauFunc
   tau=min(tau, TauMax)
   tau=max(tau,TauMin)

do itMed=1,1 ! 20  iteration on median (shut down)

ProfitDomH=400.5
ProfitDomL=0.05
ProfitDom=0.5*(ProfitDomL+ProfitDomH)
do it=1,1000 !Bisection on ProfitDomregates

!*********************
!PROBLEM OF THE FIRMS
!*********************

r=1./beta-1.
do ivf=1,800 !VFI
 
 do jj=1,N

 if (abs(jj-1)<0.001) then  ! if we are at min Z
  p(jj)=(log(beta*(1-delta)*(ValueG(jj+1)-ValueG(jj)))-log(phi*h)-z(jj))/phi ! FOC wrt innovation
  if (p(jj)<0.) then 
   p(jj)=0.
  elseif (p(jj)>1.) then
   p(jj)=1.
  endif

  c(jj)=exp(z(jj))*h*exp(phi*p(jj))
  ProfitD(jj)=ProfitDom*exp(z(jj))*((1.-Tau(jj))**real(rho))-nf-c(jj)
  EV(jj)=(1-delta)*beta*(p(jj)*ValueG(jj+1)+(1-p(jj))*ValueG(jj)) 
  ValueZ(jj)=ProfitD(jj)+EV(jj)

 elseif (abs(jj-N)<0.001) then !if we are at max Z

  p(jj)=(log(beta*(1-delta)*(ValueG(jj)-ValueG(jj-1)))-log(phi*h)-z(jj))/phi ! FOC wrt innovation
  if (p(jj)<0.) then 
   p(jj)=0.
  elseif (p(jj)>1.) then
   p(jj)=1.
  endif

  c(jj)=exp(z(jj))*h*exp(phi*p(jj))
  ProfitD(jj)=ProfitDom*exp(z(jj))*((1.-Tau(jj))**real(rho))-nf-c(jj)
  EV(jj)=(1-delta)*beta*(p(jj)*ValueG(jj)+(1-p(jj))*ValueG(jj-1)) 
  ValueZ(jj)=ProfitD(jj)+EV(jj)
 else  ! if minZ < Z < maxZ
  p(jj)=(log(beta*(1-delta)*(ValueG(jj+1)-ValueG(jj-1)))-log(phi*h)-z(jj))/phi ! FOC wrt innovation

  if (p(jj)<0.) then 
   p(jj)=0.
  elseif (p(jj)>1.) then
   p(jj)=1.
  endif

  c(jj)=exp(z(jj))*h*exp(phi*p(jj))
  ProfitD(jj)=ProfitDom*exp(z(jj))*((1.-Tau(jj))**real(rho))-nf-c(jj)
  EV(jj)=(1-delta)*beta*(p(jj)*ValueG(jj+1)+(1-p(jj))*ValueG(jj-1))
  ValueZ(jj)=ProfitD(jj)+EV(jj)
 endif
           
 enddo  
  ValueF=max(ValueZ,0.)

  test=maxval(abs(ValueF(1:N)-ValueG(1:N)))
  
  if (test<0.0000000001) exit    
   ValueG=ValueF
enddo!VFI
  IndexZbar=minloc(ValueF(1:N),1, mask=ValueF(1:N)>0.0)  !At IndexZbar firm stays in operation
  Zbar=z(IndexZbar)  !becuase IndexZbarH<IndexZbar

EntryTest=ne*(1+TauE)-(1-delta)*beta*sum(ValueF*PP)



if (EntryTest .lt. -0.00000001) THEN
  ProfitDomH=ProfitDom
  ProfitDom=0.5*(ProfitDomL+ProfitDomH)
else if (EntryTest .gt. 0.00000001) THEN
  ProfitDomL=ProfitDom
  ProfitDom=0.5*(ProfitDomL+ProfitDomH)
end if

if (abs(EntryTest)<0.00000001) exit

enddo

print *, 'entry', EntryTest, ProfitDom, it

! Identifying Marginal Type
!------------------------------
! Interpolation nodes are: IndexZbar, IndexZbar-1
do j=1,sr
  ZG(j)=z(IndexZbar-1)+(j-1)*(z(IndexZbar)-z(IndexZbar-1))/sr
  VVint(j)=ValueZ(IndexZbar-1)+(ValueZ(IndexZbar)-ValueZ(IndexZbar-1))/(z(IndexZbar)-z(IndexZbar-1))*(ZG(j)-z(IndexZbar-1)) !Notice interpolation over ValueZ, not ValueF
enddo
IzbarG=minloc(VVint(1:sr),1, mask=VVint(1:sr)>0.0)
ED=(z(IndexZbar)-ZG(IzbarG))/(z(IndexZbar)-z(IndexZbar-1))  ! Fraction of agents at IndexZbar-1 that will stay in operation

print *, 'CASE', tauE, eta, itE

!**************************
!STATIONARY DISTRIBUTION
!**************************
MeasureGuess=0.
do it=1,2500 !@L
    do t=1,N
         if (abs(t-IndexZbar)<0.00001) then
             Measure(t)=(1-delta)*(1-p(t+1))*MeasureGuess(t+1)+(1-delta)*p(t-1)*MeasureGuess(t-1)+(1-delta)*PP(t)
         elseif (t<IndexZbar-1) then
             Measure(t)=0.     
         elseif (abs(t-(IndexZbar-1))<0.0001) then
             Measure(t)=(1-delta)*ED*(1-p(t+1))*MeasureGuess(t+1)+(1-delta)*ED*PP(t)      
         elseif (abs(t-N)<0.0001) then
             Measure(t)=(1-delta)*p(t)*MeasureGuess(t)+(1-delta)*p(t-1)*MeasureGuess(t-1)+(1-delta)*(PP(t));
         else
             Measure(t)=(1-delta)*(1-p(t+1))*MeasureGuess(t+1)+(1-delta)*p(t-1)*MeasureGuess(t-1)+(1-delta)*(PP(t));
         end if
     enddo
	 Measure(1:IndexZbar-1)=0.
     test4=maxval(abs(Measure-MeasureGuess))
     MeasureGuess=Measure;
     if (test4<0.0000000000000000001) exit
enddo !@L


!*******************************************************
!COMPUTING MEDIAN AND EMPLOYMENT-WEIGHTED MEDIAN
!*******************************************************
LpFirm=((exp(z)*(1-Tau)**rho)+c+nf)*Measure
CumM(1)=Measure(1)
CumL(1)=LpFirm(1)
do i=2,N
  CumM(i)=CumM(i-1)+Measure(i)
  CumL(i)=CumL(i-1)+LpFirm(i)  
enddo
cdfM(1:N)=CumM(1:N)/CumM(N)
cdfL(1:N)=CumL(1:N)/CumL(N)
IndexM=minloc(abs(cdfM(1:N)-0.5),1) !Median Productivity
IndexL=minloc(abs(cdfL(1:N)-0.5),1) !Median Employment
Index90=minloc(abs(cdfM(1:N)-0.9),1)


!------Interpolation of Median----------------------------------------
if (cdfM(IndexM)>0.5) then   ! identifying interpolation nodes
do j=1,sr
  ZI(j)=z(IndexM-1)+(j-1)*(z(IndexM)-z(IndexM-1))/sr
  CDFint(j)=cdfM(IndexM-1)+(cdfM(IndexM)-cdfM(IndexM-1))/(z(IndexM)-z(IndexM-1))*(ZI(j)-z(IndexM-1)) !Notice interpolation over ValueZ, not ValueF
enddo
IndexMint=minloc(CDFint(1:sr),1, mask=CDFint(1:sr)>0.0)
Zmed=ZI(IndexMint)

elseif (cdfM(IndexM)<0.5) then

do j=1,sr
  ZI(j)=z(IndexM)+(j-1)*(z(IndexM+1)-z(IndexM))/sr
  CDFint(j)=cdfM(IndexM)+(cdfM(IndexM+1)-cdfM(IndexM))/(z(IndexM+1)-z(IndexM))*(ZI(j)-z(IndexM)) !Notice interpolation over ValueZ, not ValueF
enddo
IndexMint=minloc(CDFint(1:sr),1, mask=CDFint(1:sr)>0.0)
Zmed=ZI(IndexMint)

else
Zmed=z(IndexM)

endif
!-----------------------------------------------------------------------
test6=abs(Zmed-ZmedG)

ShareTax=1-cdfM(IndexZmedian+1) 
print *, 'Median', test6, itMed, Zmed, ZmedG

if (test6.le.0.00001) exit
ZmedG=Zmed
IndexZmedian=IndexM 
TauFunc=(exp(z)/exp(z(IndexZmedian)))**((-eta)/(rho-1.))
tau=1-TauFunc
tau=min(tau, TauMax)
tau=max(tau, TauMin)

enddo ! iteration on median


print *, 'VFI', test, ivf, IndexZbar

!**************************
!RECOVERING AGGREGATES
!**************************

Za=sum(exp(z)*(1-Tau)**(rho)*Measure)       !aggregate productivities, per unit of entrants
ZaP=sum(exp(z)*(1-Tau)**(rho-1)*Measure)

LFixedCostProd=sum(nf*Measure);             !ProfitDomregate labor used for fixed costs, per unit of entrant
LIpm=sum(c*Measure)
YWr=ProfitDom*(rho**rho)/((rho-1.)**(rho-1.))  !gives Y/W^rho
LpInd=((rho-1.)/rho)**rho*YWr*(exp(z)*(1-Tau)**rho)
LpInd(1:IndexZbar-1)=0.
Wpm=(rho-1.)/rho*(ZaP**(1./(rho-1.)))          !wage, assuming Me=1;
LPpm=((rho-1.)/rho)**rho*YWr*Za;               !Labor Demand in Productio per unit of Entrant
Me=(L/(LPpm+ne+LFixedCostProd+LIpm))    
Lp=LPpm*Me;
Y=Lp*(Me**(1./(rho-1.)))*(ZaP**(rho/(rho-1.))/Za)
Wage=(rho-1.)/rho*(Me**(1./(rho-1.)))*(ZaP**(1./(rho-1.)))
LI=LIpm*Me

Mass(1:N)=Measure(1:N)*Me


M=sum(Mass(1:N))
TFPpm=(ZaP**(rho/(rho-1.)))/(Za)*(Me/M)**(1./(rho-1))
AvSize=(Lp+LI+LFixedCostProd*Me)/M                     !Av.Size using all workers
TFP=Y/Lp
ExitRInv=M/Me
ProfitDomZa=sum(exp(z)*(1-Tau)**(rho)*Mass)
ProfitDomZaP=sum(exp(z)*(1-Tau)**(rho-1)*Mass)

ZaCF=sum(exp(z(1:N))*Mass(1:N));
TFPcf=(ZaCF)**(1./(rho-1.));

ValueF(1:N)=ValueF(1:N)*Wage

slopeM=log((1-cdfM(IndexL+12))/(1-cdfM(IndexL)))/log(LpInd(IndexL+12)/LpInd(IndexL))
slopeL=log((1-cdfL(IndexL+12))/(1-cdfL(IndexL)))/log(LpInd(IndexL+12)/LpInd(IndexL))

! CONDITIONING ON L>0
EmpTot=LpInd(1:N)+c(1:N)+nf ! total employment at the firm level
EmpTot(1:IndexZbar-1)=0. 

!calculating vector of employment >=250
do i=1,N

if (EmpTot(i) >=250 ) then
EmpTotTop(i) = EmpTot(i)
else
EmpTotTop(i) = 0.
endif

enddo

Lp10=sum(LpInd(1:N)*Mass(1:N),mask=EmpTot(1:N).ge.10.)
LI10=sum(c(1:N)*Mass(1:N),mask=EmpTot(1:N).ge.10.)
Lfc10=nf*sum(Mass(1:N),mask=EmpTot(1:N).ge.10.)
M10=sum(Mass(1:N),mask=EmpTot(1:N).ge.10.)
AvSize10=sum(EmpTot(1:N)*Mass(1:N),mask=EmpTot(1:N).ge.10.)/M10 !  (Lp10+LI10+Lfc10)/M10

!**************************
!INNOVATION STATISTICS
!**************************
Sales=((rho-1.)/rho)**(rho-1.)*YWr*(exp(z)*(1-Tau)**rho)
RDint=(c/Sales)


AvLE=sum(LpInd(IndexZbar:N)*PP(IndexZbar:N))/sum(PP(IndexZbar:N))  !average employment entrant (dropping aggregates)
AvLP=sum(LpInd(IndexZbar:N)*measure(IndexZbar:N))/sum(measure(IndexZbar:N)) 
Lshare90=sum(EmpTot(Index90:N)*Mass(Index90:N))/sum(EmpTot*Mass)
LshareTop = sum(EmpTot(1:N)*Mass(1:N),mask=EmpTot(1:N).ge.250.)/sum(EmpTot*Mass)
AvRDint=sum(RDint*Mass)/M
EntryR=sum(PP(IndexZbar:N))*Me/M
ExitR=(delta*M+Mass(IndexZbar)*(1-p(IndexZbar)))/M

ResultsInt=(/IndexM, IndexL, IndexZbar/)  !results for integer variables
ProfitDomIntan=Me*ne*Wage+Wage*LI  !aggregate investment in intangibles
ProfitDomIntanSh=ProfitDomIntan/Y           ! investment in intangibles share of GDP
IntanMeSh=Me*ne*Wage/Y       ! intangible investment due to entry, share of GDP
IntanPrSh=Wage*LI/Y          ! intangible investment due to process innovation, share of GDP


counter=iy*1.
print *,'CASES', counter, iy, ic, eta, TauE
print *, 'Conver Cases', test, EntryTest, test4, testM, testTE
!print *, 'INPUT', Wage, Y, Me, LI, maxval(p)
!print *, 'ENTRANT-INCUMBENT', AvLE/AvLP, Lshare90, AvSize, AvSize10

Results(counter2,1:NR)=(/counter,eta,TauE,Me,Y,TFP,TFPpm,Za*Me/M,Lp,Wage,EntryR,ExitR,M,TauE*Wage*ne/(Y/L),AvSize,AvSize10,slopeL &
                      , LI, AvRDint,AvLE/AvLP, Lshare90, LshareTop/)
Tests(counter2,1:5)=(/test, EntryTest,test4, testM, testTE/)

ValCases(1:N,counter2)=ValueF(1:N)
pCases(1:N,counter2)=p(1:N)
RDintCases(1:N,counter2)=RDint(1:n)
LpFirmCases(1:N,counter2)=LpFirm(1:n)
EmpTotCases(1:N,counter2)=EmpTot(1:n)
MassCases(1:N,counter2)=Mass(1:N)
print *, 'Resutls', Results(counter2,1:2),LshareTop

open(11,file="ResMasterAvs10indWeightedRegTr5DB.dat")
do i=1,NC*Nco
write(11,*) Results(i,1:NR)
enddo
close(11)


enddo   ! loop on cases

enddo ! loop on countries


end program ssEndo    
