
! Age Distribution and Life Cycle Dynamics, distorted SS for all countries under estimated TauE and regcoeffs



program age
implicit none

real(8), parameter :: beta=1./1.05, delta=0.0250, s=0.25, nf=3.165, L=10.0000
real(8), parameter :: rho=4., ne=31.650000 !Labor supply, elast of subst, entry cost, transition length
integer, parameter :: nty=180, N=180, AA=1000 ! number of types and processors (!), Max Age
integer, parameter :: Nco=19
real(8), parameter :: ga=4.5, ZZmin=1. , al=1-1./(rho-1.)	!parameters of pareto a is the DRS equivalente of rho=1. just to be able to use code from DRS model
real(8) :: ZZmax, cdfPar(nty), ZZ(nty), EZZ(nty), PP(nty), Z(nty) 
real(8) :: qcase(nty,Nco), PRcase(nty,Nco), FLDcase(nty,Nco), LDcase(nty,Nco)
real(8) :: AvLdAvec(AA,Nco), NAvec(AA,Nco), AvZAvec(AA,Nco), AvTauAvec(AA,Nco)
real(8), parameter :: npercmin=0.1, npercmax=0.9999, jump=0.25
integer :: OC(nty)  ! 0-1 indicator variable if producer or out of the market
integer :: IndexZmedian
real(8) :: TauFunc(nty), tau(nty), RegCoeff(Nco), eta
real(8) :: PR(nty), PRA(nty,AA), q(nty), AvZE
integer :: IndexZbar, IndexAvZE, IndexME, Index75E, Index90E
real(8) :: FLD(nty),LD(nty), Y(nty), FLDA(AA), LDA(AA), NA(AA), AvLdA(AA), AvZA(AA), ZZA(AA), TauA(AA), AvTauA(AA) ! life cycle
real(8) :: FLDAE(AA), NAE(AA), AvLdAE(AA), AvZAE(AA), ZZAE(AA), PRAE(nty,AA), AvLdAtot(AA)  ! life-cycle of cohort with expected productivity only
real(8) :: FLDAME(AA), NAEM(AA), AvLdAME(AA), AvZAME(AA), ZZAME(AA), PRAME(nty,AA)   ! life-cycle of cohort with expected productivity only
real(8) :: FLDA75E(AA), NA75E(AA), AvLdA75E(AA), AvZA75E(AA), ZZA75E(AA), PRA75E(nty,AA)   ! life-cycle of cohort with expected productivity only
real(8) :: FLDA90E(AA), NA90E(AA), AvLdA90E(AA), AvZA90E(AA), ZZA90E(AA), PRA90E(nty,AA)   ! life-cycle of cohort with expected productivity only
real(8) :: pdfE(nty), cdfE(nty), ZZI
integer :: a,i, ic

!ZZ is the support in the DRS, OC model
ZZ(1)=0.1 !(1./(1.-npercmin))**(1./ga)*ZZmin;

do i=2,nty
  ZZ(i)=ZZ(i-1)*(exp((1.-al)*jump)); ! constructed so that log[Z(i)^(1/(1-al-th))]-log[Z(i-1)^(1/(1-al-th))]=jump
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

EZZ=ZZ**(1./(1.-al))   ! ZZ^(1/(1-a)) corresponds to exp(Z) in CES model where exp(Z)^(1/(rho-1)) is TFP
Z=log(EZZ)

open (2,file="regcoeff_weighted_reg_tr5.txt")
read (2,*) RegCoeff(1:Nco)
close (2)

open (4,file="innovdist_weighted_reg_tr5.txt")
do i=1,nty
read (4,*) qcase(i,1:Nco)
enddo
 close (4)

open (5,file="massdist_weighted_reg_tr5.txt")
do i=1,nty
read (5,*) PRcase(i,1:Nco)
enddo
 close (5)

open (7,file="lpdist_weighted_reg_tr5.txt")
do i=1,nty
read (7,*) FLDcase(i,1:Nco)
enddo
 close (7)

 open (8,file="emptotdist_weighted_reg_tr5.txt") 
do i=1,nty
 read (8,*) LDcase(i,1:Nco)
enddo
 close (8)
 
do ic=1,Nco ! loop on countries

PR = PRcase(1:nty, ic)
FLD = FLDcase(1:nty,ic)
LD = LDcase(1:nty,ic)
q = qcase(1:nty, ic)

! PARAMETRIZATION OF DISTORTIONS
eta = RegCoeff(ic)  ! chile's regression coefficient
TauFunc=(exp(z))**((-eta)/(rho-1.)) !(exp(z)/exp(z(IndexZmedian)))**((-eta)/(rho-1.))
   tau=1-TauFunc
   tau=min(tau, .8)
   tau=max(tau,-.8)


! constructing indicator variable of production
do i=1,nty
  if( PR(i) >0 ) then
   OC(i) = 0  ! active producer
  else
   OC(i) = 1  ! out of market
  endif
enddo

PRA(1:nty,1)=PP(1:nty)
FLDA(1)=sum(FLD(1:nty)*PRA(1:nty,1)*(1.-OC(1:nty)))
LDA(1)=sum(LD(1:nty)*PRA(1:nty,1)*(1.-OC(1:nty)))
!FLDA(1)=sum((exp(z)*(TauFunc**rho))*PRA(1:nty,1)*(1.-OC(1:nty)))
ZZA(1)=sum(exp(z)*PRA(1:nty,1)*(1.-OC(1:nty)))
TauA(1)=sum(TauFunc**rho*PRA(1:nty,1)*(1.-OC(1:nty)))
NA(1)=sum(PRA(1:nty,1)*(1.-OC(1:nty)))
AvLdA(1)=FLDA(1)/NA(1)
AvLdAtot(1)=LDA(1)/NA(1)
AvZA(1)=ZZA(1)/NA(1)
AvTauA(1)=TauA(1)/NA(1)

do a=2,AA

  do i=1,nty

    if (abs(i-1)<0.001) then  ! boundary condition 1
     PRA(i,a)=PRA(i+1,a-1)*(1-q(i+1))*(1.-OC(i+1))+PRA(i,a-1)*(1-q(i))*(1.-OC(i))
    elseif (abs(i-nty)<0.001) then  ! boundary condition 2
     PRA(i,a)=PRA(i-1,a-1)*q(i-1)*(1.-OC(i-1))+PRA(i,a-1)*q(i)*(1.-OC(i))
    else
    PRA(i,a)=PRA(i+1,a-1)*(1-q(i+1))*(1.-OC(i+1))+PRA(i-1,a-1)*q(i-1)*(1.-OC(i-1))
    endif

   FLDA(a)=sum(FLD(1:nty)*PRA(1:nty,a)*(1.-OC(1:nty)))
   LDA(a)=sum(LD(1:nty)*PRA(1:nty,a)*(1.-OC(1:nty)))
   !FLDA(a)=sum((exp(z)*(TauFunc**rho))*PRA(1:nty,a)*(1.-OC(1:nty)))
   ZZA(a)=sum((exp(z))*PRA(1:nty,a)*(1.-OC(1:nty)))
   TauA(a)=sum((TauFunc**rho)*PRA(1:nty,a)*(1.-OC(1:nty)))
   NA(a)=sum(PRA(1:nty,a)*(1.-OC(1:nty)))
   AvLdA(a)=FLDA(a)/NA(a)
   AvLdAtot(a)=LDA(a)/NA(a)
   AvZA(a)=ZZA(a)/NA(a)
   AvTauA(a)=TauA(a)/NA(a)

 enddo

enddo

AvLdAvec(1:AA,ic) = AvLdAtot
NAvec(1:AA,ic) = NA
AvZAvec(1:AA,ic)=AvZA
AvTauAvec(1:AA,ic) = AvTauA

enddo ! loop countries


open (11,file="LdAWeightedRegTr5.pc")
do i=1,AA
write(11,*) AvLdAvec(i,1:Nco)
enddo
close (11)

open (12,file="NAWeightedRegTr5.pc")
do i=1,AA
write(12,*) NAvec(i,1:Nco)
enddo
close (12)

open (13,file="AvZAWeightedRegTr5.pc")
do i=1,AA
write(13,*) AvZAvec(i,1:Nco)
enddo
close (13)

open (13,file="AvTauAWeightedRegTr5.pc")
do i=1,AA
write(13,*) AvTauAvec(i,1:Nco)
enddo
close (13)

end program age

