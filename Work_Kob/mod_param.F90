
!!!=== Copyright (c) 2012-2020 Takashi NAKAMURA  =====

!--------------------------------------------------------------------------------
!
!  Parameter module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

MODULE mod_param
  
  implicit none
!
!-----------------------------------------------------------------------
!  Grid nesting parameters.
!-----------------------------------------------------------------------

  real(8) :: time    !(day)
  
  real(8), allocatable :: tide_data(:)
  real(8) :: tide
  real(8) :: dt_tide
  integer :: nm_tide
  
  real(8), allocatable :: Hs_data(:)
  real(8), allocatable :: Tp_data(:)
  real(8) :: Hs
  real(8) :: Tp
  real(8) :: dt_wave
  integer :: nm_wave
  
  real(8), allocatable :: Pair_data (:)
  real(8), allocatable :: Uwind_data(:)
  real(8), allocatable :: Vwind_data(:)
  real(8), allocatable :: Tair_data (:)
  real(8), allocatable :: Qair_data (:)
  real(8), allocatable :: rain_data (:)
  real(8), allocatable :: swrad_data(:)
  real(8), allocatable :: cloud_data(:)
  real(8) :: Pair
  real(8) :: Uwind
  real(8) :: Vwind
  real(8) :: Tair
  real(8) :: Qair
  real(8) :: rain
  real(8) :: swrad
  real(8) :: cloud
  real(8) :: dt_air
  integer :: nm_air
  
  real(8), allocatable :: dlwrad_data(:)
  real(8) :: dlwrad
  real(8) :: dt_dlwrad
  integer :: nm_dlwrad
  
  real(8), allocatable :: offtemp_data(:)
  real(8) :: offtemp
  real(8) :: dt_offtemp
  integer :: nm_offtemp
  
  real(8), allocatable :: PFD_time(:)
  real(8), allocatable :: PFD_data(:)
  integer :: N_PFD
  integer :: i_PFD = 1
  real(8), allocatable :: WQ_time(:)
  real(8), allocatable :: TA_data(:)
  real(8), allocatable :: DIC_data(:)
  real(8), allocatable :: DO_data(:)
  integer :: N_WQ
  integer :: i_WQ = 1


  real(8) :: PFDsurf    
  real(8) :: tau        
  real(8) :: pCO2air    
  real(8) :: U10        
  
  real(8) :: sspH      
  real(8) :: ssfCO2    
  real(8) :: ssWarg    
  real(8) :: ssCO2flux 
  real(8) :: ssO2flux  
  real(8) :: PFDbott   

  real(8) :: fvol_cre      ! volume flux through the reef crest (m3 m-2 s-1)
  real(8) :: fvol_cha      ! volume flux through the channel(m3 m-2 s-1)
  real(8) :: fvol_pre      ! Precipitation volume flux (m s-1)
!  real(8) :: dw_lwradi     ! Downward longwave radiation (W m-2)
  real(8) :: ereef         ! sea surface elevation on the reef flat (m)
  real(8), parameter :: z_crest = -0.15d0 ! reef crest position (m)
  real(8), parameter :: kvol_cre   = 5.0d-2 ! reef crest conductivity
  real(8), parameter :: kvol_cha = 1.0d-1 ! Channel conductivety
  real(8) :: wave_setup   ! wave setup (m)

  integer :: iTemp                  ! Temperature
  integer :: iSalt                  ! Salinity
  integer :: iSedi                  ! suspended sediment

#if defined ORGANIC_MATTER
  integer, parameter :: N_phyt = 3  ! Number of functional groups of phytoplankton
                                    !  1: dinoflagellate
                                    !  2: diatom
                                    !  3: coccolithophorids
  integer, parameter :: N_zoop = 1  ! Number of functional groups of zooplankton
                                    !  1: General Zooplankton 
  integer, parameter :: N_dom  = 2  ! Number of functional groups of Dissolved organoic matter
                                    !  1: Labile DOM
                                    !  2: Refractory DOM
  integer, parameter :: N_pom  = 2  ! Number of functional groups of Particulate organoic matter
                                    !  1: Detritus
                                    !  2: Coarse POM (leaf litter, etc.)
  integer, parameter :: N_pim  = 1  ! Number of functional groups of Particulate inorganoic matter
                                    !  1: coccolith (CaCO3)
#endif
  integer :: idCrl1                 ! Coral coverage
  integer :: idCrl2                 ! Coral2 coverage
  integer :: idSgrs                 ! Seagrass coverage
  integer :: idAlga                 ! Algal coverage
  integer :: idSand                 ! Sand coverage
  
  integer :: idd13C                 ! d13C of DIC

  integer, allocatable :: idbio(:)  ! Biological tracers
  integer :: iTIC_                  ! Total inorganic carbon
  integer :: iTAlk                  ! Total alkalinity
  integer :: iOxyg                  ! Dissolved oxygen concentration
#if defined ORGANIC_MATTER
  integer :: iDOC(N_dom)            ! Dissolved organic C-concentration
  integer :: iPOC(N_pom)            ! Particulate organic C-concentration
  integer :: iPhyt(N_phyt)          ! Phytoplankton1 density
  integer :: iZoop(N_zoop)          ! Zooplankton density
  integer :: iPIC(N_pim)            ! Particulate inorganic C-concentration
# if defined FOODWEB_TESTMODE
  integer :: iDeadPhyt(N_phyt)      ! Dead Phytoplankton density
# endif
#endif
#if defined CARBON_ISOTOPE
  integer :: iT13C                  ! Corbon 13 of total inorganic carbon
# if defined ORGANIC_MATTER
  integer :: iDO13C(N_dom)          ! Dissolved organic 13C-concentration
  integer :: iPO13C(N_pom)          ! Particulate organic 13C-concentration
  integer :: iPhyt13C(N_phyt)       ! Phytoplankton1 13C-concentration
  integer :: iZoop13C(N_zoop)       ! Zooplankton 13C-concentration
  integer :: iPI13C(N_pim)          ! Particulate inorganic 13C-concentration
# endif
#endif
#if defined NUTRIENTS
  integer :: iNO3_                  ! Nitrate concentration
!  integer :: iNO2_                  ! Nitrite concentration
  integer :: iNH4_                  ! Ammonium concentration
  integer :: iPO4_                  ! Ammonium concentration
# if defined ORGANIC_MATTER
  integer :: iDON(N_dom)            ! Dissolved organic N-concentration
  integer :: iPON(N_pom)            ! Particulate organic N-concentration
  integer :: iDOP(N_dom)            ! Dissolved organic P-concentration
  integer :: iPOP(N_pom)            ! Particulate organic P-concentration
# endif
# if defined NITROGEN_ISOTOPE
  integer :: i15NO3                  ! Nitrogen isotope concentration in Nitrate 
!  integer :: i15NO2                  ! Nitrogen isotope concentration in Nitrite
  integer :: i15NH4                  ! Nitrogen isotope concentration in Ammonium
#  if defined ORGANIC_MATTER
  integer :: iDO15N(N_dom)          ! Dissolved organic 15N-concentration
  integer :: iPO15N(N_pom)          ! Particulate organic 15N-concentration
  integer :: iPhyt15N(N_phyt)       ! Phytoplankton 15N-concentration
  integer :: iZoop15N(N_zoop)       ! Zooplankton 15N-concentration
#  endif
# endif
#endif
#if defined COT_STARFISH
  integer :: iCOTe                  ! Eggs of crown-of-thorns starfish
  integer :: iCOTl                  ! Larvae of crown-of-thorns starfish
#endif

  real(8), allocatable :: dz(:,:,:)    
  real(8), allocatable :: C(:,:,:,:,:)
  real(8), allocatable :: dC_dt(:,:,:,:)
  real(8), allocatable :: Co(:,:,:)
  
  real(8), allocatable :: p_coral(:,:,:)
  real(8), allocatable :: p_sgrass(:,:)
  real(8), allocatable :: p_algae(:,:)
  real(8), allocatable :: p_sand(:,:)

contains

!!! **********************************************************************
!!!  Set initial conditions
!!! **********************************************************************

  subroutine initialize_params(LBi, UBi, LBj, UBj, N, Nid)

    use mod_geochem
    
    implicit none
! input parameters
    integer, intent(in) :: LBi, UBi, LBj, UBj, N
    integer, intent(out) :: Nid
    real(8)  R13C
    integer i,j,k
!  Local variable declarations
!
    integer :: ic

    real(8) :: d0, T0, S0
    real(8) :: TAlk0, TIC_0, Oxyg0, Sed1_0
    real(8) :: DOC1_0, DOC2_0, POC1_0, POC2_0, Phyt1_0, Phyt2_0, Phyt3_0, Zoop1_0, PIC1_0
    real(8) :: d13C_TIC0, d13C_DOC1_0, d13C_DOC2_0, d13C_POC1_0, d13C_POC2_0
    real(8) :: d13C_Phyt1_0, d13C_Phyt2_0, d13C_Phyt3_0, d13C_Zoop1_0, d13C_PIC1_0
    real(8) :: NO3_0, NH4_0, PO4_0
    real(8) :: DON1_0, DON2_0, PON1_0, PON2_0, DOP1_0, DOP2_0, POP1_0, POP2_0
    real(8) :: d15N_NO3_0, d15N_NH4_0, d15N_DON1_0, d15N_DON2_0, d15N_PON1_0, d15N_PON2_0
    real(8) :: d15N_Phyt1_0, d15N_Phyt2_0, d15N_Phyt3_0, d15N_Zoop1_0
    real(8) :: p_coral1_0, p_coral2_0, p_algae_0, p_sgrass_0, p_sand_0

    namelist/initial/d0, T0, S0
    namelist/initial/TAlk0, TIC_0, Oxyg0, Sed1_0
    namelist/initial/DOC1_0, DOC2_0, POC1_0, POC2_0, Phyt1_0, Phyt2_0, Phyt3_0, Zoop1_0, PIC1_0
    namelist/initial/d13C_TIC0, d13C_DOC1_0, d13C_DOC2_0, d13C_POC1_0, d13C_POC2_0
    namelist/initial/d13C_Phyt1_0, d13C_Phyt2_0, d13C_Phyt3_0, d13C_Zoop1_0, d13C_PIC1_0
    namelist/initial/NO3_0, NH4_0, PO4_0
    namelist/initial/DON1_0, DON2_0, PON1_0, PON2_0, DOP1_0, DOP2_0, POP1_0, POP2_0
    namelist/initial/d15N_NO3_0, d15N_NH4_0, d15N_DON1_0, d15N_DON2_0, d15N_PON1_0, d15N_PON2_0
    namelist/initial/d15N_Phyt1_0, d15N_Phyt2_0, d15N_Phyt3_0, d15N_Zoop1_0
    namelist/initial/p_coral1_0, p_coral2_0, p_algae_0, p_sgrass_0, p_sand_0
  
    read (*, nml=initial)
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
    ic=0

    i=1
    iTemp=ic+i
    i=i+1
    iSalt=ic+i

    i=i+1
    iSedi=ic+i

    i=i+1
    iTIC_=ic+i
    i=i+1
    iTAlk=ic+i
    i=i+1
    iOxyg=ic+i  !  4
#if defined ORGANIC_MATTER
    DO j=1,N_dom
      i=i+1
      iDOC(j)=ic+i
    END DO
    DO j=1,N_pom
      i=i+1
      iPOC(j)=ic+i
    END DO
    DO j=1,N_phyt
      i=i+1
      iPhyt(j)=ic+i
    END DO
    DO j=1,N_zoop
      i=i+1
      iZoop(j)=ic+i
    END DO
    DO j=1,N_pim
      i=i+1
      iPIC(j)=ic+i
    END DO
# if defined FOODWEB_TESTMODE
    DO j=1,N_phyt
      i=i+1
      iDeadPhyt(j)=ic+i
    END DO
# endif
#endif
#if defined CARBON_ISOTOPE
    i=i+1
    iT13C=ic+i  ! +1
# if defined ORGANIC_MATTER
    DO j=1,N_dom
      i=i+1
      iDO13C(j)=ic+i
    END DO
    DO j=1,N_pom
      i=i+1
      iPO13C(j)=ic+i
    END DO
    DO j=1,N_phyt
      i=i+1
      iPhyt13C(j)=ic+i
    END DO
    DO j=1,N_zoop
      i=i+1
      iZoop13C(j)=ic+i
    END DO
    DO j=1,N_pim
      i=i+1
      iPI13C(j)=ic+i
    END DO
# endif
#endif
#if defined NUTRIENTS
    i=i+1
    iNO3_=ic+i
!    i=i+1
!    iNO2_=ic+i
    i=i+1
    iNH4_=ic+i
    i=i+1
    iPO4_=ic+i
# if defined ORGANIC_MATTER
    DO j=1,N_dom
      i=i+1
      iDON(j)=ic+i
    END DO
    DO j=1,N_pom
      i=i+1
      iPON(j)=ic+i
    END DO
    DO j=1,N_dom
      i=i+1
      iDOP(j)=ic+i
    END DO
    DO j=1,N_pom
      i=i+1
      iPOP(j)=ic+i
    END DO
# endif
# if defined NITROGEN_ISOTOPE
    i=i+1
    i15NO3=ic+i  ! +1
!    i=i+1
!    i15NO2=ic+i  ! +1
    i=i+1
    i15NH4=ic+i  ! +1
#  if defined ORGANIC_MATTER
    DO j=1,N_dom
      i=i+1
      iDO15N(j)=ic+i
    END DO
    DO j=1,N_pom
      i=i+1
      iPO15N(j)=ic+i
    END DO
    DO j=1,N_phyt
      i=i+1
      iPhyt15N(j)=ic+i
    END DO
    DO j=1,N_zoop
      i=i+1
      iZoop15N(j)=ic+i
    END DO
#  endif
# endif
#endif
#if defined COT_STARFISH
    i=i+1
    iCOTe=ic+i
    i=i+1
    iCOTl=ic+i
#endif


!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------

    Nid=i
!---------------------------------------------------------------------


   allocate( dz(LBi:UBi,LBj:UBj,N) )

   allocate( C(LBi:UBi, LBj:UBj, N, 1, Nid) , &
             dC_dt(LBi:UBi, LBj:UBj, N, Nid)   )
   allocate( Co(LBi:UBi, LBj:UBj, Nid) )
   
   allocate( p_coral (2,LBi:UBi, LBj:UBj) , &
             p_sgrass(LBi:UBi, LBj:UBj) , &
             p_sand  (LBi:UBi, LBj:UBj) , &
             p_algae (LBi:UBi, LBj:UBj)   )   


!-----------------------------------------------------------------------
!  Set initial conditions
!-----------------------------------------------------------------------

    do j=LBj,UBj
      do i=LBi,UBi
        do k=1,N
        
          dz(i,j,k)=d0 !(m)
        
          C(i,j,k,1,iTemp) = T0       
          C(i,j,k,1,iSalt) = S0

          C(i,j,k,1,iSedi) = Sed1_0    !Sediment concentration (g m-3) 0.e0, 1.e0

          C(i,j,k,1,iTIC_) = TIC_0     !DIC  (umol kg-1) 1915.0d0 1930.0d0, 1700.0d0, 2100.0d0, 2300.0d0, 2500.0d0, 2700.0d0, 3000.0d0 
          C(i,j,k,1,iTAlk) = TAlk0      !TA  (umol kg-1)  2232.0d0 2275.0d0, 2500.0d0, 2150.0d0, 2000.0d0, 1800.0d0, 1700.0d0, 1600.0d0 

!          C(i,j,k,1,iOxyg) = Oxyg0      !DO  (umol L-1)
          C(i,j,k,1,iOxyg) = O2satu(C(i,j,k,1,iTemp)+273.15d0, C(i,j,k,1,iSalt))

! Offshore values
          Co(i,j,iTemp) = T0   !27.0d0 32.0d0
          Co(i,j,iSalt) = S0

          Co(i,j,iSedi) = 0.0d0    !Sediment concentration (g m-3) 0.e0, 1.e0

          Co(i,j,iTIC_) = TIC_0     !DIC  (umol kg-1) 1915.0d0 1930.0d0, 1700.0d0, 2100.0d0, 2300.0d0, 2500.0d0, 2700.0d0, 3000.0d0 
          Co(i,j,iTAlk) = TAlk0     !TA  (umol kg-1)  2232.0d0 2275.0d0, 2500.0d0, 2150.0d0, 2000.0d0, 1800.0d0, 1700.0d0, 1600.0d0 

          Co(i,j,iOxyg) = O2satu(C(i,j,k,1,iTemp)+273.15d0, C(i,j,k,1,iSalt))
#if defined ORGANIC_MATTER
          C(i,j,k,1,iDOC(1)) = DOC1_0       !DOC  (umol L-1) 
          C(i,j,k,1,iDOC(2)) = DOC2_0       !DOC  (umol L-1) 
          C(i,j,k,1,iPOC(1)) = POC1_0       !POC  (umol L-1) 
          C(i,j,k,1,iPOC(2)) = POC2_0       !POC  (umol L-1) 
          C(i,j,k,1,iPhyt(1)) = Phyt1_0     !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          C(i,j,k,1,iPhyt(2)) = Phyt2_0     !Phytoplankton2 0.3(umolC L-1) 							diatom
          C(i,j,k,1,iPhyt(3)) = Phyt3_0     !Phytoplankton2 0.3(umolC L-1) 							diatom
          C(i,j,k,1,iZoop(1)) = Zoop1_0     !Zooplankton (umol L-1)1.3
          C(i,j,k,1,iPIC(1)) =  PIC1_0      !Particulated inorganic carbon (umol L-1)1.3
# if defined FOODWEB_TESTMODE
          C(i,j,k,1,iDeadPhyt(1)) = 0.0d0   !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          C(i,j,k,1,iDeadPhyt(2)) = 0.0d0   !Phytoplankton2 0.3(umolC L-1) 							diatom
          C(i,j,k,1,iDeadPhyt(3)) = 0.0d0   !Phytoplankton2 0.3(umolC L-1) 							diatom
# endif
! Offshore values
          Co(i,j,iDOC(1)) = DOC1_0       !DOC  (umol L-1) 
          Co(i,j,iDOC(2)) = DOC2_0       !DOC  (umol L-1) 
          Co(i,j,iPOC(1)) = POC1_0       !POC  (umol L-1) 
          Co(i,j,iPOC(2)) = POC2_0       !POC  (umol L-1) 
          Co(i,j,iPhyt(1)) = Phyt1_0       !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          Co(i,j,iPhyt(2)) = Phyt2_0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          Co(i,j,iPhyt(3)) = Phyt3_0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          Co(i,j,iZoop(1)) = Zoop1_0      !Zooplankton (umol L-1)1.3
          Co(i,j,iPIC(1)) =  PIC1_0      !Particulated inorganic carbon (umol L-1)1.3
# if defined FOODWEB_TESTMODE
          Co(i,j,iDeadPhyt(1)) = 0.0d0   !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          Co(i,j,iDeadPhyt(2)) = 0.0d0   !Phytoplankton2 0.3(umolC L-1) 							diatom
          Co(i,j,iDeadPhyt(3)) = 0.0d0   !Phytoplankton2 0.3(umolC L-1) 							diatom
# endif         
#endif
#if defined CARBON_ISOTOPE
          R13C=R13C_fromd13C(d13C_TIC0)
!          Ci(nDI13C,i,j,k)=R13C/(1.+R13C)*Ci(nDIC,i,j,k) !DI13C (umol kg-1)
          C(i,j,k,1,iT13C) =R13C*C(i,j,k,1,iTIC_) !DI13C (umol kg-1) 
! Offshore values
          Co(i,j,iT13C) =R13C*C(i,j,k,1,iTIC_) !DI13C (umol kg-1) 
# if defined ORGANIC_MATTER
          R13C=R13C_fromd13C(d13C_DOC1_0)
          C(i,j,k,1,iDO13C(1)) = R13C*C(i,j,k,1,iDO13C(1))       !DOC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_DOC2_0)
          C(i,j,k,1,iDO13C(2)) = R13C*C(i,j,k,1,iDO13C(2))      !DOC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_POC1_0)
          C(i,j,k,1,iPO13C(1)) = R13C*C(i,j,k,1,iPO13C(1))        !POC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_POC2_0)
          C(i,j,k,1,iPO13C(2)) = R13C*C(i,j,k,1,iPO13C(2))        !POC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_Phyt1_0)
          C(i,j,k,1,iPhyt13C(1)) = R13C*C(i,j,k,1,iPhyt13C(1))     !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          R13C=R13C_fromd13C(d13C_Phyt2_0)
          C(i,j,k,1,iPhyt13C(2)) = R13C*C(i,j,k,1,iPhyt13C(2))       !Phytoplankton2 0.3(umolC L-1) 							diatom
          R13C=R13C_fromd13C(d13C_Phyt3_0)
          C(i,j,k,1,iPhyt13C(3)) = R13C*C(i,j,k,1,iPhyt13C(3))       !Phytoplankton2 0.3(umolC L-1) 							diatom
          R13C=R13C_fromd13C(d13C_Zoop1_0)
          C(i,j,k,1,iZoop13C(1)) = R13C*C(i,j,k,1,iZoop13C(1))       !Zooplankton (umol L-1)1.3
          R13C=R13C_fromd13C(d13C_PIC1_0)
          C(i,j,k,1,iPI13C(1)) = R13C*C(i,j,k,1,iPI13C(1))       !Particulated inorganic carbon (umol L-1)1.3
! Offshore values
          R13C=R13C_fromd13C(d13C_DOC1_0)
          Co(i,j,iDO13C(1)) = R13C*Co(i,j,iDO13C(1))       !DOC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_DOC2_0)
          Co(i,j,iDO13C(2)) = R13C*Co(i,j,iDO13C(2))       !DOC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_POC1_0)
          Co(i,j,iPO13C(1)) = R13C*Co(i,j,iPO13C(1))        !POC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_POC2_0)
          Co(i,j,iPO13C(2)) = R13C*Co(i,j,iPO13C(2))        !POC  (umol L-1) 
          R13C=R13C_fromd13C(d13C_Phyt1_0)
          Co(i,j,iPhyt13C(1)) = R13C*Co(i,j,iPhyt13C(1))       !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          R13C=R13C_fromd13C(d13C_Phyt2_0)
          Co(i,j,iPhyt13C(2)) = R13C*Co(i,j,iPhyt13C(2))       !Phytoplankton2 0.3(umolC L-1) 							diatom
          R13C=R13C_fromd13C(d13C_Phyt3_0)
          Co(i,j,iPhyt13C(3)) = R13C*Co(i,j,iPhyt13C(3))       !Phytoplankton2 0.3(umolC L-1) 							diatom
          R13C=R13C_fromd13C(d13C_Zoop1_0)
          Co(i,j,iZoop13C(1)) = R13C*Co(i,j,iZoop13C(1))       !Zooplankton (umol L-1)1.3
          R13C=R13C_fromd13C(d13C_PIC1_0)
          Co(i,j,iPI13C(1)) =  R13C*Co(i,j,iPI13C(1))       !Particulated inorganic carbon (umol L-1)1.3
# endif
#endif
#if defined NUTRIENTS            
!!!!!!!!!!!!
          C(i,j,k,1,iNO3_) = NO3_0      !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
          C(i,j,k,1,iPO4_) = PO4_0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
          C(i,j,k,1,iNH4_) = NH4_0      !NH4  (umol L-1) 
! Offshore values
          Co(i,j,iNO3_) = NO3_0      !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
          Co(i,j,iPO4_) = PO4_0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
          Co(i,j,iNH4_) = NH4_0       !NH4  (umol L-1) 

# if defined ORGANIC_MATTER
          C(i,j,k,1,iDON(1)) = DON1_0     !DON  (umol L-1) 
          C(i,j,k,1,iDON(2)) = DON2_0     !DON  (umol L-1) 
          C(i,j,k,1,iPON(1)) = PON1_0     !PON  (umol L-1) 
          C(i,j,k,1,iPON(2)) = PON2_0    !PON  (umol L-1) 
          C(i,j,k,1,iDOP(1)) = DOP1_0     !DOP  (umol L-1) 
          C(i,j,k,1,iDOP(2)) = DOP2_0    !DOP  (umol L-1) 
          C(i,j,k,1,iPOP(1)) = POP1_0      !POP  (umol L-1) 
          C(i,j,k,1,iPOP(2)) = POP2_0     !POP  (umol L-1) 
! Offshore values
          Co(i,j,iDON(1)) = DON1_0     !DON  (umol L-1) 
          Co(i,j,iDON(2)) = DON2_0     !DON  (umol L-1) 
          Co(i,j,iPON(1)) = PON1_0     !PON  (umol L-1) 
          Co(i,j,iPON(2)) = PON2_0    !PON  (umol L-1) 
          Co(i,j,iDOP(1)) = DOP1_0     !DOP  (umol L-1) 
          Co(i,j,iDOP(2)) = DOP2_0    !DOP  (umol L-1) 
          Co(i,j,iPOP(1)) = POP1_0      !POP  (umol L-1) 
          Co(i,j,iPOP(2)) = POP2_0     !POP  (umol L-1) 
# endif
# if defined NITROGEN_ISOTOPE            
!!!!!!!!!!!!
          C(i,j,k,1,i15NO3) =  0.2d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
          C(i,j,k,1,i15NH4) =  0.26d0       !NH4  (umol L-1) 
! Offshore values
          Co(i,j,i15NO3) =  0.2d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
          Co(i,j,i15NH4) =  0.26d0       !NH4  (umol L-1) 
#  if defined ORGANIC_MATTER
          C(i,j,k,1,iDO15N(1)) =   1.0d0 !3.0d0       !DON  (umol L-1) 
          C(i,j,k,1,iDO15N(2)) =  10.0d0 !3.0d0       !DON  (umol L-1) 
          C(i,j,k,1,iPO15N(1)) =  0.06d0 !0.05d0      !PON  (umol L-1) 
          C(i,j,k,1,iPO15N(2)) =  0.6d0 !0.05d0      !PON  (umol L-1) 
          C(i,j,k,1,iPhyt15N(1)) = 0.0d0     !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          C(i,j,k,1,iPhyt15N(2)) = 0.0d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          C(i,j,k,1,iPhyt15N(3)) = 0.0d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          C(i,j,k,1,iZoop15N(1)) = 0.0d0       !Zooplankton (umol L-1)1.3
! Offshore values
          Co(i,j,iDO15N(1)) =   1.0d0 !3.0d0       !DON  (umol L-1) 
          Co(i,j,iDO15N(2)) =  10.0d0 !3.0d0       !DON  (umol L-1) 
          Co(i,j,iPO15N(1)) =  0.06d0 !0.05d0      !PON  (umol L-1) 
          Co(i,j,iPO15N(2)) =  0.6d0 !0.05d0      !PON  (umol L-1) 
          Co(i,j,iPhyt15N(1)) = 0.0d0     !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          Co(i,j,iPhyt15N(2)) = 0.0d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          Co(i,j,iPhyt15N(3)) = 0.0d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          Co(i,j,iZoop15N(1)) = 0.0d0       !Zooplankton (umol L-1)1.3
#  endif
# endif
#endif
#if defined COT_STARFISH
          C(i,j,k,1,iCOTe) =  0.0d0       !COTS eggs (umol L-1) 
          C(i,j,k,1,iCOTl) =  0.0d0       !COTS larvae (umol L-1) 
#endif
        enddo
          
#if defined CHAMBER_SITE4
        p_coral(1,i,j) = 0.208d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE5
        p_coral(1,i,j) = 0.265d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE6
        p_coral(1,i,j) = 0.37d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE7
        p_coral(1,i,j) = 0.231d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE9
        p_coral(1,i,j) = 0.0d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.412d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 0.0d0
#elif defined CHAMBER_SITE10
        p_coral(1,i,j) = 0.0d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.486d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 0.0d0

#else
        p_coral(1,i,j) = p_coral1_0
        p_coral(2,i,j) = p_coral2_0
        p_sand(i,j)  = p_sand_0
#endif
        p_algae(i,j) = p_algae_0
        p_sgrass(i,j)= p_sgrass_0
        
        
      enddo
    enddo
    
    ereef = 0.0d0

    return
  end subroutine initialize_params

END MODULE mod_param

