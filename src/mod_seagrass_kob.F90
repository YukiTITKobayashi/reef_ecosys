
!!!=== Copyright (c) 2012-2020 Takashi NAKAMURA  =====
!!!=== Improved by Kobayashi =====
#include "cppdefs.h"


!!!*************** SEAGRASS ***********************************************


MODULE mod_seagrass

    implicit none


    integer, parameter :: Nsg = 1    !! Number of seagrass groups
#if defined SGRASS_KOB
    real(8), parameter :: g_R_per_Pg = 0.2d0    !coefficient of activeness for light-attenuated respiration
#endif

!Parameters not relating to seagrass
#if defined SGRASS_KOB
    real(8), parameter :: SGD(:,:,:)    !C(:,:,:,:)からk=1:Nの次元を除く
    real(8), parameter :: depth(:,:)
    real(8), parameter :: dist_coastline(:,:,:)
    real(8), parameter :: dist_rivermouth(:,:,:)
#endif

    TYPE T_SGRASS
        real(8), pointer :: Pg(:,:,:)   !ng: No. of seagrass group , i: x grid, j: y grid
        real(8), pointer :: R (:,:,:) 
        real(8), pointer :: QC(:,:,:)   !海草が使用可能な溶存炭素量と思われる，単位umol?
#if defined SGRASS_KOB
        real(8), pointer :: R_dark(:,:,:)
        real(8), pointer :: R_light(:,:,:)
        real(8), pointer :: F_sap(:,:,:)        ![m s-1] Sap flow inside seagrass driven by salinity difference
#endif
#if defined CARBON_ISOTOPE
!  13C isotope
        real(8), pointer :: Q13C(:,:,:)
#endif
#if defined NUTRIENTS         
        real(8), pointer :: QN(:,:,:)
        real(8), pointer :: QP(:,:,:)
#endif
#if defined SGRAASS_SIZE_DYNAMICS_KOB
        !保留．dzを使ってもっとうまく書けるはず
        !real(8), pointer :: Leaf_length(:,:,:)      ![m]
        !integer, pointer :: div_leaf(:,:,:)         ![m]
        !integer, pointer :: part_leaf(:,:,:)
        !real(8), pointer, allocatable :: I_d(:,:,:,:)              ![umol m-2 s-1]   arg4=leaf_pos
        !real(8), pointer :: I_d_ave(:,:,:)
        !real(8), pointer, allocatable :: Pg_part(:,:,:,:)          ![mmolC/m2/s] arg4=leaf_pos
        real(8), pointer :: growth(:,:,:)  ! growth rate
        real(8), pointer :: mort(:,:,:)  ! mortaliny rate
#endif

    END TYPE T_SGRASS

    TYPE (T_SGRASS), allocatable :: SGRASS(:)

CONTAINS

!!! **********************************************************************
!!!  set initial conditions
!!! **********************************************************************

subroutine initialize_seagrass(ng, Ngrids, LBi, UBi, LBj, UBj)

    USE mod_geochem
    
    implicit none
! input parameters
    integer, intent(in) :: ng, Ngrids, LBi, UBi, LBj, UBj
    real(8)  R13C
    integer i,j,n

    IF (ng.eq.1) allocate ( SGRASS(Ngrids) )

    allocate( SGRASS(ng)%Pg(Nsg,LBi:UBi,LBj:UBj)     )
    allocate( SGRASS(ng)%R (Nsg,LBi:UBi,LBj:UBj)     )
    allocate( SGRASS(ng)%QC(Nsg,LBi:UBi,LBj:UBj)     )
#if defined SGRASS_KOB
    allocate( SGRASS(ng)%Leaf_length(Nsg,LBi:UBi,LBj:UBj)      )
    !allocate( SGRASS(ng)%div_leaf(Nsg,LBi:UBi,LBj:UBj)      )
    !allocate( SGRASS(ng)%part_leaf(Nsg,LBi:UBi,LBj:UBj)      )
    !allocate( SGRASS(ng)%I_d(Nsg,LBi:UBi,LBj:UBj,part_leaf)      )
    !allocate( SGRASS(ng)%I_d_ave(Nsg,LBi:UBi,LBj:UBj)      )
    !allocate( SGRASS(ng)%Pg_part(Nsg,LBi:UBi,LBj:UBj,part_leaf)      )
    allocate( SGRASS(ng)%R_dark(Nsg,LBi:UBi,LBj:UBj)      )
    allocate( SGRASS(ng)%R_light(Nsg,LBi:UBi,LBj:UBj)     )
    allocate( SGRASS(ng)%F_sap(Nsg,LBi:UBi,LBj:UBj)       )
#endif
#if defined CARBON_ISOTOPE
    allocate( SGRASS(ng)%Q13C(Nsg,LBi:UBi,LBj:UBj)   )
#endif
#if defined NUTRIENTS         
    allocate( SGRASS(ng)%QN(Nsg,LBi:UBi,LBj:UBj)     )
    allocate( SGRASS(ng)%QP(Nsg,LBi:UBi,LBj:UBj)     )
#endif

!------------------------------------------
!  Set initial conditions
    do j=LBj,UBj
        do i=LBi,UBi
            do n=1,Nsg
!            seagrass internal conditions
                SGRASS(ng)%Pg(n,i,j) = 0.0d0
                SGRASS(ng)%R (n,i,j) = 0.0d0
                SGRASS(ng)%QC(n,i,j) = 15.0d0  !!!�Ă��Ƃ�
#if defined SGRASS_KOB
                SGRASS(ng)%Leaf_length(n,i,j) = 0.8d0
                !SGRASS(ng)%Leaf_length(n,i,j) = 8
                !SGRASS(ng)%div_leaf(n,i,j)    = Leaf_length/part_leaf
                SGRASS(ng)%R_dark(n,i,j)  = 0.0d0
                SGRASS(ng)%R_light(n,i,j) = 0.0d0
                SGRASS(ng)%F_sap(n,i,j)   = 0.0d0
#endif
#if defined CARBON_ISOTOPE
                R13C = R13C_fromd13C(-15.0d0)
!               c13CH2O (n,i,j)=R13C/(1.+R13C)*CH2O(n,i,j)
                SGRASS(ng)%Q13C(n,i,j) = R13C * SGRASS(ng)%QC(n,i,j)
#endif
#if defined NUTRIENTS         
                SGRASS(ng)%QN(n,i,j) = 1.5d0  !!!�Ă��Ƃ�
                SGRASS(ng)%QP(n,i,j) = 0.1d0  !!!�Ă��Ƃ�
#endif
            enddo

#if defined SGRASS_KOB
            dist_coastline(i,j,1) = 
#endif

        enddo
    enddo
    
    RETURN
    
END SUBROUTINE initialize_seagrass

!!! **********************************************************************
!!!  Main program of seagrass model 
!!! **********************************************************************

SUBROUTINE seagrass           &
!   input parameters
    ( ng, n, i, j    &   ! ng: nested grid number; n: seagrass compartment; i,j: position
    , PFD            &   ! Photon flux density (umol m-2 s-1)
    , rho_sw         &   ! Density of seawater (g cm-3)
    , DICamb         &   ! DIC (umol kg-1)
    , DOamb          &   ! DO (umol L-1)
#if defined NUTRIENTS         
    , NH4amb         &   ! NH4 concentration (umol L-1)
#endif
#if defined CARBON_ISOTOPE
    , DI13Camb       &   ! 13C of DIC (umol kg-1)
#endif
!   output parameters
    , DICuptake      &   ! DIC uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
    , DOuptake       &   ! DO  uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
#if defined NUTRIENTS         
    , NO3uptake      &   ! NO3 uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
    , NH4uptake      &   ! NH4 uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
    , PO4uptake      &   ! PO4 uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
#endif
#if defined CARBON_ISOTOPE
    , DI13Cuptake    &   ! DI13C uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
#endif
    )
END SUBROUTINE seagrass
    
    !-----------------------------------------------------------------------
!作業メモ：DICの出元を探して，表層水，底層水，間隙水のプロファイルを取得したい


END MODULE mod_seagrass