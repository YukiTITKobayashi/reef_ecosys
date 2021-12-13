
!!!=== Copyright (c) 2017-2020 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** box model MODULE ************************************

MODULE mod_reef_flow
  implicit none
  TYPE T_REEF
    real(8), pointer :: Wir(:,:)
    real(8), pointer :: Dir(:,:)
    real(8), pointer :: Lir(:,:)
    real(8), pointer :: Wrc(:,:)
    real(8), pointer :: Drc(:,:)
    real(8), pointer :: Lrc(:,:)
    real(8), pointer :: Wch(:,:)
    real(8), pointer :: Dch(:,:)
    real(8), pointer :: Lch(:,:)
    real(8), pointer :: Slp(:,:)
    
    real(8), pointer :: Air(:,:)
    real(8), pointer :: Qrc(:,:) 
    real(8), pointer :: Qch(:,:) 
    real(8), pointer :: el(:,:) 
  END TYPE T_REEF
  TYPE (T_REEF), allocatable :: REEF(:)
CONTAINS

!!! **********************************************************************
!!!  set initial conditions for reef hydrodynamics
!!! **********************************************************************

  SUBROUTINE initialize_reef_flow(ng, Ngrids, LBi, UBi, LBj, UBj)

    implicit none
! input parameters
    integer, intent(in) :: ng, Ngrids, LBi, UBi, LBj, UBj
    integer i,j,n

    IF (ng.eq.1) allocate ( REEF(Ngrids) )

    allocate( REEF(ng)%Wir(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Dir(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Lir(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Wrc(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Drc(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Lrc(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Wch(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Dch(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Lch(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Slp(LBi:UBi,LBj:UBj)  )
    
    allocate( REEF(ng)%Air(LBi:UBi,LBj:UBj)  )
    
    allocate( REEF(ng)%Qrc(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%Qch(LBi:UBi,LBj:UBj)  )
    allocate( REEF(ng)%el (LBi:UBi,LBj:UBj)  )

!  S initial conditions
    do j=LBj,UBj
      do i=LBi,UBi
        REEF(ng)%Wir(i,j)=800.0d0   ! Inner reef width (m)
        REEF(ng)%Dir(i,j)=1.0d0     ! Inner reef depth (m)
        REEF(ng)%Lir(i,j)=700.0d0   ! Inner reef length (m)
        REEF(ng)%Wrc(i,j)=750.0d0   ! Reef crest width (m)
        REEF(ng)%Drc(i,j)=0.0d0     ! Reef crest depth (m)
        REEF(ng)%Lrc(i,j)=100.0d0   ! Reef crest length (m)
        REEF(ng)%Wch(i,j)=REEF(ng)%Wir(i,j)-REEF(ng)%Wrc(i,j)   ! Channel width (m)
        REEF(ng)%Dch(i,j)=0.8d0               ! Channel depth (m)
        REEF(ng)%Lch(i,j)=REEF(ng)%Lrc(i,j)   ! Channel length (m)
        REEF(ng)%Slp(i,j)=0.3   ! Reef slope
        
        REEF(ng)%Air(i,j)=REEF(ng)%Wir(i,j)*REEF(ng)%Lir(i,j)   ! Inner reef area (m2)

        REEF(ng)%Qrc(i,j)=0.0d0   ! Volume flux per unit width through the reef crest (m s-1)
        REEF(ng)%Qch(i,j)=0.0d0   ! Volume flux per unit width through the channel (m s-1)
        REEF(ng)%el (i,j)=0.5d0   ! Seawater elevation inside the reef (m)
      enddo
    enddo

    RETURN
  END SUBROUTINE initialize_reef_flow


! **********************************************************************
!  Mass balance of water column
! **********************************************************************

  SUBROUTINE reef_flow &
!   input parameters
    ( ng, i, j       & ! ng: nested grid number; i,j: position
    , dt             & ! Time step (sec)
    , Hs_o           & ! Significant wave hight at offshore (m)
    , Ts             & ! Significant Wave period (s)
    , el_o           & ! Sea surface elevation at offshore (m)
    )
!
!-----------------------------------------------------------------------
!                                                                       
!                   ___________ _____                                   
!                  |           :     :                                  
!                  |           :CHAN-:Wch                               
!                  |           : NEL<-- Qch                             
!                  |           : Dch :                                  
!        SHORELINE |           :_____:                                  
!                  |   INNER   |     |     OFFSHORE                     
!               Wir|   REEF    |REEF |                                  
!                  |           |CREST|Wrc                               
!                  |    Dir,   |    <-- Qrc                             
!                  |    el     | Dcr |                                  
!                  |           |     |                                  
!                  |___________|____ |                                  
!                        Lir     Lrc                                    
!                                                                       
!      Plan view of an idealized reef-channel system.                   
!-----------------------------------------------------------------------
!
    implicit none

! input parameters
    integer, intent(in ) :: ng, i, j        
    real(8), intent(in ) :: dt         
    real(8), intent(in ) :: Hs_o
    real(8), intent(in ) :: Ts
    real(8), intent(in ) :: el_o

!-----------------------------------------------------------------------
    real(8), parameter :: rho = 1024.0d0   ! Seawater density (kg m-3)
    real(8), parameter :: g   = 9.80665d0  ! Gravitational acceleration (m s-2)
    real(8), parameter :: pi  = 3.14159265359d0  ! Circle ratio
    
    real(8), parameter :: dx_o  = 50.0d0   ! Surf zone length (m)

!-----------------------------------------------------------------------
    real(8) :: d_o          ! Water depth (m)
    real(8) :: d_i          ! Water depth (m)
    real(8) :: del          ! difference of elevation between inner and outer reef (m)
    
    d_o = el_o + REEF(ng)%Drc(i,j)
!    d_o = el_o + REEF(ng)%Slp(i,j)*dx_o
    d_i = REEF(ng)%el (i,j) + REEF(ng)%Drc(i,j)
    del = REEF(ng)%el (i,j) - el_o

    CALL momentum_eq      &
!     input parameters
      ( dt                &   ! Time step (sec)
      , REEF(ng)%Lrc(i,j) &   ! Distance (m)
!      , dx_o &   ! Distance (m)
      , Hs_o              &   ! Significant wave hight at offshore (m)
      , Ts                &   ! Significant Wave period (s)
      , d_o               &   ! Depth at offshore side edge (m)
      , d_i               &   ! Depth at inshore side edge (m)
!      , del               &   ! Differnce of elevation (m)
!     input and output parameters
      , REEF(ng)%Qrc(i,j) &   ! Volume flux per unit width (m2 s-1) (o->i is positive, m s-1)
      )
    
    d_o = el_o + REEF(ng)%Dch(i,j)
    d_i = REEF(ng)%el (i,j) + REEF(ng)%Dch(i,j)

    CALL momentum_eq      &
!     input parameters
      ( dt                &   ! Time step (sec)
      , REEF(ng)%Lch(i,j) &   ! Distance (m)
!      , dx_o &   ! Distance (m)
      , Hs_o              &   ! Significant wave hight at offshore (m)
      , Ts                &   ! Significant Wave period (s)
      , d_o               &   ! Depth at offshore side edge (m)
      , d_i               &   ! Depth at inshore side edge (m)
!      , del               &   ! Differnce of elevation (m)
!     input and output parameters
      , REEF(ng)%Qch(i,j) &   ! Volume flux per unit width (m2 s-1) (o->i is positive, m s-1)
      )
    
    ! Mass balance equation
    REEF(ng)%el(i,j) = REEF(ng)%el(i,j) + (                           &
                         REEF(ng)%Qrc(i,j)*REEF(ng)%Wrc(i,j)          &
                        +REEF(ng)%Qch(i,j)*REEF(ng)%Wch(i,j)          &
                       )/REEF(ng)%Air(i,j) * dt
    write(98,*) el_o, REEF(ng)%el (i,j),REEF(ng)%Qrc(i,j),REEF(ng)%Qch(i,j)

  END SUBROUTINE reef_flow
    
! **********************************************************************
!  Momentum equation
! **********************************************************************

  SUBROUTINE momentum_eq &
!   input parameters
    ( dt             & ! Time step (sec)
    , dx             & ! Distance (m)
    , Hs_o           & ! Significant wave hight at offshore (m)
    , Ts             & ! Significant Wave period (s)
    , d_o            & ! Depth at offshore side edge (m)
    , d_i            & ! Depth at inshore side edge (m)
!    , del            & ! Differnce of elevation (m)

!   input and output parameters
    , q              & ! Volume flux per unit width (m2 s-1) (o->i is positive, m s-1)
    )
!-----------------------------------------------------------------------
!
    implicit none

! input parameters
    real(8), intent(in ) :: dt         
    real(8), intent(in ) :: dx
    real(8), intent(in ) :: Hs_o
    real(8), intent(in ) :: Ts
    real(8), intent(in ) :: d_o
    real(8), intent(in ) :: d_i
!    real(8), intent(in ) :: del
! input and output parameters
    real(8), intent(inout) :: q

!-----------------------------------------------------------------------
    real(8), parameter :: rho = 1024.0d0   ! Seawater density (kg m-3)
    real(8), parameter :: g   = 9.80665d0  ! Gravitational acceleration (m s-2)
    real(8), parameter :: pi  = 3.14159265359d0  ! Circle ratio
    real(8), parameter :: Cd  = 0.02d0  ! Bottom drag coefficient (Lowe et al., 2009)
    real(8), parameter :: dcrt= 0.01d0  ! Critical depth for Weting and Drying (m)

    real(8) :: Hs_i         ! Significant wave hight (m)
    real(8) :: k            ! Wave number
    real(8) :: cp           ! Wave phase verocity (m s-1)
    real(8) :: cg           ! Wave group verocity (m s-1)
    real(8) :: L            ! Wave length (m)
    real(8) :: Ew_o         ! Wave energy (J)
    real(8) :: Ew_i         ! Wave energy (J)
    real(8) :: Sxx_o        ! Radiation stress in the reef (N m-1)
    real(8) :: Sxx_i        ! Radiation stress in the reef (N m-1)
    real(8) :: del          ! difference of elevation between inner and outer reef (m)
    real(8) :: dSxx         ! difference of radiation stress between inner and outer reef (N m-2)
    real(8) :: d_m          ! mean depth (m)
    real(8) :: tau          ! Bottom shear stress (N m-2)
    real(8) :: u            ! mean velocity (N m-2)
    real(8) :: d_o2         ! Depth at offshore side edge (m)
    real(8) :: d_i2         ! Depth at inshore side edge (m)
          
!-----------------------------------------------------------------------

    if (d_o <= dcrt) then
      d_o2  = dcrt
      Sxx_o = 0.0d0
    else
      d_o2  = d_o
      CALL wavelength_from_T_h(Ts, d_o, L)  ! wave length
!      L = 0.5d0*g*Ts*Ts*pi  ! wave length of Deep water wave
      k = 2.0d0*pi/L
      Ew_o = 0.125d0*rho*g*Hs_o*Hs_o
      Sxx_o = Ew_o * (2.0d0*k*d_o/(sinh(2*k*d_o)) + 0.5)
    endif
    
    if (d_i <= dcrt) then
      d_i2   = dcrt
      Sxx_i = 0.0d0
    else
      d_i2 = d_i
      CALL wavelength_from_T_h(Ts, d_i, L)  ! wave length
      k = 2.0d0*pi/L
!      Hs_i = min(0.8d0*d_i, Hs_o)
      Hs_i = min(0.8d0*d_o, Hs_o)
      Ew_i = 0.125d0*rho*g*Hs_i*Hs_i
      Sxx_i = Ew_i * (2.0d0*k*d_i/(sinh(2*k*d_i)) + 0.5)
    endif
    
    del = d_i2 - d_o2
    dSxx = Sxx_i - Sxx_o
    d_m  = 0.5d0*(d_i2 + d_o2)
!    if (d_m <= 0.01d0) then
!      u  = 0.0d0
!    else
!      u = q/d_m
!    endif
    u = q/d_m
    tau = rho*Cd*abs(u)*u
    
    ! Momentum equation
    q = q + (-g*d_m*del/dx -dSxx/rho/dx -tau/rho)*dt
!    q = q + (-g*d_m*del/dx -dSxx/rho/dx)*dt
    
    write(99,*) q, d_i2,d_o2,Sxx_i,Sxx_o,u, d_m, del, tau

  END SUBROUTINE momentum_eq
    
  SUBROUTINE wavelength_from_T_h &
!   input parameters
    ( T              & ! Wave period (s)
    , h              & ! Depth (m)
!   output parameters
    , L              & ! Wave length
    )
!-----------------------------------------------------------------------
!
    implicit none

! input parameters
    real(8), intent(in ) :: T         
    real(8), intent(in ) :: h
! input and output parameters
    real(8), intent(out) :: L

!-----------------------------------------------------------------------
    real(8), parameter :: pi  = 3.14159265359d0  ! Circle ratio
    real(8), parameter :: g   = 9.80665d0  ! Gravitational acceleration (m s-2)
    real(8), parameter :: err = 0.01d0  ! calculation error (m)
    
    real(8) :: L2
!-----------------------------------------------------------------------
!
!   c= L/T=sqrt( (g*L/2/pi*tanh(2*pi*h/L)) ) 
    
    ! Initial wave length
    L = 0.5d0*g*T*T*pi  ! wave length of Deep water wave
    L2 = sqrt(0.5d0*g*L/pi*tanh(2.0d0*pi*h/L))*T
    do while (abs(L-L2)>err)
      L=L2
      L2 = sqrt(0.5d0*g*L/pi*tanh(2.0d0*pi*h/L))*T
    end do
    L=L2

  END SUBROUTINE wavelength_from_T_h

END MODULE mod_reef_flow


