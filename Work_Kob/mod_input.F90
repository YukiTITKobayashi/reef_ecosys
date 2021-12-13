
!!!=== Copyright (c) 2012-2020 Takashi NAKAMURA  =====

!--------------------------------------------------------------------------------
!
!              Input module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

module mod_input

contains
! **********************************************************************
!  Read timeseries data file
! **********************************************************************

  subroutine read_timeseies
  
    USE mod_param
    
    implicit none
    
    integer :: nm_data, ios
    integer i,j
!
! ----- READ Tide data (m) ---------------------------------------------

    open(77,file='./input/tide2010_0820-0930.dat')
    dt_tide=1.d0  !1.0 hour interval
      
!  ---- count data number -----
    nm_data = 0
    do
      read(77,*,iostat=ios) 
      if(ios==-1) exit
      nm_data = nm_data+1
    end do
    allocate( tide_data(nm_data) )
    nm_tide = nm_data
    rewind(77)
    do i=1, nm_data
      read(77,*) tide_data(i)
    end do
    close(77)

! ----- READ meterological dataset --------------------------------

    open(77,file='./input/air2010_0820-1001.dat')
    dt_air=1.d0  !1.0 hour interval
      
!  ---- count data number -----
    nm_data = 0
    do
      read(77,*,iostat=ios) 
      if(ios==-1) exit
      nm_data = nm_data+1
    end do
    allocate( Pair_data (nm_data))
    allocate( Uwind_data(nm_data))
    allocate( Vwind_data(nm_data))
    allocate( Tair_data (nm_data))
    allocate( Qair_data (nm_data))
    allocate( rain_data (nm_data))
    allocate( swrad_data(nm_data))
    allocate( cloud_data(nm_data))
    nm_air = nm_data
    rewind(77)

    do i=1, nm_data
      read(77,*) Pair_data (i), Uwind_data(i), Vwind_data(i),    &
                 Tair_data (i), Qair_data (i), rain_data (i),    &
                 swrad_data(i), cloud_data(i)
    end do
    close(77)

#ifdef LONGWAVE_IN
      
! ----- READ downward longwave radiation (W m-2) ---------------------------------------------

    open(77,file='./input/dwlw2010_0820-1001.dat')
    dt_dlwrad=1.d0  !1.0 hour interval
      
!  ---- count data number -----
    nm_data = 0
    do
      read(77,*,iostat=ios)
      if(ios==-1) exit
      nm_data = nm_data+1
    end do
    allocate( dlwrad_data(nm_data))
    nm_dlwrad = nm_data
    rewind(77)

    do i=1, nm_data
      read(77,*) dlwrad_data(i)
    end do
    close(77)
#endif

! --------- READ Wave data ---------------------------------------------

!    open(77,file='./input/wave2010_0820-0905test.dat')
    open(77,file='./input/wave2010_0820-0905.dat')
    dt_wave=0.25d0  !15 min interval

!  ---- count data number -----
    nm_data = 0
    do
      read(77,*,iostat=ios)
      if(ios==-1) exit
      nm_data = nm_data+1
    end do
    allocate( Hs_data(nm_data), Tp_data(nm_data) )
    nm_wave = nm_data
    rewind(77)

    do i=1, nm_data
      read(77,*) Hs_data(i), Tp_data(i)
    end do
    close(77)

! --------- READ Offshre temperature data ---------------------------------------------

    open(77,file='./input/offshoretemp.dat')
    dt_offtemp=3.0d0  !3 hour interval

!  ---- count data number -----
    nm_data = 0
    do
      read(77,*,iostat=ios) 
      if(ios==-1) exit 
      nm_data = nm_data+1
    end do
    allocate( offtemp_data(nm_data) )
    nm_offtemp = nm_data
    rewind(77)

    do i=1, nm_data
      read(77,*) offtemp_data(i)
    end do
    close(77)
 
!-----------READ river data---------------------------------------------


    return

  end subroutine read_timeseies


! **********************************************************************
!  Read chamber condition data file
! **********************************************************************

  subroutine read_chambercondition
  
    USE mod_param
    
    implicit none
    
    integer, parameter :: N_rep = 5    ! Repeat number of same light condition
    integer :: N_data
    integer :: ios
    
    real(8), parameter :: c1 = -2.9339d-3
    real(8), parameter :: c2 = 10.326d0
    real(8), parameter :: c3 = -8787.5d0
    integer :: i,j
!
      
! ===== READ PPFD data ======================================================
!         time (hour), PPFD (umol m-2 s-1)
#if defined CHAMBER_SITE4
    open(77,file='./input/pfd_04.txt')
#elif defined CHAMBER_SITE5
    open(77,file='./input/pfd_05.txt')
#elif defined CHAMBER_SITE6
    open(77,file='./input/pfd_06.txt')
#elif defined CHAMBER_SITE7
    open(77,file='./input/pfd_07.txt')
#elif defined CHAMBER_SITE9
    open(77,file='./input/pfd_09.txt')
#elif defined CHAMBER_SITE10
    open(77,file='./input/pfd_10.txt')
#endif
! ----- Count data -----
    N_PFD = 0
    do
      read(77,*,iostat=ios)
      if(ios==-1) exit
      N_PFD = N_PFD + 1
    end do
    allocate( PFD_time(N_PFD*N_rep), PFD_data(N_PFD*N_rep) )   ! Repeat same light condition for 5 times
    rewind(77) 
! ----- Read data -----
    do j=1, N_rep
      do i=1, N_PFD
        read(77,*) PFD_time(i+N_PFD*(j-1)), PFD_data(i+N_PFD*(j-1))
        if(j >= 2) then
          PFD_time(i+N_PFD*(j-1)) = PFD_time(i+N_PFD*(j-1)) + PFD_time(N_PFD)*dble(j-1)
        end if
      end do
      rewind(77) 
    end do
    N_PFD = N_PFD*N_rep
    close(77)

! ===== READ TA and DIC data =================================================
!         time (hour), TA & DIC (umol kg-1)
#if defined CHAMBER_SITE4
    open(77,file='./input/site04.txt')
#elif defined CHAMBER_SITE5
    open(77,file='./input/site05.txt')
#elif defined CHAMBER_SITE6
    open(77,file='./input/site06.txt')
#elif defined CHAMBER_SITE7
    open(77,file='./input/site07.txt')
#elif defined CHAMBER_SITE9
    open(77,file='./input/site09.txt')
#elif defined CHAMBER_SITE10
    open(77,file='./input/site10.txt')
#endif
! ----- Count data -----
    N_WQ = 0
    do
      read(77,*,iostat=ios)
      if(ios==-1) exit
      N_WQ = N_WQ + 1
    end do
    allocate( WQ_time(N_WQ), TA_data(N_WQ), DIC_data(N_WQ), DO_data(N_WQ) )
    rewind(77) 
! ----- Read data -----
    do i=1, N_WQ
      read(77,*) WQ_time(i), TA_data(i), DIC_data(i)
      WQ_time(i) = WQ_time(i) + 4.0d0*24.0d0
      
      if( DIC_data(i) > -0.25d0*c2/c1) then
        DO_data(i) = c1*DIC_data(i)*DIC_data(i) + c2*DIC_data(i) + c3
      else
        DO_data(i) = -(c2*c2-4.0d0*c1*c3)*0.25/c1
      end if
    end do
    close(77)
    
    return

  end subroutine read_chambercondition

!
! **********************************************************************
!  Set environmental condition
! **********************************************************************

  subroutine setdata(nSetting)
!   Setting of condition (nSetting)
!
!   nSetting = 1: Stable condition
!              2: Closed chamber condition
!              3: Constant flow condition
!              4: Reef simulation condition
!              5: Incubation chamber condition
  
    USE mod_param
    
    implicit none
    
    integer, intent(in) :: nSetting

!  -- Set Tide data ---------------------------------------------
    tide=0.0d0
!    tide=lin_interpol(time,tide_data,dt_tide,nm_tide)
! -- Set wind data -----------------------------------------------

    Uwind=0.
    Vwind=0.
!    Uwind=lin_interpol(time,Uwind_data,dt_air,nm_air)
!    Vwind=lin_interpol(time,Vwind_data,dt_air,nm_air)


    if (nSetting .eq. 5) then
    
      CALL lin_interpol2(time*24.0d0,PFD_time,PFD_data,N_PFD, PFDsurf, i_PFD)
    
    else

! -- Set solar radiation data -----------------------------------------------

!      swrad=lin_interpol(time,swrad_data,dt_air,nm_air)   !data from file

!      swrad=solar_radi(time,1400.d0,0.d0)!9.d0/24.d0)                !Artificial solar radiation

!      swrad=light_and_dark(time, 140.0d0, 10./60./24., 0.5/60./24.)   !light and dark method 10 min interval
!      swrad=light_and_dark(time, 350.0d0, 30./60./24., 0.5/60./24.)   !light and dark method 1 hour interval

!      swrad=light_and_dark2(time, 0.0d0, 3.0d0/24.0d0, 21.0d0/24.0d0, 0.5d0/60.0d0/24.0d0)   !light and dark method 1 hour interval
!      swrad=light_and_dark2(time, 100.0d0, 3.0d0/24.0d0, 21.0d0/24.0d0, 0.5d0/60.0d0/24.0d0)   !light and dark method 1 hour interval
!      swrad=light_and_dark2(time, 250.0d0, 3.0d0/24.0d0, 21.0d0/24.0d0, 0.5d0/60.0d0/24.0d0)   !light and dark method 1 hour interval
!      swrad=light_and_dark2(time, 500.0d0, 3.0d0/24.0d0, 21.0d0/24.0d0, 0.5d0/60.0d0/24.0d0)   !light and dark method 1 hour interval
!      swrad=light_and_dark2(time, 1000.0d0, 3.0d0/24.0d0, 21.0d0/24.0d0, 0.5d0/60.0d0/24.0d0)   !light and dark method 1 hour interval

!      swrad=light_and_dark3(time, 3.0d0/24.0d0, 21.0d0/24.0d0, 0.5d0/60.0d0/24.0d0)   !light and dark method 1 hour interval

!      swrad=short_radi(time, 0.0d0, 1.0d0, 24.0d0, 25.0d0, 50.0d0, 1) !shortwave radiation by Zillman equation
      swrad=short_radi(time, 0.0d0, 77.0d0, 24.0d0, 27.0d0, 50.0d0, 2) !shortwave radiation around 3/21 by Zillman equation
!      swrad=short_radi(time, 0.0d0, 77.0d0, 24.0d0, 27.0d0, 50.0d0, 2)*0.25d0 !shortwave radiation around 3/21 by Zillman equation
!      swrad=short_radi(time, 0.0d0, 77.0d0, 24.0d0, 27.0d0, 50.0d0, 2)*0.5d0 !shortwave radiation around 3/21 by Zillman equation
    
!   Convert solar radiation (W m-2) to photosynthetic photon flux 
!   density: 2.1 umol m-2 s-1 per W m-2 (Britton and Dodd 1976)
!   Sea surface albedo: 0.07

      PFDsurf = 2.1d0 * swrad*(1.0d0-0.07d0)
    end if

!----- Temperature controle -----------------------------------------

!    if(time > 5.0d0) then
!
!      do j=1,Jm
!        do i=1,Im
!          do k=1,N
!            C(i,j,k,1,iTemp) = 27.0d0   !27.0d0 32.0d0
!          enddo
!        enddo
!      enddo
!    
!    end if
      
!----- Flux calculation -----------------------------------------

!    tau = 1024*0.01*0.0**2. *0.5 !densSW*Cd*Ub**2    (0 cm s-1)
!    tau = 1024*0.01*0.02**2. *0.5  !densSW*Cd*Ub**2  (2 cm s-1)
    tau = 1024*0.01*0.01**2. *0.5  !densSW*Cd*Ub**2  (1 cm s-1)

!    tau = 1024*0.14*0.02**2. *0.5  !densSW*Cd*Ub**2  (2 cm s-1)
!    tau = 1024*0.14*0.05**2. *0.5  !densSW*Cd*Ub**2  (5 cm s-1)
!    tau = 1024*0.14*0.10**2. *0.5  !densSW*Cd*Ub**2  (10 cm s-1)

!    tau = 1024*0.01*0.2**2. *0.5 !densSW*Cd*Ub**2   (20 cm s-1)
    
    pCO2air = 370.0d0 !(uatm)
    
    U10 = SQRT( Uwind*Uwind + Vwind*Vwind )
    
    fvol_pre =0.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
#if defined USE_HEAT
!--- Set heat data ------------------------------------------------
    Pair = lin_interpol(time,Pair_data ,dt_air,nm_air)   !data from file
    Tair = lin_interpol(time,Tair_data ,dt_air,nm_air)   !data from file
    Qair = lin_interpol(time,Qair_data ,dt_air,nm_air)   !data from file
    rain = lin_interpol(time,rain_data ,dt_air,nm_air)   !data from file (mm h-1)
    rain = rain/3600000.0d0 ! mm h-1 -> m s-1
    cloud= lin_interpol(time,cloud_data,dt_air,nm_air)   !data from file
#endif
#if defined REEF_FLOW
    Hs = lin_interpol(time,Hs_data ,dt_wave,nm_wave)   !data from file
    Tp = lin_interpol(time,Tp_data, dt_wave,nm_wave)   !data from file
    
    Co(1,1,iTemp) = lin_interpol(time,offtemp_data, dt_offtemp,nm_offtemp)   !data from file
#endif
      
    return

  end subroutine setdata

  real(8) function lin_interpol(time,dataset,dt_data,datamax)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  linear interpolation between input data.            *
! *              time:     progress time (day)                         *
! *              dataset(i): data set                                  *
! *              dt_data: data interval (hour)                         *
! *                                                                    *
! **********************************************************************
!
    implicit none

    integer i
!    real lin_interpol
    real(8) dataset(datamax)
    real(8) time,dt_data
    integer datamax

    i=int(time*24.e0/dt_data)+1
    
    if(i.ge.datamax) i=datamax-1

    lin_interpol=dataset(i)+(dataset(i+1)-dataset(i))      &
          *(time*24.e0/dt_data-real(i-1))

    return

  end function lin_interpol
    
! **********************************************************************

  subroutine lin_interpol2(time,t_dat,d_dat,N_dat, out_dat, i_dat)
!
    implicit none

!    real lin_interpol
    real(8), intent(in ) :: time
    real(8), intent(in ) :: t_dat(N_dat)
    real(8), intent(in ) :: d_dat(N_dat)
    integer, intent(in ) :: N_dat
    real(8), intent(out) :: out_dat
    integer, intent(inout) :: i_dat

    real(8) :: dt
    integer ::  i,j


    do i=i_dat, N_dat-1
      if( t_dat(i+1)>=time .and. t_dat(i+1)>t_dat(i) ) exit
    end do
    if(i+1.gt.N_dat) then
      out_dat = d_dat(i)
    else
      dt = t_dat(i+1)-t_dat(i)
      out_dat = d_dat(i)+(d_dat(i+1)-d_dat(i)) * (time-t_dat(i))/dt
    end if
    i_dat = i

    return

  end subroutine lin_interpol2

  real(8) function solar_radi(time,PFD,start_time)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
! *              time:       progress time (day)                       *
! *              start_time: start time (day)                          *
! *                 ex.) 0:00-> start_time=0; 12:00-> start_time=0.5   *
! *                                                                    *
! **********************************************************************

    implicit none
    
    real(8), parameter :: pi = 3.141592654d0

    real(8) time, PFD,start_time

!    solar_radi=PFD/4.57*sin(2.*pi*(time+start_time-0.25))
    solar_radi=PFD/1.82d0*sin(2.*pi*(time+start_time-0.25))
                    !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)

    if(solar_radi.lt.0.) then
      solar_radi=0.
    endif

    return

  end function solar_radi


  real(8) function short_radi( time,start_time,start_day, lat   &
                              ,temp,Hum, mode )
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Caliculate short wave radiation (W/m2).             *
! *                 by Zillman equation                                *
! *                                                                    *
! *              time:       progress time (day)                       *
! *              start_time: start time (day)                          *
! *                 ex.) 0:00-> start_time=0; 12:00-> start_time=0.5   *
! *              start_day:  start day (day)                           *
! *                 ex.) Jan. 1-> start_day=1; Jun 22-> start_day=174  *
! *                                                                    *
! *              lat: latitude (degree)                                *
! *              Hum: rerative humidity (%)                            *
! *              temp: temperature (degree C)                          *
! *                                                                    *
! *              Hum: rerative humidity (%)                            *
! *                                                                    *
! *              mode: mode =1 -> day time progress                    *
! *                    mode =2 -> loop same day condition              *
! *                                                                    *
! **********************************************************************
!
    implicit none
    
    real(8), parameter :: pi = 3.141592654d0
    real(8), parameter :: Sc = 1353.0d0        !Solar constant (W m-2)

    real(8) :: time, start_time,lat,temp, Hum
    integer :: mode
    real(8) :: start_day
    real(8) :: yDay !day of year
    real(8) :: cosZ
    real(8) :: DA    !DA=2.*pi*yDay/365.
    real(8) :: e_vep !vapor pressure (Pa)
    real(8) :: dec   !declination (radian)
    real(8) :: HA    !hour angle (radian)
    real(8) :: lati
    integer :: i

    e_vep=611.*10.**(7.5*(temp/(temp+273.15-35.86))) * Hum/100.

    HA=(0.5-time)*2.*pi

    if(mode .eq. 1) then
      yDay=start_day+Int(time) ! day time progress
    else
      yDay=start_day            ! light condition releated same day.
    endif
    
    yDay=mod(yDay-1.,365.)+1.

    dec=23.44*cos((172.-yDay)*2.*pi/365.) *pi/180.
!    DA=2.*pi*yDay/365.
!    dec=0.006918-0.399912*cos(DA)   +0.070257*sin(DA)
!               -0.006758*cos(2.*DA)+0.000907*sin(2.*DA)
!               -0.002697*cos(3.*DA)+0.001480*sin(3.*DA)
!
    lati=lat *pi/180. !degree -> radian

    cosZ=sin(lati)*sin(dec)+cos(lati)*cos(dec)*cos(HA)



    if(cosZ.lt.0.) then
      short_radi=0.
    else
      short_radi=Sc*cosZ**2./((cosZ+2.7)*e_vep*1.e-5 + 1.085*cosZ +0.10)
    endif

    return

  end function short_radi


  real(8) function light_and_dark(time,PFD,interval,transperi)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
! *              time:       progress time (day)                       *
! *              interval: interval of light & dark (day)              *
! *              transperi: transition period between light & dark (day)*
! *                                                                    *
! **********************************************************************
!
    implicit none
    
    real(8), parameter :: pi = 3.141592654d0

    real(8) :: time, PFD,interval, transperi
    real(8) :: Imax   !light intensity (W/m2)

!          !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)
    Imax = PFD/4.57 ! for Al-Horani et al. (2003) experiment
!    Imax = 350./4.57 ! for Kuhl et al. (1995) experiment                      

!  Dark->light->Dark->light

!    if(mod(time,interval*2.).lt.interval) then
!      if(mod(time,interval*2.).lt.transperi) then
!        light_and_dark=
!               Imax/2.*(Cos(pi*mod(time,interval*2.)/transperi)+1.)
!      else
!        light_and_dark=0.
!      endif
!    else
!      if(mod(time,interval*2.).lt.interval+transperi) then
!        light_and_dark=
!              Imax/2.*(-Cos(pi*(mod(time,interval*2.)-interval)
!                                 /transperi)+1.)
!      else
!        light_and_dark= Imax
!      endif    
!    endif


!  light->Dark->light->Dark

    if(mod(time,interval*2.).lt.interval) then
      if(mod(time,interval*2.).lt.transperi) then
        light_and_dark=                                         &
              Imax/2.*(-Cos(pi*(mod(time,interval*2.)-interval) &
                                 /transperi)+1.)
      else
        light_and_dark= Imax
      endif
    else
      if(mod(time,interval*2.).lt.interval+transperi) then
        light_and_dark=                                            &
               Imax/2.*(Cos(pi*mod(time,interval*2.)/transperi)+1.)
      else
        light_and_dark=0.

      endif    
    endif

    return

  end function light_and_dark

  real(8) function light_and_dark2(time,PFD,interval1,interval2,transperi)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
! *              time:       progress time (day)                       *
! *              interval: interval of light & dark (day)              *
! *              transperi: transition period between light & dark (day)*
! *                                                                    *
! **********************************************************************
!
    implicit none
    
    real(8), parameter :: pi = 3.141592654d0

    real(8) :: time, PFD,interval1,interval2, transperi
    real(8) :: Imax   !light intensity (W/m2)

!            !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)
    Imax = PFD/2.1d0/(1.0d0-0.07d0) ! for Takahashi et al. (2004) experiment

!  Dark->light->Dark->light

    if(mod(time,interval1+interval2)<=interval2) then
      if(mod(time,interval1+interval2)<=transperi) then
        light_and_dark2=                                         &
              Imax/2.*(Cos(pi*(mod(time,interval1+interval2)-interval1) &
                                 /transperi)+1.)
      else
        light_and_dark2= 20.0d0/2.1d0/(1.0d0-0.07d0)
      endif
    else
      if(mod(time,interval1+interval2)<=interval1+transperi) then
        light_and_dark2=                                            &
               Imax/2.*(-Cos(pi*mod(time,interval1+interval2)/transperi)+1.)
      else
        light_and_dark2=Imax

      endif    
    endif

    return

  end function light_and_dark2

  real(8) function light_and_dark3(time,interval1,interval2,transperi)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
! *              time:       progress time (day)                       *
! *              interval: interval of light & dark (day)              *
! *              transperi: transition period between light & dark (day)*
! *                                                                    *
! **********************************************************************
!
    implicit none
    
    real(8), parameter :: pi = 3.141592654d0

    real(8) :: time, interval1,interval2, transperi
    real(8) :: Imax   !light intensity (W/m2)
    real(8) :: PFD

!            !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)
    if (time <= 1.1d0) then
      PFD = 0.0d0
    else if (time <= 2.1d0) then
      PFD = 0.0d0
    else if (time <= 3.1d0) then
      PFD = 100.0d0
    else if (time <= 4.1d0) then
      PFD = 250.0d0
    else if (time <= 5.1d0) then
      PFD = 500.0d0
    else
      PFD = 1000.0d0
    endif
    
    Imax = PFD/2.1d0/(1.0d0-0.07d0) ! for Takahashi et al. (2004) experiment

!  Dark->light->Dark->light

    if(mod(time,interval1+interval2)<=interval2) then
      if(mod(time,interval1+interval2)<=transperi) then
        light_and_dark3=                                         &
              Imax/2.*(Cos(pi*(mod(time,interval1+interval2)-interval1) &
                                 /transperi)+1.)
      else
        light_and_dark3= 20.0d0/2.1d0/(1.0d0-0.07d0)
      endif
    else
      if(mod(time,interval1+interval2)<=interval1+transperi) then
        light_and_dark3=                                            &
               Imax/2.*(-Cos(pi*mod(time,interval1+interval2)/transperi)+1.)
      else
        light_and_dark3=Imax

      endif    
    endif

    return

  end function light_and_dark3
      
end module mod_input

