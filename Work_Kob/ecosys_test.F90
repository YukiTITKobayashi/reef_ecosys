
!!!=== Copyright (c) 2012-2020 Takashi NAKAMURA  =====

#include "cppdefs.h"


PROGRAM ecosys_test
! **********************************************************************
! *                                                                    *
! *   Test program of mod_reef_ecosys                                  *
! *                                                                    *
! **********************************************************************
!
  USE mod_param
  USE mod_reef_ecosys
#if defined USE_HEAT
  USE mod_heat
#endif
  USE mod_input
  USE mod_output
  USE mod_geochem
#if defined REEF_FLOW
  USE mod_reef_flow
#endif

  implicit none
  
  integer, parameter :: Im = 1
  integer, parameter :: Jm = 1
  integer, parameter :: N  = 1

  real(8), parameter :: dt = 1.0d0
  integer, parameter :: isplitc  = 20
  integer, parameter :: isplitsed  = 1
  integer, parameter :: Ngrids = 1

  integer :: i,j,k,id, Nid
  integer :: istep, iprint
  integer :: nheat
  
  integer :: ipcl =1       ! Step of the protocol for setting 5 (Incubation chamber condition simulated Nakamura & Nakamori (2009) experiments)
  integer :: inight = 1 ! Timer for setting 6
  integer :: iflag = 1 ! Timer for setting 6
  integer :: iclose = 0 ! Timer for setting 6
  
#if defined CARBON_ISOTOPE
  real(8) :: R13C
#endif
!  For Output      
  real(8), parameter :: OUTPUT_INTERVAL = 5.0d0     ! Output interval (min)
  real(8), save :: dsec = 0.d0 !sec
  
  real(8) :: h1,h2
  real(8) :: fvol_rc, fvol_ch

  real(8) :: Tmax, time0
  integer :: nSetting

  namelist/ecosys_config/Tmax, time0
  namelist/ecosys_config/nSetting

  read (*, nml=ecosys_config)

!----- Open output files -------------------------

  CALL files_open
      
!----- Set initial conditions -------------------------

  CALL initialize_params(1, Im, 1, Jm, N, Nid)

!----- Import data -----------------------------------

  if (nSetting .eq. 4) then
    CALL read_timeseies
  else if (nSetting .eq. 5) then
    CALL read_chambercondition
  endif

!----- Set initial conditions -------------------------

#if defined USE_HEAT
  CALL initialize_heat(1, Ngrids, 1, Im, 1, Jm)
#endif
#if defined REEF_FLOW
  CALL initialize_reef_flow(1, Ngrids, 1, Im, 1, Jm)
#endif
  CALL initialize_reef_ecosys(1, Ngrids, 1, Im, 1, Jm)

!  call Coral_iniSizeDis
!  call Coral_Size2Cover

  time=time0
  
  istep=0
  iprint=0
      
!----- Write data labels -------------------------
  CALL write_env_lavel(10)        
#if defined CORAL_TESTMODE
  CALL write_crl_his_lavel(11)
  CALL write_crl_his_lavel(12)
  CALL write_crl_ave_lavel(21)
  CALL write_crl_ave_lavel(22)
# if defined CORAL_ZOOXANTHELLAE
  CALL write_zox_his_lavel(31)
  CALL write_zox_his_lavel(32)
#  if defined CORAL_PHOTOINHIBITION
  CALL write_zphot_his_lavel(41)
#  endif
# endif
#endif
#if defined ECOSYS_TESTMODE
  CALL write_ecosys_his_lavel(40)
#endif


!----- Main loop -------------------------------------------

  do istep=1, int(24.*60.*60./dt * Tmax)  +1

!------ Set environmental parameters ----------------------------

    CALL setdata(nSetting)
        
#if defined REEF_FLOW
    if (nSetting .eq. 4) then
        

!----- Reef hydrodynamics model ----------------------------------------

      h1 = REEF(1)%el(1,1)+REEF(1)%Dir(1,1)
    
      CALL reef_flow         &
!      input parameters
              (1, 1, 1        &   ! ng: nested grid number; i,j: position
              ,dt             &   ! Time step (sec)
              ,Hs             &   ! Significant wave hight at offshore (m)
              ,Tp             &   ! Significant Wave period (s)
              ,tide           &   ! Sea surface elevation at offshore (m)
               )
        
      h2 = (REEF(1)%el(1,1)+REEF(1)%Dir(1,1))
      dz(1,1,:) =h2/N
    end if
#endif
        
#if defined USE_HEAT

!----- Heat & mass balance model ----------------------------------------

    CALL heat_mass_balance    &
!      input parameters
              (1,1,1          &   ! ng: nested grid number; i,j: position
              ,N              &   ! Number of vertical grid (following ROMS vertical grid)
              ,dt             &   ! Time step (sec)
              ,dz(1,1,:)      &   ! dz(N): vertical grid size (m)
              ,swrad          &   ! Surface shortwave radiation (W m-2)
              ,Tair           &   ! air temperature (oC)
              ,Pair           &   ! atm pressure (hPa)
              ,Qair           &   ! vapor pressur (hPa)
              ,U10            &   ! wind speed (m s-1)
              ,rain           &   ! Precipitation volume flux (m s-1)
#ifdef LONGWAVE_IN
              ,dlwrad         &   ! Downward longwave radiation (W m-2)
#endif

              ,C(1,1,:,1,iTemp)     &   ! Tmp(N): Temperature (oC)
              ,C(1,1,:,1,iSalt)     &   ! Sal(N): Salinity (PSU)

!          output parameters
              ,dC_dt(1,1,:,iTemp)   &   ! dDIC_dt(N): dDIC/dt (umol kg-1 s-1)
               )
#endif

!----- Ecosystem model ----------------------------------------

    CALL reef_ecosys          &
!      input parameters
              (1, 1, 1        &   ! ng: nested grid number; i,j: position
              ,N              &   ! Number of vertical grid (following ROMS vertical grid)
              ,isplitc        &   ! Internal loop counts of coral polyp model
              ,isplitsed      &   ! Internal loop counts of sediment ecosystem model
              ,dt             &   ! Time step (sec)
              ,dz(1,1,:)      &   ! dz(N): vertical grid size (m)
              ,PFDsurf        &   ! Sea surface photon flux density (umol m-2 s-1)
              ,tau            &   ! bottom shear stress (N m-2)
              ,pCO2air        &   ! Air CO2 pertial pressure (uatm)
              ,U10            &   ! wind speed (m s-1)
#ifdef CORAL_POLYP
              ,p_coral(:,1,1) &   ! Coral coverage (0-1)
#endif
#ifdef SEAGRASS
              ,p_sgrass(1,1)  &   ! seagrass coverage (0-1)
#endif
#ifdef MACROALGAE
              ,p_algae(1,1)   &   ! algal coverage (0-1)
#endif
#ifdef SEDIMENT_ECOSYS
              ,p_sand(1,1)    &   ! sediment coverage (0-1)
#endif

              ,C(1,1,:,1,iTemp)     &   ! Tmp(N): Temperature (oC)
              ,C(1,1,:,1,iSalt)     &   ! Sal(N): Salinity (PSU)
              ,C(1,1,:,1,iTIC_)     &   ! DIC(N): Total dissolved inorganic carbon (DIC: umol kg-1)
              ,C(1,1,:,1,iTAlk)     &   ! TA (N): Total alkalinity (TA: umol kg-1)
              ,C(1,1,:,1,iOxyg)     &   ! DOx(N): Dissolved oxygen (umol L-1)
#if defined ORGANIC_MATTER
              ,C(1,1,:,1,iDOC(:))       &   ! DOC(N): Dissolved organic carbon (DOC: umol L-1)
              ,C(1,1,:,1,iPOC(:))       &   ! POC(N): Particulate organic carbon (POC: umol L-1)
              ,C(1,1,:,1,iPhyt(:))      &   ! PHY(N): phytoplankton (umol C L-1)
              ,C(1,1,:,1,iZoop(:))      &   ! ZOO(N): zooplankton (umol C L-1)
              ,C(1,1,:,1,iPIC(:))       &   ! PIC(N): Particulate inorganic carbon (PIC: umol L-1), coccolith (CaCO3)
# if defined FOODWEB_TESTMODE
              ,C(1,1,:,1,iDeadPhyt(:))  &   ! DeadPHY(Nphy): Dead phytoplankton (umol C L-1)
# endif
#endif
#if defined CARBON_ISOTOPE
              ,C(1,1,:,1,iT13C)         &   ! DI13C(N): 13C of DIC (umol kg-1)
# if defined ORGANIC_MATTER
              ,C(1,1,:,1,iDO13C(:))     &   ! DO13C (N): 13C of Labile Dissolved organic carbon (LDOC: umol L-1)
              ,C(1,1,:,1,iPO13C(:))     &   ! PO13C (N): 13C of Particulate organic carbon (POC: umol L-1)
              ,C(1,1,:,1,iPhyt13C(:))   &   ! PHY13C(N): 13C of phytoplankton1 (umol C L-1), dinoflagellate
              ,C(1,1,:,1,iZoop13C(:))   &   ! ZOO13C(N): 13C of zooplankton (umol C L-1)
              ,C(1,1,:,1,iPI13C(:))     &   ! PI13C (N): 13C of Particulate inorganic carbon (PIC: umol L-1), coccolith (CaCO3)
# endif
#endif
#if defined NUTRIENTS         
              ,C(1,1,:,1,iNO3_)     &   ! NO3(N): NO3 (umol L-1)
!              ,C(1,1,:,1,iNO2_)     &   ! NO2(N): NO2 (umol L-1)
              ,C(1,1,:,1,iNH4_)     &   ! NH4(N): NH4 (umol L-1)
              ,C(1,1,:,1,iPO4_)     &   ! PO4(N): PO4 (umol L-1)
# if defined ORGANIC_MATTER
              ,C(1,1,:,1,iDON(:))     &   ! DON(N): Dissolved organic nitrogen (DON: umol L-1)
              ,C(1,1,:,1,iPON(:))     &   ! PON(N): Particulate organic nitrogen (PON: umol L-1)
              ,C(1,1,:,1,iDOP(:))     &   ! DOP(N): Dissolved organic phosporius (DOP: umol L-1)
              ,C(1,1,:,1,iPOP(:))     &   ! POP(N): Particulate organic phosporius (POP: umol L-1)
# endif
# if defined NITROGEN_ISOTOPE
             ,C(1,1,:,1,i15NO3)     &   ! NO3_15N(N): 15N of NO3 (umol L-1)
!             ,C(1,1,:,1,i15NO2)     &   ! NO2_15N(N): 15N of NO2 (umol L-1)
             ,C(1,1,:,1,i15NH4)     &   ! NH4_15N(N): 15N of NH4 (umol L-1)
#  if defined ORGANIC_MATTER
             ,C(1,1,:,1,iDO15N(:))     &   ! DO15N (N): 15N of Labile Dissolved organic nitrogen (LDON: umol L-1)
             ,C(1,1,:,1,iPO15N(:))     &   ! PO15N (N): 15N of Particulate organic nitrogen (PON: umol L-1)
             ,C(1,1,:,1,iPhyt15N(:))   &   ! PHY15N(N): 15N of phytoplankton1 (umol C L-1), dinoflagellate
             ,C(1,1,:,1,iZoop15N(:))   &   ! ZOO15N(N): 15N of zooplankton (umol C L-1)
#  endif
# endif
#endif
#if defined COT_STARFISH         
              ,C(1,1,:,1,iCOTe)     &   ! COTe(N): COT starfish egg (umol L-1)
              ,C(1,1,:,1,iCOTl)     &   ! COTl(N): COT starfish larvae (umol L-1)
#endif
!      output parameters
              ,dC_dt(1,1,:,iTIC_)   &   ! dDIC_dt(N): dDIC/dt (umol kg-1 s-1)
              ,dC_dt(1,1,:,iTAlk)   &   ! dTA_dt (N): dTA/dt (umol kg-1 s-1)
              ,dC_dt(1,1,:,iOxyg)   &   ! dDOx_dt(N): dDO/dt (umol L-1 s-1)
#if defined ORGANIC_MATTER
              ,dC_dt(1,1,:,iDOC(1):iDOC(N_dom))      &   ! dDOC_dt(N): dDOC/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iPOC(1):iPOC(N_pom))      &   ! dPOC_dt(N): dPOC/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iPhyt(1):iPhyt(N_phyt))   &   ! dPHY_dt(N): dPHY/dt (umol C L-1 s-1)
              ,dC_dt(1,1,:,iZoop(1):iZoop(N_zoop))   &   ! dZOO_dt(N): dZOO/dt (umol C L-1 s-1)
              ,dC_dt(1,1,:,iPIC(1):iPIC(N_pim))      &   ! dPIC_dt(N): dPIC/dt (umol L-1 s-1)
# if defined FOODWEB_TESTMODE
              ,dC_dt(1,1,:,iDeadPhyt(1):iDeadPhyt(N_phyt))   &   ! dDeadPHY/dt  (umol L-1 s-1)
# endif
#endif
#if defined CARBON_ISOTOPE
              ,dC_dt(1,1,:,iT13C)   &   ! dDI13C_dt(N): dDI13C/dt (umol kg-1 s-1)
# if defined ORGANIC_MATTER     
              ,dC_dt(1,1,:,iDO13C(1):iDO13C(N_dom))    &   ! dDO13C_dt (N): dDO13C/dt  (umol L-1 s-1) 
              ,dC_dt(1,1,:,iPO13C(1):iPO13C(N_pom))    &   ! dPO13C_dt (N): dPO13C/dt  (umol L-1 s-1) 
              ,dC_dt(1,1,:,iPhyt13C(1):iPhyt13C(N_phyt))  &   ! dPHY13C_dt(N): dPHY13C/dt  (umol L-1 s-1)  
              ,dC_dt(1,1,:,iZoop13C(1):iZoop13C(N_zoop))  &   ! dZOO13C_dt(N): dZOO13C/dt  (umol L-1 s-1)  
              ,dC_dt(1,1,:,iPI13C(1):iPI13C(N_pim))    &   ! dPI13C_dt (N): dPI13C/dt  (umol L-1 s-1)  
# endif
#endif
#if defined NUTRIENTS 
              ,dC_dt(1,1,:,iNO3_)   &   ! dNO3_dt(N): dNO3/dt (umol L-1 s-1)
!              ,dC_dt(1,1,:,iNO2_)   &   ! dNO2_dt(N): dNO2/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iNH4_)   &   ! dNH4_dt(N): dNH4/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iPO4_)   &   ! dPO4_dt(N): dPO4/dt (umol L-1 s-1)
# if defined ORGANIC_MATTER
              ,dC_dt(1,1,:,iDON(1):iDON(N_dom))      &   ! dDON_dt(N): dDON/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iPON(1):iPON(N_pom))      &   ! dPON_dt(N): dPON/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iDOP(1):iDOP(N_dom))      &   ! dDOP_dt(N): dDOP/dt (umol L-1 s-1)
              ,dC_dt(1,1,:,iPOP(1):iPOP(N_pom))      &   ! dPOP_dt(N): dPOP/dt (umol L-1 s-1)
# endif
# if defined NITROGEN_ISOTOPE
              ,C_dt(1,1,:,i15NO3)        &   ! dNO3_15N_dt(N): dNO3_15N/dt (umol L-1 s-1)
!              ,C_dt(1,1,:,i15NO2)        &   ! dNO2_15N_dt(N): dNO2_15N/dt (umol L-1 s-1)
              ,C_dt(1,1,:,i15NH4)        &   ! dNH4_15N_dt(N): dNH4_15N/dt (umol L-1 s-1)
#  if defined ORGANIC_MATTER
              ,C_dt(1,1,:,iDO15N(1):iDO15N(N_dom))     &   ! dDO15N_dt (N): dDO15N/dt  (umol L-1 s-1) 
              ,C_dt(1,1,:,iPO15N(1):iPO15N(N_pom))     &   ! dPO15N_dt (N): dPO15N/dt  (umol L-1 s-1) 
              ,C_dt(1,1,:,iPhyt15N(1):iPhyt15N(N_phyt))   &   ! dPHY15N_dt(N): dPHY1_15N/dt  (umol L-1 s-1)  
              ,C_dt(1,1,:,iZoop15N(1):iZoop15N(N_zoop))   &   ! dZOO15N_dt(N): dZOO_15N/dt  (umol L-1 s-1)  
#  endif
# endif
#endif
#if defined COT_STARFISH         
              ,dC_dt(1,1,:,iCOTe)   &   ! dCOTe/dt(N): (umol L-1 s-1)
              ,dC_dt(1,1,:,iCOTl)   &   ! dCOTl/dt(N): (umol L-1 s-1)
#endif
              ,sspH           &   ! sea surface pH
              ,ssfCO2         &   ! sea surface fCO2 (uatm)
              ,ssWarg         &   ! sea surface aragonite saturation state
              ,ssCO2flux      &   ! sea surface CO2 flux (mmol m-2 s-1)
              ,ssO2flux       &   ! sea surface O2 flux (mmol m-2 s-1)
              ,PFDbott        &   ! Bottom photon flux density (umol m-2 s-1)
               )                    
!

    do k=1,N
        

!---------- for Stable condition ---------------------------------
              
      if (nSetting .eq. 1) then
          
        ! nothing to calculate
        
        if (time >= 5.0d0 ) then
!          C(1,1,k,1,iTemp) = 27.0d0
        end if

!---------- for Closed Chamber condition -------------------------

      else if (nSetting .eq. 2) then
        
        C(1,1,k,1,iTemp)=C(1,1,k,1,iTemp)+0.
        C(1,1,k,1,iSalt)=C(1,1,k,1,iSalt)+0.
        C(1,1,k,1,iSedi)=C(1,1,k,1,iSedi)+0.
        
        do id=4,Nid
          C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt
        end do

!---------- Constant Flow condition ------------------------------------

      else if (nSetting .eq. 3) then









!---------- Reef condition ------------------------------------

      else if (nSetting .eq. 4) then
#if defined REEF_FLOW
        fvol_rc = REEF(1)%Qrc(1,1)*REEF(1)%Wrc(1,1)
        fvol_ch = REEF(1)%Qch(1,1)*REEF(1)%Wch(1,1)
        do id=1,Nid
          C(1,1,k,1,id) = C(1,1,k,1,id)*h1/h2                          &
                        +(0.5d0*(ABS(fvol_rc)+fvol_rc)* Co(1,1,id)     &  !  (t unit) m s-1
                         -0.5d0*(ABS(fvol_rc)-fvol_rc)* C(1,1,k,1,id)  &  !  (t unit) m s-1
                         +0.5d0*(ABS(fvol_ch)+fvol_ch)* Co(1,1,id)     &  !  (t unit) m s-1
                         -0.5d0*(ABS(fvol_ch)-fvol_ch)* C(1,1,k,1,id)  &  !  (t unit) m s-1
                         )/h2/REEF(1)%Air(1,1)*dt
          C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt
        end do
#endif

!---------- Incubation chamber condition ------------------------------------

      else if (nSetting .eq. 5) then
        if( ipcl == 1) then
          C(1,1,k,1,iTIC_) = DIC_data(1)
          C(1,1,k,1,iTAlk) = TA_data(1)
          C(1,1,k,1,iOxyg) = DO_data(1)
!          C(1,1,k,1,iOxyg) = O2satu(C(1,1,k,1,iTemp)+273.15d0, C(1,1,k,1,iSalt))
#if defined CARBON_ISOTOPE
          R13C=R13C_fromd13C(0.7d0)
          C(1,1,k,1,iT13C) =R13C*C(1,1,k,1,iTIC_) !DI13C (umol kg-1) 
#endif
          ipcl = 2
        end if
      
        if (time >= 4.0d0 ) then
          do id=4,Nid
            C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt
          end do
          if (WQ_time(ipcl)-15.0d0/60.0d0 < time*24.0d0 .and. time*24.0d0 < WQ_time(ipcl)  ) then
            C(1,1,k,1,iTIC_) = DIC_data(ipcl)
            C(1,1,k,1,iTAlk) = TA_data(ipcl) 
            C(1,1,k,1,iOxyg) = DO_data(ipcl)
!            C(1,1,k,1,iOxyg) = O2satu(C(1,1,k,1,iTemp)+273.15d0, C(1,1,k,1,iSalt))
#if defined CARBON_ISOTOPE
            R13C=R13C_fromd13C(0.7d0)
            C(1,1,k,1,iT13C) =R13C*C(1,1,k,1,iTIC_) !DI13C (umol kg-1) 
#endif
          else if (time*24.0d0 >= WQ_time(ipcl)  ) then
            ipcl = ipcl +1
            if(ipcl>N_WQ) ipcl = N_WQ
          end if
        end if

!---------- Flume simulation ------------------------------------

      else if (nSetting .eq. 6) then
      
        if (aint(time)+6.0d0/24.0d0 <= time .and.  &
            time <= aint(time)+18.0d0/24.0d0     ) then
          if (inight == 1 ) then
            iflag = 1
            inight = 0
          end if
        else
          if (inight == 0 ) then
            iflag = 1
            inight = 1
          end if
        end if
        
        if (7.0d0+ 6.0d0/24.0d0 < time ) then
          if (iclose == 0 ) then
            iflag = 1
            iclose = 1
          end if
        end if
        if (7.0d0+12.0d0/24.0d0 < time ) then
          if (iclose == 1 ) then
            iflag = 1
            iclose = 2
          end if
        end if
        if (7.0d0+18.0d0/24.0d0 < time ) then
          if (iclose == 2 ) then
            iflag = 1
            iclose = 3
          end if
        end if
        if (8.0d0+ 0.0d0/24.0d0 < time ) then
          if (iclose == 3 ) then
            iflag = 1
            iclose = 4
          end if
        end if
        if (8.0d0+ 6.0d0/24.0d0 < time ) then
          if (iclose == 4 ) then
            iflag = 1
            iclose = 0
          end if
        end if
        
        if( inight == 1 .and. iflag == 1 ) then  ! Night treatment ~pH -0.1 (18:00-06:00)
#if defined FLUME_AMBIENT
          C(1,1,k,1,iTIC_) = 2090.0d0
#elif defined FLUME_HPCO2
          C(1,1,k,1,iTIC_) = 2269.0d0
#endif
          C(1,1,k,1,iTAlk) = 2340.0d0
          C(1,1,k,1,iOxyg) = O2satu(C(1,1,k,1,iTemp)+273.15d0, C(1,1,k,1,iSalt))
#if defined CARBON_ISOTOPE
          R13C=R13C_fromd13C(0.7d0)
          C(1,1,k,1,iT13C) =R13C*C(1,1,k,1,iTIC_) !DI13C (umol kg-1) 
#endif
          iflag = 0
        end if
        
        if( inight == 0 .and. iflag == 1 ) then
#if defined FLUME_AMBIENT
          C(1,1,k,1,iTIC_) = 2033.0d0
#elif defined FLUME_HPCO2
          C(1,1,k,1,iTIC_) = 2243.0d0
#endif
          C(1,1,k,1,iTAlk) = 2340.0d0
          C(1,1,k,1,iOxyg) = O2satu(C(1,1,k,1,iTemp)+273.15d0, C(1,1,k,1,iSalt))
#if defined CARBON_ISOTOPE
          R13C=R13C_fromd13C(0.7d0)
          C(1,1,k,1,iT13C) =R13C*C(1,1,k,1,iTIC_) !DI13C (umol kg-1) 
#endif
          iflag = 0
        end if
        
        
        if( iclose >= 1 ) then  ! Closed
!          tau = 1024*0.14*0.02**2. *0.5  !densSW*Cd*Ub**2  (2 cm s-1)
!          tau = 1024*0.14*0.05**2. *0.5  !densSW*Cd*Ub**2  (5 cm s-1)
          tau = 1024*0.14*0.10**2. *0.5  !densSW*Cd*Ub**2  (10 cm s-1)
          do id=4,Nid
            C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt
          end do
        end if

      end if
          
!------------------------------------------------------------------------------

!          depsed(i,j)=0.
!          radi(i,j,k)=-swrad

          
    enddo


!    p_coral(1,i,j)=p_coral(1,i,j)                         &
!          +(g_coral(1,i,j)-m_coral(1,i,j))*p_coral(1,i,j)
!
!    if (p_coral(1,i,j) .lt. 1.e-5) then
!      p_coral(1,i,j)=0.
!    end if


!    time=dt*istep/86400.
    time=time+dt/86400.

!------- Print section --------------------------------------

    if(time .ge. dsec/86400.) then
      dsec=dsec+OUTPUT_INTERVAL*60.
    
      write(*,*) 'Time (day): ', time  ! Output for standard out
      CALL write_env_data(10)

    endif

  enddo
      
!----- End loop --------------------------------------


!----- Close output files --------------------------------------

  CALL files_close

END PROGRAM ecosys_test
!----------------------------------------------------------------------!

!     End of main program

!-----------------------------------------------------------------------

