      program CO2system_test
! **********************************************************************
! *                                                                    *
! *   Test program of CO2system.f                                      *
! *                                                                    *
! **********************************************************************
!
      USE mod_geochem
      
      implicit none

      real(8) T,S
      real(8) pH, d11BT, d11B_BOH4
      real(8) :: DOsatu, ssO2flux, U10, DO    
!
      integer i,j
!      integer istep
!
      open(52,file='./output/debug2.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(51,file='./output/result2.txt') 
      
      T=28.0d0+273.15
      S=34.5
      U10=1.0d0
      
      DO=300.0d0
      
!      do i=1,100
!        pH = 7.0d0+0.03*dble(i)
!        d11BT= 39.5d0
!        d11B_BOH4= d11B_BOH4_frompHd11BT(pH,d11BT, T, S)
!        write(51,*) pH, d11B_BOH4       
!      enddo

!      write(*,*) densSW(25.0d0,35.0d0),densW(25.0d0)

      DOsatu=O2satu(T,S)
      ssO2flux = Flux_O2(DO, DOsatu, U10, T, S )  ! sea to air is positive


      write(*,*) T, DOsatu, ssO2flux
      
      end program CO2system_test

!-----------------------------------------------------------------------


