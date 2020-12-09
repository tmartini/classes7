MODULE modTH
IMPLICIT NONE

 public :: InitProcess_TH      
 public :: EvalXSec_PP_TH,EvalXSec_PP_TBH         ! hadronic
 public :: EvalAmp_QB_TH, EvalAmp_QbarBbar_TbarH  ! t-channel
 public :: EvalAmp_QQB_THBbar,EvalAmp_QQB_TbarHB  ! s-channel
  
 public :: EvalXSec_PP_TWMH,EvalXSec_PP_TBWPH, EvalXSec_PP_TWH     ! hadronic
 public :: EvalAmp_GB_TWMH,EvalAmp_GBB_TBWPH      ! tw-channel
 
 private

 real(8), parameter :: GeV=1d0/100d0 ! we are using units of 100GeV, i.e. Lambda=10 is 1TeV 
 real(8), parameter :: Gf = 1.16639d-5/GeV**2        ! Fermi constant
 real(8), parameter :: vev = 1.0d0/sqrt(Gf*sqrt(2.0d0))

 real(8) :: M_Reso=125d0*GeV
 

 CONTAINS 


 

SUBROUTINE InitProcess_TH(m_Higgs)
implicit none
real(8) :: m_Higgs


  M_Reso = m_Higgs
  
RETURN
END SUBROUTINE
  
  
SUBROUTINE EvalXSec_PP_TH(Mom,TTBHcoupl,TopDecays,Channel,Res)
implicit none
real(8) :: Mom(1:4,1:9),Res
complex(8) :: TTBHcoupl(1:2)
integer :: Channel! 0=s+t channel, 1=t channel, 2=s channel
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: eta1,eta2,Etot,Pztot,MatElSq_GG,MatElSq_QQB,MatElSq_QBQ
real(8) :: x1,x2,PDFScale,Collider_Energy,E_CMS
real(8) :: NNpdf(1:2,-6:7),m_ferm,LO_Res_Unpol(-6:6,-6:6)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'
      

      
      Collider_Energy = -(Mom(1,inLeft)+Mom(1,inRight))
      if( TopDecays.eq.0 ) then
          Etot = Mom(1,Hbos)+Mom(1,t)+Mom(1,qout)
          Pztot= Mom(4,Hbos)+Mom(4,t)+Mom(4,qout)
      else
          Etot = Mom(1,Hbos) + Mom(1,qout) + Mom(1,b)+Mom(1,lep)+Mom(1,nu)
          Pztot= Mom(4,Hbos) + Mom(4,qout) + Mom(4,b)+Mom(4,lep)+Mom(4,nu)
      endif
      x1 = (Etot+Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      x2 = (Etot-Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      E_CMS = dsqrt(x1*x2)*Collider_Energy
      PDFScale = 0.25d0*( m_Top + m_Reso ) * 100d0

      Mom(1:4,inLeft)  = x1 * Mom(1:4,inLeft)
      Mom(1:4,inRight) = x2 * Mom(1:4,inRight)
      
      call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
      call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))
      
      Res = 0d0
      if( Channel.eq.1 .or. Channel.eq.0 ) then
          call EvalAmp_QB_TH(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)      
          Res =   LO_Res_Unpol(Up_,Bot_)   * ( NNpdf(1,+2)*NNpdf(2,+5)  +  NNpdf(1,+4)*NNpdf(2,+5) )   &
                + LO_Res_Unpol(Bot_,Up_)   * ( NNpdf(1,+5)*NNpdf(2,+2)  +  NNpdf(1,+5)*NNpdf(2,+4) )   &
                + LO_Res_Unpol(ADn_,Bot_)  * ( NNpdf(1,-1)*NNpdf(2,+5)  +  NNpdf(1,-3)*NNpdf(2,+5) )   &
                + LO_Res_Unpol(Bot_,ADn_)  * ( NNpdf(1,+5)*NNpdf(2,-1)  +  NNpdf(1,+5)*NNpdf(2,-3) )
      endif
      if( Channel.eq.2 .or. Channel.eq.0 ) then
          call EvalAmp_QQB_THBBAR(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
          Res = Res &
                + LO_Res_Unpol(Up_,ADn_)   *  NNpdf(1,+2)*NNpdf(2,-1)     &
                + LO_Res_Unpol(ADn_,Up_)   *  NNpdf(1,-1)*NNpdf(2,+2)     &
                + LO_Res_Unpol(Chm_,AStr_) *  NNpdf(1,+4)*NNpdf(2,-3)     &
                + LO_Res_Unpol(AStr_,Chm_) *  NNpdf(1,-3)*NNpdf(2,+4)  
      endif
      Res = Res/x1/x2/(2d0*E_CMS**2)
            
!     restore incoming momenta (in all-outgoing convention)
      Mom(1,1:2) = -0.5d0*Collider_Energy
      Mom(4,1)   = -0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
      Mom(4,2)   = +0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
  
RETURN
END SUBROUTINE



      
SUBROUTINE EvalXSec_PP_TBH(Mom,TTBHcoupl,TopDecays,Channel,Res)
implicit none
real(8) :: Mom(1:4,1:13),Res
complex(8) :: TTBHcoupl(1:2)
integer :: Channel! 0=s+t channel, 1=t channel, 2=s channel
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: eta1,eta2,Etot,Pztot,MatElSq_GG,MatElSq_QQB,MatElSq_QBQ
real(8) :: x1,x2,PDFScale,Collider_Energy,E_CMS
real(8) :: NNpdf(1:2,-6:7),m_ferm,LO_Res_Unpol(-6:6,-6:6)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'
      

      
      Collider_Energy = -(Mom(1,inLeft)+Mom(1,inRight))
      if( TopDecays.eq.0 ) then
          Etot = Mom(1,Hbos)+Mom(1,t)+Mom(1,qout)
          Pztot= Mom(4,Hbos)+Mom(4,t)+Mom(4,qout)
      else
          Etot = Mom(1,Hbos) + Mom(1,qout) + Mom(1,b)+Mom(1,lep)+Mom(1,nu)
          Pztot= Mom(4,Hbos) + Mom(4,qout) + Mom(4,b)+Mom(4,lep)+Mom(4,nu)
      endif
      x1 = (Etot+Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      x2 = (Etot-Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      E_CMS = dsqrt(x1*x2)*Collider_Energy
      PDFScale = 0.25d0*( m_Top + m_Reso ) * 100d0

      Mom(1:4,inLeft)  = x1 * Mom(1:4,inLeft)
      Mom(1:4,inRight) = x2 * Mom(1:4,inRight)
      
      call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
      call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))

      Res = 0d0
      if( Channel.eq.1 .or. Channel.eq.0 ) then
          call EvalAmp_QbarBbar_TbarH(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
          Res = + LO_Res_Unpol(Dn_,ABot_)  * ( NNpdf(Dn_,1)*NNpdf(ABot_,2)  + NNpdf(Str_,1)*NNpdf(ABot_,2) )    &
                + LO_Res_Unpol(ABot_,Dn_)  * ( NNpdf(ABot_,1)*NNpdf(Dn_,2)  + NNpdf(ABot_,1)*NNpdf(Str_,2) )    &
                + LO_Res_Unpol(AUp_,ABot_) * ( NNpdf(AUp_,1)*NNpdf(ABot_,2) + NNpdf(AChm_,1)*NNpdf(ABot_,2))    &
                + LO_Res_Unpol(ABot_,AUp_) * ( NNpdf(ABot_,1)*NNpdf(AUp_,2) + NNpdf(ABot_,1)*NNpdf(AChm_,2))
      endif
      if( Channel.eq.2 .or. Channel.eq.0 ) then
          call EvalAmp_QQB_TBARHB(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
          Res = Res &
                + LO_Res_Unpol(AUp_,Dn_)   *  NNpdf(1,-2)*NNpdf(2,+1)     &
                + LO_Res_Unpol(Dn_,AUp_)   *  NNpdf(1,+1)*NNpdf(2,-2)     &
                + LO_Res_Unpol(AChm_,Str_) *  NNpdf(1,-4)*NNpdf(2,+3)     &
                + LO_Res_Unpol(Str_,AChm_) *  NNpdf(1,+3)*NNpdf(2,-4)  
      endif
      Res = Res/x1/x2/(2d0*E_CMS**2)

!     restore incoming momenta (in all-outgoing convention)
      Mom(1,1:2) = -0.5d0*Collider_Energy
      Mom(4,1)   = -0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
      Mom(4,2)   = +0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
  
RETURN
END SUBROUTINE


SUBROUTINE EvalXSec_PP_TWH(Mom,TTBHcoupl,TopDecays,Res)
implicit none
real(8) :: Mom(1:4,1:11),Res, Rest, Restb
complex(8) :: TTBHcoupl(1:2)
integer :: TopDecays! 0=stable, 1=di-leptonic

call EvalXSec_PP_TWMH(Mom,TTBHcoupl,TopDecays,Rest)
call EvalXSec_PP_TBWPH(Mom,TTBHcoupl,TopDecays,Restb)

Res=Rest+Restb
      
  
RETURN
END SUBROUTINE


SUBROUTINE EvalXSec_PP_TWMH(Mom,TTBHcoupl,TopDecays,Res)
implicit none
real(8) :: Mom(1:4,1:11),Res
complex(8) :: TTBHcoupl(1:2)
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: Etot,Pztot
real(8) :: x1,x2,PDFScale,Collider_Energy,E_CMS
real(8) :: NNpdf(1:2,-6:7),LO_Res_Unpol(-6:6,-6:6)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9, lepW=10, nuW=11
include 'includeVars.F90'
      

      
      Collider_Energy = -(Mom(1,inLeft)+Mom(1,inRight))
      if( TopDecays.eq.0 ) then
          Etot = Mom(1,Hbos)+Mom(1,t)+Mom(1,qout)
          Pztot= Mom(4,Hbos)+Mom(4,t)+Mom(4,qout)
      else
          Etot = Mom(1,Hbos) + Mom(1,lepW) + Mom(1,nuW) + Mom(1,b)+Mom(1,lep)+Mom(1,nu)
          Pztot= Mom(4,Hbos) + Mom(4,lepW) + Mom(4,nuW) + Mom(4,b)+Mom(4,lep)+Mom(4,nu)
      endif
      x1 = (Etot+Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      x2 = (Etot-Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      E_CMS = dsqrt(x1*x2)*Collider_Energy
      PDFScale = 0.25d0*( m_Top + m_Reso ) * 100d0

      Mom(1:4,inLeft)  = x1 * Mom(1:4,inLeft)
      Mom(1:4,inRight) = x2 * Mom(1:4,inRight)
      
      call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
      call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))
      
      Res = 0d0
      call EvalAmp_GB_TWMH(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
      Res = Res &
            + LO_Res_Unpol(0,Bot_)   *  NNpdf(1,0)*NNpdf(2,5)     &
            + LO_Res_Unpol(Bot_,0)   *  NNpdf(1,5)*NNpdf(2,0)

      Res = Res/x1/x2/(2d0*E_CMS**2)
            
!     restore incoming momenta (in all-outgoing convention)
      Mom(1,1:2) = -0.5d0*Collider_Energy
      Mom(4,1)   = -0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
      Mom(4,2)   = +0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
  
RETURN
END SUBROUTINE



      
SUBROUTINE EvalXSec_PP_TBWPH(Mom,TTBHcoupl,TopDecays,Res)
implicit none
real(8) :: Mom(1:4,1:11),Res
complex(8) :: TTBHcoupl(1:2)
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: Etot,Pztot
real(8) :: x1,x2,PDFScale,Collider_Energy,E_CMS
real(8) :: NNpdf(1:2,-6:7),LO_Res_Unpol(-6:6,-6:6)
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9, lepW=10, nuW=11
include 'includeVars.F90'
      

      
      Collider_Energy = -(Mom(1,inLeft)+Mom(1,inRight))
      if( TopDecays.eq.0 ) then
          Etot = Mom(1,Hbos)+Mom(1,t)+Mom(1,qout)
          Pztot= Mom(4,Hbos)+Mom(4,t)+Mom(4,qout)
      else
          Etot = Mom(1,Hbos) + Mom(1,lepW) + Mom(1,nuW) + Mom(1,b)+Mom(1,lep)+Mom(1,nu)
          Pztot= Mom(4,Hbos) + Mom(4,lepW) + Mom(4,nuW) + Mom(4,b)+Mom(4,lep)+Mom(4,nu)
      endif
      x1 = (Etot+Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      x2 = (Etot-Pztot * sign(1d0,Mom(4,inRight)) )/Collider_Energy
      E_CMS = dsqrt(x1*x2)*Collider_Energy
      PDFScale = 0.25d0*( m_Top + m_Reso ) * 100d0

      Mom(1:4,inLeft)  = x1 * Mom(1:4,inLeft)
      Mom(1:4,inRight) = x2 * Mom(1:4,inRight)
      
      call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
      call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))
      
      Res = 0d0
      call EvalAmp_GBB_TBWPH(Mom,TTBHcoupl,TopDecays,LO_Res_Unpol)
      Res = Res &
            + LO_Res_Unpol(0,ABot_)   *  NNpdf(1,0)*NNpdf(2,-5)     &
            + LO_Res_Unpol(ABot_,0)   *  NNpdf(1,-5)*NNpdf(2,0)

      Res = Res/x1/x2/(2d0*E_CMS**2)
            
!     restore incoming momenta (in all-outgoing convention)
      Mom(1,1:2) = -0.5d0*Collider_Energy
      Mom(4,1)   = -0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
      Mom(4,2)   = +0.5d0*Collider_Energy * sign(1d0,Mom(4,2))
  
RETURN
END SUBROUTINE



      
 
 
! t-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing
SUBROUTINE EvalAmp_QB_TH(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'



! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,qout)*MomExt(1,t)-MomExt(2,qout)*MomExt(2,t)-MomExt(3,qout)*MomExt(3,t)-MomExt(4,qout)*MomExt(4,t)
      p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
      p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)
      MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
      MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
      MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
      MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
      MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,qout)/2d0/p4Dp5
      MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
      MomExtFlat(7,1:4)=MomExt(1:4,qout)

    ! use different flattened momenta for top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
!          ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,b)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)
      ENDIF


! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp(:,:,:)=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call tdecay(5,6,8,9,10,za,zb,decay_amp)
      ENDIF
            
      call ubhtdamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Up_,Bot_,1:2))
      call ubhtdamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Bot_,Up_,1:2))        
      call ubhtdamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ADn_,Bot_,1:2))
      call ubhtdamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Bot_,ADn_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Up_,Bot_)  = cdabs(LOAmp(Up_,Bot_,1))**2  + cdabs(LOAmp(Up_,Bot_,2))**2
      LO_Res_Unpol(Bot_,Up_)  = cdabs(LOAmp(Bot_,Up_,1))**2  + cdabs(LOAmp(Bot_,Up_,2))**2
      LO_Res_Unpol(Chm_,Bot_) = LO_Res_Unpol(Up_,Bot_)
      LO_Res_Unpol(Bot_,Chm_) = LO_Res_Unpol(Bot_,Up_)

      LO_Res_Unpol(ADn_,Bot_) = cdabs(LOAmp(ADn_,Bot_,1))**2 + cdabs(LOAmp(ADn_,Bot_,2))**2
      LO_Res_Unpol(Bot_,ADn_) = cdabs(LOAmp(Bot_,ADn_,1))**2 + cdabs(LOAmp(Bot_,ADn_,2))**2
      LO_Res_Unpol(AStr_,Bot_)= LO_Res_Unpol(ADn_,Bot_)
      LO_Res_Unpol(Bot_,AStr_)= LO_Res_Unpol(Bot_,ADn_)      
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
! t-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing   
SUBROUTINE EvalAmp_QbarBbar_TbarH(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, qout=5, b=6,W=7,lep=8,nu=9
include 'includeVars.F90'



! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,qout)*MomExt(1,t)-MomExt(2,qout)*MomExt(2,t)-MomExt(3,qout)*MomExt(3,t)-MomExt(4,qout)*MomExt(4,t)
   p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
   p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)
   MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
   MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
   MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
   MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
   MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,qout)/2d0/p4Dp5
   MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
   MomExtFlat(7,1:4)=MomExt(1:4,qout)

    ! use different flattened momenta for anti-top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,b)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)

      ENDIF



! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call atdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
      call dbbarhtbaruamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Dn_,ABot_,1:2))
      call dbbarhtbaruamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ABot_,Dn_,1:2))
      call dbbarhtbaruamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,TTBHcoupl,LOAmp(AUp_,ABot_,1:2))
      call dbbarhtbaruamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ABot_,AUp_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Dn_,ABot_)  = cdabs(LOAmp(Dn_,ABot_,1))**2  + cdabs(LOAmp(Dn_,ABot_,2))**2
      LO_Res_Unpol(ABot_,Dn_)  = cdabs(LOAmp(ABot_,Dn_,1))**2  + cdabs(LOAmp(ABot_,Dn_,2))**2
      LO_Res_Unpol(Str_,ABot_) = LO_Res_Unpol(Dn_,ABot_)
      LO_Res_Unpol(ABot_,Str_) = LO_Res_Unpol(ABot_,Dn_)
      
      LO_Res_Unpol(AUp_,ABot_) = cdabs(LOAmp(AUp_,ABot_,1))**2 + cdabs(LOAmp(AUp_,ABot_,2))**2
      LO_Res_Unpol(ABot_,AUp_) = cdabs(LOAmp(ABot_,AUp_,1))**2 + cdabs(LOAmp(ABot_,AUp_,2))**2
      LO_Res_Unpol(AChm_,ABot_)= LO_Res_Unpol(AUp_,ABot_)
      LO_Res_Unpol(ABot_,AChm_)= LO_Res_Unpol(ABot_,AUp_)
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
   
   

! s-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing   
SUBROUTINE EvalAmp_QQB_THBBAR(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, bout=5, bdk=6,W=7,lep=8,nu=9
include 'includeVars.F90'


! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,bout)*MomExt(1,t)-MomExt(2,bout)*MomExt(2,t)-MomExt(3,bout)*MomExt(3,t)-MomExt(4,bout)*MomExt(4,t)
      p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
      p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)

      MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
      MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
      MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
      MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
      MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,bout)/2d0/p4Dp5
      MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
      MomExtFlat(7,1:4)=MomExt(1:4,bout)

    ! use different flattened momenta for top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
!          ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,bdk)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)
      ENDIF


! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp(:,:,:)=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call tdecay(5,6,8,9,10,za,zb,decay_amp)
      ENDIF
      call udbar_htbbaramp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Up_,ADn_,1:2))
      call udbar_htbbaramp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(ADn_,Up_,1:2))       
      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Up_,ADn_)  = cdabs(LOAmp(Up_,ADn_,1))**2  + cdabs(LOAmp(Up_,ADn_,2))**2
      LO_Res_Unpol(ADn_,Up_)  = cdabs(LOAmp(ADn_,Up_,1))**2  + cdabs(LOAmp(ADn_,Up_,2))**2
      LO_Res_Unpol(Chm_,AStr_)  = LO_Res_Unpol(Up_,ADn_) 
      LO_Res_Unpol(AStr_,Chm_)  = LO_Res_Unpol(ADn_,Up_) 
           
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2

   
RETURN
END SUBROUTINE EvalAmp_QQB_THBBAR




! s-channel 
! MARKUS: changed MomExt(1:4, inLeft/inRight) to all outgoing
SUBROUTINE EvalAmp_QQB_TBARHB(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:9),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol(-6:6,-6:6),s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),decay_amp(1:2),TTBHcoupl(1:2)
real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
integer :: j
integer, parameter :: inLeft=1,inRight=2,Hbos=3,t=4, bout=5, bdk=6,W=7,lep=8,nu=9
include 'includeVars.F90'


! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,bout)*MomExt(1,t)-MomExt(2,bout)*MomExt(2,t)-MomExt(3,bout)*MomExt(3,t)-MomExt(4,bout)*MomExt(4,t)
   p4Dp7=MomExt(1,lep)*MomExt(1,t)-MomExt(2,lep)*MomExt(2,t)-MomExt(3,lep)*MomExt(3,t)-MomExt(4,lep)*MomExt(4,t)
   p2Dp3=-MomExt(1,inRight)*MomExt(1,Hbos)+MomExt(2,inRight)*MomExt(2,Hbos)+MomExt(3,inRight)*MomExt(3,Hbos)+MomExt(4,inRight)*MomExt(4,Hbos)
   MomExtFlat(1,1:4)=-MomExt(1:4,inLeft)
   MomExtFlat(2,1:4)=-MomExt(1:4,inRight)
   MomExtFlat(3,1:4)=-m_Reso**2/2d0/p2Dp3*MomExt(1:4,inRight)
   MomExtFlat(4,1:4)=MomExt(1:4,Hbos)-MomExtFlat(3,1:4)
   MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,bout)/2d0/p4Dp5
   MomExtFlat(6,1:4)=MomExt(1:4,t)-MomExtFlat(5,1:4)
   MomExtFlat(7,1:4)=MomExt(1:4,bout)

    ! use different flattened momenta for anti-top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,lep)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,t)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,bdk)
         MomExtFlatDK(9,1:4)=MomExt(1:4,lep)
         MomExtFlatDK(10,1:4)=MomExt(1:4,nu)

      ENDIF



! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call atdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
      call ubard_Htbarbamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(AUp_,Dn_,1:2))
      call ubard_Htbarbamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,TTBHcoupl,LOAmp(Dn_,AUp_,1:2))


      ! coupling factors in decay incl in tdecay function
      LOAmp(:,:,:) = LOAmp(:,:,:) * 2d0*gwsq/vev*ci

      LO_Res_Unpol(Dn_,AUp_)  = cdabs(LOAmp(Dn_,AUp_,1))**2  + cdabs(LOAmp(Dn_,AUp_,2))**2
      LO_Res_Unpol(AUp_,Dn_)  = cdabs(LOAmp(AUp_,Dn_,1))**2  + cdabs(LOAmp(AUp_,Dn_,2))**2
      LO_Res_Unpol(Str_,AChm_)  = LO_Res_Unpol(Dn_,AUp_) 
      LO_Res_Unpol(AChm_,Str_)  = LO_Res_Unpol(AUp_,Dn_)       
      
      ColFac=9d0   
      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:) * ColFac * SpinAvg * QuarkColAvg**2
   
RETURN
END SUBROUTINE   
   
!tw channel   
!Till: all outgoing   
SUBROUTINE EvalAmp_GB_TWMH(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:11),LO_Res_UnPol(-6:6,-6:6)
complex(8) :: LOAmpab=0d0, LOAmpba=0d0

complex(8) :: UBbD(4),UBbD2(4),Eg(4),Eg2(4),UBtD(4),ECw(4)
complex(8) :: TTBHcoupl(1:2)
integer, parameter :: iL=1,iR=2, H=3, t=4, w=5, bt=6,Wt=7,lept=8,nut=9, lepW=10,nuW=11
integer :: i,j,bhel,thel,gpol,wpol,theld=-1,wpold=-1
! real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
include 'includeVars.F90'

!outgoing->incoming
      MomExt(:,iL:iR)=-MomExt(:,iL:iR)

      LO_Res_UnPol=0d0
      
      do bhel=-1,1,2
        call ubarSpi_Dirac(dcmplx(MomExt(:,iR)),0d0,bhel,UBbD) 
        call ubarSpi_Dirac(dcmplx(MomExt(:,iL)),0d0,bhel,UBbD2)  
        
        
        do gpol=-1,1,2  

          Eg=pol_mless(dcmplx(MomExt(:,iL)),gpol,.false.)
          Eg2=pol_mless(dcmplx(MomExt(:,iR)),gpol,.false.)
                          
                                
          if( TOPDECAYS .EQ. 0 ) then 
            theld=+1
            wpold=+1
          endif     
          do wpol=-1,wpold,1
            if( TOPDECAYS .EQ. 0 ) then 
              ECw=pol_mass(dcmplx(MomExt(:,w)),wpol,.true.)                
            elseif( TOPDECAYS .NE. 0 ) then
              call WDecay(Wm_,(/MomExt(1:4,lepW),MomExt(1:4,nuW)/),ECw)
            endif          
            do thel=-1,theld,2   
               
              if( TOPDECAYS .EQ. 0 ) then 
                call ubarSpi_Dirac(dcmplx(MomExt(:,t)),M_Top,thel,UBtD) 
              elseif( TOPDECAYS .NE. 0 ) then
                call TopDecay(Top_,(/MomExt(1:4,bt),MomExt(1:4,lept),MomExt(1:4,nut)/),UBtD)
              endif  
                  


              call gb_twmHamp(MomExt(1:4,iL),Eg,MomExt(1:4,iR),UBbD,MomExt(1:4,t),UBtD,MomExt(1:4,w),ECw,MomExt(1:4,H),TTBHcoupl,LOAmpab)
                        
              LO_Res_Unpol(0,6) = LO_Res_Unpol(0,6) + cdabs(LOAmpab)**2

                
                
              call gb_twmHamp(MomExt(1:4,iR),Eg2,MomExt(1:4,iL),UBbD2,MomExt(1:4,t),UBtD,MomExt(1:4,w),ECw,MomExt(1:4,H),TTBHcoupl,LOAmpba)


              LO_Res_Unpol(6,0) = LO_Res_Unpol(6,0) + cdabs(LOAmpba)**2
                
            enddo
          enddo 
        enddo
        

      enddo

        LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:)*1d0/2d0*1d0/2d0*1d0/8d0*(4d0/3d0) !* ColFac * SpinAvg * QuarkColAvg**2
      

RETURN
END SUBROUTINE EvalAmp_GB_TWMH

!tw channel   
!Till: all outgoing
SUBROUTINE EvalAmp_GBB_TBWPH(MomExt,TTBHcoupl,TopDecays,LO_Res_Unpol)
implicit none
integer :: TopDecays! 0=stable, 1=di-leptonic
real(8) :: MomExt(1:4,1:11),LO_Res_UnPol(-6:6,-6:6)
complex(8) :: LOAmpab=0d0, LOAmpba=0d0

complex(8) :: VbD(4),Eg(4),VbD2(4),Eg2(4),VtD(4),ECw(4)
complex(8) :: TTBHcoupl(1:2)
integer, parameter :: iL=1,iR=2, H=3, t=4, w=5, bt=6,Wt=7,lept=8,nut=9, lepW=10,nuW=11
integer :: i,j,bhel,thel,gpol,wpol,theld=-1,wpold=-1
! real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
include 'includeVars.F90'

!outgoing->incoming
      MomExt(:,iL:iR)=-MomExt(:,iL:iR)
      
      LO_Res_UnPol=0d0
      
      do bhel=-1,1,2
        call vSpi_Dirac(dcmplx(MomExt(:,iR)),0d0,bhel,VbD)       
        call vSpi_Dirac(dcmplx(MomExt(:,iL)),0d0,bhel,VbD2)  
                                                                           
        do gpol=-1,1,2  

          Eg=pol_mless(dcmplx(MomExt(:,iL)),gpol,.false.)
          Eg2=pol_mless(dcmplx(MomExt(:,iR)),gpol,.false.)
                        
                                 
          if( TOPDECAYS .EQ. 0 ) then  
            theld=+1
            wpold=+1
          endif    
          do wpol=-1,wpold,1         

            if( TOPDECAYS .EQ. 0 ) then  
              ECw=pol_mass(dcmplx(MomExt(:,w)),wpol,.true.) 
            elseif( TOPDECAYS .NE. 0 ) then
              call WDecay(Wp_,(/MomExt(1:4,lepW),MomExt(1:4,nuW)/),ECw)               
            endif  
                      
            do thel=-1,theld,2   
    
              if( TOPDECAYS .EQ. 0 ) then  
                call vSpi_Dirac(dcmplx(MomExt(:,t)),M_Top,thel,VtD) 
              elseif( TOPDECAYS .NE. 0 ) then
                call TopDecay(ATop_,(/MomExt(1:4,bt),MomExt(1:4,lept),MomExt(1:4,nut)/),VtD)
              endif  

              call gbb_tbwpHamp(MomExt(1:4,iL),Eg,MomExt(1:4,iR),VbD,MomExt(1:4,t),VtD,MomExt(1:4,w),ECw,MomExt(1:4,H),TTBHcoupl,LOampab)
              
              LO_Res_Unpol(0,-6) = LO_Res_Unpol(0,-6) + cdabs(LOAmpab)**2
                
                
              call gbb_tbwpHamp(MomExt(1:4,iR),Eg2,MomExt(1:4,iL),VbD2,MomExt(1:4,t),VtD,MomExt(1:4,w),ECw,MomExt(1:4,H),TTBHcoupl,LOampba)

              LO_Res_Unpol(-6,0) = LO_Res_Unpol(-6,0) + cdabs(LOAmpba)**2
                
            enddo
          enddo 
        enddo
      enddo

      LO_Res_Unpol(:,:) = LO_Res_Unpol(:,:)*1d0/2d0*1d0/2d0*1d0/8d0*(4d0/3d0) !* ColFac * SpinAvg * QuarkColAvg**2

RETURN
END SUBROUTINE EvalAmp_GBB_TBWPH   
   
   
   
! t-channel   
      SUBROUTINE ubhtdamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production u(p1)+b(p2)->H(p3)+t(p4)+d(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2),TTBHcoupl(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),KL,KR
        real(8)    :: s24,s34,s15,mt,mw
        include 'includeVars.F90'

        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top

! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))        
        
        ampw(1) = 1/( - mw**2 + s24)/( - mw**2 + s15)*za(p5,e4)*zb(p1,&
     & p2)*mw**2 + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)/(zb(k4,e4)&
     & )*za(p5,e3)*zb(p2,k4)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + s24)&
     & /( - mw**2 + s15)/(zb(k4,e4))*za(p5,k3)*zb(p2,k4)*zb(k3,p1)*&
     & mt**2
        
        ampw(2) = 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)*za(p5,e3)*&
     & zb(p2,e4)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + &
     & s15)*za(p5,k3)*zb(p2,e4)*zb(k3,p1)*mt + 1/( - mw**2 + s24)/( - &
     & mw**2 + s15)/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*mw**2*mt
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p5,e4&
     & )*zb(p1,p2)*mt*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15&
     & )/(zb(k4,e4))*za(p5,p1)*zb(p1,p2)*zb(p1,k4)*mt*vev*KL - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s15)/(zb(k4,e4))*za(p5,p2)*zb(p1,&
     & p2)*zb(p2,k4)*mt*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p5,p1&
     & )*zb(p1,p2)*zb(p1,e4)*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s15)*za(p5,p2)*zb(p1,p2)*zb(p2,e4)*vev*KL - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s15)/(za(k4,e4))*za(p5,k4)*zb(p1,p2)*&
     & mt**2*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      RETURN
      END SUBROUTINE
        
    
! t-channel    
      SUBROUTINE dbbarhtbaruamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production d(p1)+bbar(p2)->H(p3)+tbar(p4)+u(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
! modification of amplitude in MCFM and hep-ph:/1302.3856
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),TTBHcoupl(1:2),KL,KR
        real(8)    :: s24,s34,s15,mt,mw
        include 'includeVars.F90'

        s24=s(p2,k4)+s(p2,e4)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s15=s(p1,p5)
        mw=M_W
        mt=m_Top

! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))
        
        ampw(1) = 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)*za(p2,e4)*&
     & za(p5,e3)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + &
     & s15)*za(p2,e4)*za(p5,k3)*zb(k3,p1)*mt - 1/( - mw**2 + s24)/( - &
     & mw**2 + s15)/(zb(k4,e4))*za(p2,p5)*zb(p1,k4)*mw**2*mt
        
        ampw(2) =  - 1/( - mw**2 + s24)/( - mw**2 + s15)*za(p2,p5)*zb(&
     & p1,e4)*mw**2 + 1d0/2d0/( - mw**2 + s24)/( - mw**2 + s15)/(za(k4,&
     & e4))*za(p2,k4)*za(p5,e3)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + &
     & s24)/( - mw**2 + s15)/(za(k4,e4))*za(p2,k4)*za(p5,k3)*zb(k3,p1)*&
     & mt**2
        
        ampt(1) = 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p2,p5)*&
     & za(e4,p2)*zb(p2,p1)*vev*KR + 1d0/2d0/( - mt**2 + s34)/( - mw**2&
     &  + s15)*za(p2,p5)*za(e4,p5)*zb(p5,p1)*vev*KR + 1d0/2d0/( - mt**2&
     &  + s34)/( - mw**2 + s15)/(zb(k4,e4))*za(p2,p5)*zb(p1,k4)*mt**2*&
     & vev*KL
        
        ampt(2) = 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)*za(p2,p5)*&
     & zb(p1,e4)*mt*vev*KL + 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s15)&
     & /(za(k4,e4))*za(p2,p5)*za(k4,p2)*zb(p2,p1)*mt*vev*KR + 1d0/2d0/(&
     &  - mt**2 + s34)/( - mw**2 + s15)/(za(k4,e4))*za(p2,p5)*za(k4,p5)&
     & *zb(p5,p1)*mt*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)

      RETURN        
      END SUBROUTINE 
    


    
! s-channel    
      SUBROUTINE udbar_htbbaramp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production u(p1)+dbar(p2)->H(p3)+t(p4)+bbar(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),TTBHcoupl(1:2),KL,KR
        real(8)    :: s12,s34,s45,mt,mw
        include 'includeVars.F90'    


        s45=s(k4,p5)+s(e4,p5)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s12=s(p1,p2)
        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))

        
        ampw(1) = 1/( - mw**2 + s45)/( - mw**2 + s12)*za(p2,e4)*zb(p1,&
     & p5)*mw**2 + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)/(zb(k4,e4)&
     & )*za(p2,e3)*zb(p5,k4)*zb(e3,p1)*mt**2 + 1d0/2d0/( - mw**2 + s45)&
     & /( - mw**2 + s12)/(zb(k4,e4))*za(p2,k3)*zb(p5,k4)*zb(k3,p1)*&
     & mt**2
        
        ampw(2) = 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)*za(p2,e3)*&
     & zb(p5,e4)*zb(e3,p1)*mt + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + &
     & s12)*za(p2,k3)*zb(p5,e4)*zb(k3,p1)*mt + 1/( - mw**2 + s45)/( - &
     & mw**2 + s12)/(za(k4,e4))*za(p2,k4)*zb(p1,p5)*mw**2*mt
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p2,e4&
     & )*zb(p1,p5)*mt*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12&
     & )/(zb(k4,e4))*za(p2,p1)*zb(p1,p5)*zb(p1,k4)*mt*vev*KL - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s12)/(zb(k4,e4))*za(p2,p5)*zb(p1,&
     & p5)*zb(p5,k4)*mt*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p2,p1&
     & )*zb(p1,p5)*zb(p1,e4)*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s12)*za(p2,p5)*zb(p1,p5)*zb(p5,e4)*vev*KL - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s12)/(za(k4,e4))*za(p2,k4)*zb(p1,p5)*&
     & mt**2*vev*KR
        
        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      RETURN  
      END SUBROUTINE
        


! s-channel    
      SUBROUTINE ubard_Htbarbamp(p1,p2,e3,k3,k4,e4,p5,za,zb,s,mdecay,TTBHcoupl,amp)
! amplitude for production ubar(p1)+d(p2)->H(p3)+tbar(p4)+b(p5)
! allowing for scalar & pseudoscalar couplings of Higgs to top
        implicit none
        integer    :: p1,p2,p3,e3,k3,e4,p5,k4
        complex(8) :: za(:,:),zb(:,:),mdecay(1:2)
        real(8)    :: s(:,:)
        complex(8) :: amp(2),ampw(2),ampt(2),TTBHcoupl(1:2),KL,KR
        real(8)    :: s45,s34,s12,mt,mw
        include 'includeVars.F90'    


        s45=s(k4,p5)+s(e4,p5)+s(e4,k4)
        s34=s(p1,p2)+s(p1,p5)+s(p2,p5)
        s12=s(p1,p2)

        mw=M_W
        mt=m_Top
        
! there is a factor of -2 relative to ttbH
        KL=-mt/vev*(TTBHcoupl(1)-(0d0,1d0)*TTBHcoupl(2))
        KR=-mt/vev*(TTBHcoupl(1)+(0d0,1d0)*TTBHcoupl(2))

       
        ampw(1) = 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)*za(p1,e3)*&
     & za(p5,e4)*zb(e3,p2)*mt + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + &
     & s12)*za(p1,k3)*za(p5,e4)*zb(k3,p2)*mt + 1/( - mw**2 + s45)/( - &
     & mw**2 + s12)/(zb(k4,e4))*za(p1,p5)*zb(p2,k4)*mw**2*mt
        
        ampw(2) = 1/( - mw**2 + s45)/( - mw**2 + s12)*za(p1,p5)*zb(p2,&
     & e4)*mw**2 + 1d0/2d0/( - mw**2 + s45)/( - mw**2 + s12)/(za(k4,e4)&
     & )*za(p1,e3)*za(p5,k4)*zb(e3,p2)*mt**2 + 1d0/2d0/( - mw**2 + s45)&
     & /( - mw**2 + s12)/(za(k4,e4))*za(p1,k3)*za(p5,k4)*zb(k3,p2)*&
     & mt**2
        
        ampt(1) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p1,p5&
     & )*za(e4,p1)*zb(p1,p2)*vev*KR - 1d0/2d0/( - mt**2 + s34)/( - &
     & mw**2 + s12)*za(p1,p5)*za(e4,p5)*zb(p5,p2)*vev*KR - 1d0/2d0/( - &
     & mt**2 + s34)/( - mw**2 + s12)/(zb(k4,e4))*za(p1,p5)*zb(p2,k4)*&
     & mt**2*vev*KL
        
        ampt(2) =  - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12)*za(p1,p5&
     & )*zb(p2,e4)*mt*vev*KL - 1d0/2d0/( - mt**2 + s34)/( - mw**2 + s12&
     & )/(za(k4,e4))*za(p1,p5)*za(k4,p1)*zb(p1,p2)*mt*vev*KR - 1d0/2d0&
     & /( - mt**2 + s34)/( - mw**2 + s12)/(za(k4,e4))*za(p1,p5)*za(k4,&
     & p5)*zb(p5,p2)*mt*vev*KR                

        amp(1)=(ampw(1)+ampt(1))*mdecay(1)
        amp(2)=(ampw(2)+ampt(2))*mdecay(2)
        
      RETURN  
      END SUBROUTINE 
    
    
    
    
    
    
    
    SUBROUTINE TDECAY(k4,e4,b,ep,nu,za,zb,dkamp)
! top decay routine, taken from MCFM, see hep-ph:/1204.1513
       implicit none
       integer :: k4,e4,b,ep,nu
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp
       include 'includeVars.F90'
 
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = za(b,nu)*zb(ep,e4)
       dkamp(2) = m_top * za(b,nu)*zb(ep,k4)/za(e4,k4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * (4d0*dsqrt(2d0)*m_W**2*GF)        
   
     END SUBROUTINE
   
    
    SUBROUTINE ATDECAY(k4,e4,bbar,em,nubar,za,zb,dkamp)
! anti-top decay routine, taken from MCFM, see hep-ph:/1204.1513
       implicit none
       integer :: k4,e4,bbar,em,nubar
       complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
       real(8) :: NWAFactor_Top,NWAFactor_W
       complex(8) :: WProp
       include 'includeVars.F90'
  
   ! if one flattens the top wrt to e, then amp(2) = 0
       dkamp(1) = -m_top * zb(bbar,nubar)*za(em,k4)/zb(e4,k4)
       dkamp(2) = -zb(bbar,nubar)*za(em,e4)
   
       NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top) 
       NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
       WProp = (0d0,-1d0)*NWAFactor_W
   
       dkamp = dkamp * WProp * NWAFactor_Top * (4d0*dsqrt(2d0)*m_W**2*GF)
   
   
     END SUBROUTINE 
   

 


    SUBROUTINE convert_to_MCFM(p,pout)
      implicit none
! converts from (E,px,py,pz) to (px,py,pz,E)
      real(8) :: p(1:4),tmp(1:4)
      real(8), optional :: pout(1:4)

      if( present(pout) ) then
          pout(1)=p(2)  
          pout(2)=p(3)  
          pout(3)=p(4) 
          pout(4)=p(1)  
      else
          tmp(1)=p(1)
          tmp(2)=p(2)
          tmp(3)=p(3)
          tmp(4)=p(4)

          p(1)=tmp(2)  
          p(2)=tmp(3) 
          p(3)=tmp(4)  
          p(4)=tmp(1)  
      endif  
      
    END SUBROUTINE




subroutine spinoru(N,p,za,zb,s)
!---Calculate spinor products      
!---taken from MCFM & modified by R. Rontsch, May 2015
!---extended to deal with negative energies ie with all momenta outgoing                                                                
!---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,                                                                                  
!---za(i,j)*zb(j,i)=s(i,j)                      
      implicit none
      real(8) :: p(:,:),two
      integer, parameter :: mxpart=14
      complex(8):: c23(N),f(N),rt(N),za(:,:),zb(:,:),czero,cone,ci
      real(8)   :: s(:,:)
      integer i,j,N
      
      if (size(p,1) .ne. N) then
         print *, "spinorz: momentum mismatch"
         stop
      endif
      two=2d0
      czero=dcmplx(0d0,0d0)
      cone=dcmplx(1d0,0d0)
      ci=dcmplx(0d0,1d0)
      

!---if one of the vectors happens to be zero this routine fails.                                                                                                                
      do j=1,N
         za(j,j)=czero
         zb(j,j)=za(j,j)

!-----positive energy case                                                                                                                                                      
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(p(j,4)+p(j,1))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
!-----negative energy case                                                                                                                                                      
            rt(j)=dsqrt(-p(j,4)-p(j,1))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=ci
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

         if (abs(s(i,j)).lt.1d-5) then
         zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
         else
         zb(i,j)=-dcmplx(s(i,j))/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

    end subroutine spinoru



      subroutine gb_twmHamp(Pg,Eg,Pb,UBb,Pt,UBt,Pw,ECw,Ph,TTBHcoupl,amp)
! amplitude for production g(Pg)+b(Pb)->t(Pt)+wm(Pw)+h(Ph)
! allowing for scalar & pseudoscalar couplings of Higgs to top (Spinors and Pol-Vectors in Dirac representation!)
        implicit none
        real(8)    :: Pg(4),Pb(4),Pt(4),Pw(4),Ph(4)
        complex(8) :: UBb(4),Eg(4),UBt(4),ECw(4)
        complex(8) :: amp,ampts,ampws,amptt1,amptt2,ampwt
        real(8)    :: mt,mw,v, gs
        complex(8)    :: kap,kapt,a1WW,a2WW,a4WW   
        complex(8) TTBHcoupl(1:2)
        include 'includeVars.F90'            


        mw=M_W
        mt=m_Top
        v=vev
        
        gs=dsqrt(4d0*Pi*alphas);
        
        kap=TTBHcoupl(1)
        kapt=TTBHcoupl(2)
        a1WW=(2d0,0d0)
        a2WW=(0d0,0d0)
        a4WW=(0d0,0d0)
        

        ampts = -((gs*mw*(-(DCONJG(UBb(2))*(Eg(1)*((Pb(1) + Pg(1))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   Eg(4)*(-((Pb(1) + Pg(1))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Pb(1) + Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))) -&
     &   DCONJG(UBb(4))*(Eg(1)*((Pb(1) + Pg(1))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   Eg(4)*(-((Pb(1) + Pg(1))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (Pb(1) + Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))))) -&
     &   DCONJG(UBb(1))*(-((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(-((Pb(1) + Pg(1))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   + Eg(1)*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + Eg(4)*(-((Pb(2) + Pg(2)&
     &   + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))))) -&
     &   DCONJG(UBb(3))*(-((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(-((Pb(1) + Pg(1))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(1) -&
     &   ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) +&
     &   dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(1) + ECw(4))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   + Eg(1)*(-((Pb(2) + Pg(2) + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + Eg(4)*(-((Pb(2) + Pg(2)&
     &   + dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*(mt*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(ECw(1) - ECw(4))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + mt*(ECw(1) - ECw(4))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (Pb(1) +&
     &   Pg(1))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) +&
     &   Pg(4))*(mt*(ECw(1) + ECw(4))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + mt*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (ECw(1) + ECw(4))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))))))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) + Pg(2))**2 +&
     &   (Pb(3) + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mt**2 - (Ph(1) + Pt(1))**2 + (Ph(2) +&
     &   Pt(2))**2 + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)))

        ampws = (dcmplx(0d0,-1d0)*gs*mw*(DCONJG(UBb(1))*(Eg(1)*((Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) -&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)))) +&
     &   DCONJG(UBb(3))*(Eg(1)*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) +&
     &   Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1))&
     &   + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) -&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   Pg(3)))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)))) +&
     &   DCONJG(UBb(2))*((Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1)&
     &   + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(1)*(-((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4))) + (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) + (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)))) +&
     &   DCONJG(UBb(4))*((Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pb(1) + Pg(1))*((a1WW*mw**2*(ECw(1)&
     &   + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(4) +&
     &   Pg(4))*((a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) +&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(2) + Pg(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(4)*((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4)) - (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) - (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4))) +&
     &   Eg(1)*(-((Pb(2) + Pg(2) - dcmplx(0d0,1d0)*(Pb(3) + Pg(3)))*((a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(1) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(2) + (a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*UBt(3) + (a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) + dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(4))) + (Pb(1) +&
     &   Pg(1))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(4)) + (Pb(4)&
     &   + Pg(4))*((a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(1) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))*UBt(2) +&
     &   (a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))))*UBt(3) + (a1WW*mw**2*(ECw(1) -&
     &   ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*UBt(4))))))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) + Pg(2))**2 +&
     &   (Pb(3) + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mw**2 - (Ph(1) + Pw(1))**2 + (Ph(2) +&
     &   Pw(2))**2 + (Ph(3) + Pw(3))**2 + (Ph(4) + Pw(4))**2))

        amptt1 = (dcmplx(0d0,-1d0)*gs*mt*mw*(DCONJG(UBb(1))*((ECw(1) + ECw(4))*((Pb(4) -&
     &   Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(1) + ECw(4))*((mt - Pb(1) +&
     &   Pw(1))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) + Pw(1))*(-(kapt*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kapt*(Pg(4) - Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) -&
     &   Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3)&
     &   - Eg(4)*UBt(4)) - dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4)))))) + DCONJG(UBb(3))*((ECw(1) +&
     &   ECw(4))*((Pb(4) - Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(1) + ECw(4))*((mt - Pb(1) +&
     &   Pw(1))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) + dcmplx(0d0,1d0)*kapt*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4) -&
     &   Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) + Pw(1))*(-(kapt*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kapt*(Pg(4) - Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) -&
     &   Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3)&
     &   - Eg(4)*UBt(4)) - dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4)))))) + DCONJG(UBb(2))*((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(4) - Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt +&
     &   Pb(1) - Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(1) - ECw(4))*((Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(1) -&
     &   ECw(4))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt +&
     &   Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) +&
     &   Pw(1))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4)))))) + DCONJG(UBb(4))*((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(4) - Pw(4))*(-(kap*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt +&
     &   Pb(1) - Pw(1))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)))) + (ECw(1) - ECw(4))*((Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - dcmplx(0d0,1d0)*(mt + Pb(1) -&
     &   Pw(1))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(-(kap*(mt + Pg(1) -&
     &   Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3))) +&
     &   dcmplx(0d0,1d0)*kapt*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kapt*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kapt*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + kap*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) +&
     &   kap*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(4)&
     &   - Pw(4))*(kapt*(mt + Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) -&
     &   Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) - dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt&
     &   - Pg(1) + Pt(1))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) - Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + dcmplx(0d0,1d0)*(Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(kapt*((mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) - (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) +&
     &   (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3)&
     &   - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) + kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3)&
     &   - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))) + dcmplx(0d0,1d0)*(ECw(1) -&
     &   ECw(4))*((Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(kapt*(mt +&
     &   Pg(1) - Pt(1))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) + dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) + kap*(dcmplx(0d0,1d0)*Pg(2) - Pg(3) -&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) -&
     &   Eg(4)*UBt(2) + Eg(1)*UBt(4)) + dcmplx(0d0,1d0)*kap*(mt - Pg(1) + Pt(1))*(Eg(1)*UBt(1)&
     &   + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(Pg(4) -&
     &   Pt(4))*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) -&
     &   kapt*(Pg(2) + dcmplx(0d0,1d0)*(Pg(3) + dcmplx(0d0,1d0)*Pt(2) - Pt(3)))*(Eg(1)*UBt(2) +&
     &   Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (mt - Pb(1) +&
     &   Pw(1))*(-(kapt*(Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3))) - dcmplx(0d0,1d0)*kap*(mt + Pg(1) - Pt(1))*(Eg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) + kapt*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   kap*(dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) - kapt*(mt - Pg(1) +&
     &   Pt(1))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4)) -&
     &   dcmplx(0d0,1d0)*kap*(Pg(4) - Pt(4))*(Eg(1)*UBt(2) + Eg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) - (Pb(4) - Pw(4))*(kapt*((mt + Pg(1) -&
     &   Pt(1))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) -&
     &   (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*UBt(1) +&
     &   Eg(4)*UBt(3) + (Eg(2) + dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + (Pg(4) - Pt(4))*(Eg(1)*UBt(2)&
     &   + Eg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*UBt(3) - Eg(4)*UBt(4))) +&
     &   kap*((dcmplx(0d0,1d0)*Pg(2) + Pg(3) - dcmplx(0d0,1d0)*Pt(2) - Pt(3))*(Eg(4)*UBt(1) +&
     &   Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) + Eg(1)*UBt(3)) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) + Eg(1)*UBt(4)) +&
     &   (mt - Pg(1) + Pt(1))*(dcmplx(0d0,1d0)*Eg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(2)*UBt(3) +&
     &   Eg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*UBt(4))))))))/(sqrt2*v**2*(mt**2 - (Pg(1) -&
     &   Pt(1))**2 + (Pg(2) - Pt(2))**2 + (Pg(3) - Pt(3))**2 + (Pg(4) - Pt(4))**2)*(mt**2 -&
     &   (Pb(1) - Pw(1))**2 + (Pb(2) - Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))
     
     amptt2 = (gs*mw*(-(DCONJG(UBb(2))*((ECw(1) - ECw(4))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2)&
     &   - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) - ECw(4))*((mt -&
     &   Pb(1) + Pw(1))*(mt*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) +&
     &   mt*(Ph(4) + Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))))))&
     &   - DCONJG(UBb(4))*((ECw(1) - ECw(4))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2)&
     &   - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) - ECw(4))*((mt -&
     &   Pb(1) + Pw(1))*(mt*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) +&
     &   mt*(Ph(4) + Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) - (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(2) - dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) - Pw(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   - DCONJG(UBb(1))*((ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) +&
     &   ECw(4))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) - (Pb(4) - Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3)&
     &   + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(1) + ECw(4))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2)&
     &   + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) +&
     &   (mt - Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))))&
     &   - DCONJG(UBb(3))*((ECw(2) + dcmplx(0d0,1d0)*ECw(3))*(-((Pb(4) - Pw(4))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))) + (Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(1) +&
     &   ECw(4))*((Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) + (Pb(4) - Pw(4))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) -&
     &   (mt + Pb(1) - Pw(1))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) + (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*((mt - Pb(1) + Pw(1))*(mt*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) +&
     &   (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*Eg(1)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1)&
     &   + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + Eg(4)*(-(mt*(Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3))) + mt*(Ph(4) + Pt(4))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4)))) - (Pb(4) - Pw(4))*(mt*Eg(1)*((Ph(2) + Pt(2) - dcmplx(0d0,1d0)*(Ph(3)&
     &   + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) + Pt(4))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt - Ph(1) - Pt(1))*(kapt*UBt(2)&
     &   - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) + Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) +&
     &   kapt*UBt(4))) + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1)&
     &   - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(mt*Eg(1)*((Ph(4) +&
     &   Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4))))) +&
     &   (ECw(1) + ECw(4))*((mt - Pb(1) + Pw(1))*(mt*Eg(4)*((Ph(4) + Pt(4))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (mt + Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) +&
     &   kapt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) -&
     &   dcmplx(0d0,1d0)*kap*UBt(4))) + mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(1)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) +&
     &   (Pb(2) - Pw(2) + dcmplx(0d0,1d0)*(Pb(3) - Pw(3)))*(mt*Eg(1)*((Ph(2) + Pt(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) - (Ph(4) +&
     &   Pt(4))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (mt + Ph(1) +&
     &   Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) - mt*Eg(4)*((dcmplx(0d0,1d0)*Ph(2)&
     &   + Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) +&
     &   (mt - Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(mt*(mt - Ph(1) - Pt(1))*(kapt*UBt(1) -&
     &   dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) + Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) -&
     &   kapt*mt*UBt(3)) + (Ph(2) + Pt(2) + dcmplx(0d0,1d0)*(Ph(3) +&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) - kapt*mt*UBt(4)))) + (Pb(4) -&
     &   Pw(4))*(mt*Eg(1)*((Ph(4) + Pt(4))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (mt +&
     &   Ph(1) + Pt(1))*(dcmplx(0d0,-1d0)*kap*UBt(1) + kapt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4))) +&
     &   mt*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pt(2) + Pt(3))*(kap*UBt(1) + dcmplx(0d0,1d0)*kapt*UBt(3)) + (mt -&
     &   Ph(1) - Pt(1))*(kapt*UBt(2) - dcmplx(0d0,1d0)*kap*UBt(4)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,-1d0)*kap*UBt(2) + kapt*UBt(4))) + Eg(4)*(mt*(mt - Ph(1) -&
     &   Pt(1))*(kapt*UBt(1) - dcmplx(0d0,1d0)*kap*UBt(3)) + (Ph(4) +&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*mt*UBt(1) - kapt*mt*UBt(3)) + (Ph(2) + Pt(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pt(3)))*(dcmplx(0d0,1d0)*kap*mt*UBt(2) -&
     &   kapt*mt*UBt(4))))))))/(sqrt2*v**2*(mt**2 - (Ph(1) + Pt(1))**2 + (Ph(2) + Pt(2))**2&
     &   + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)*(mt**2 - (Pb(1) - Pw(1))**2 + (Pb(2) -&
     &   Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))
     
        ampwt = (gs*mw*(dcmplx(0d0,-1d0)*DCONJG(UBb(1))*((a1WW*mw**2*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) +&
     &   (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3))&
     &   - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) -&
     &   ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*(Eg(2)*Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) +&
     &   Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) - Eg(2)*Pt(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) -&
     &   Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) - Eg(4)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) - (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) +&
     &   Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1))&
     &   + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) +&
     &   Pw(4)))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) -&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) +&
     &   Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4)))) - dcmplx(0d0,1d0)*DCONJG(UBb(3))*((a1WW*mw**2*(ECw(2)&
     &   + dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) +&
     &   (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3))&
     &   - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) -&
     &   ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))*(Eg(2)*Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) +&
     &   Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) - Eg(2)*Pt(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) -&
     &   Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) - Eg(4)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) - (a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) +&
     &   Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1))&
     &   + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) +&
     &   Pw(4)))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) -&
     &   (a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) +&
     &   Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,-1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) -&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 - dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 - dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) - dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4)))) - DCONJG(UBb(4))*(dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(1)&
     &   - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) +&
     &   (a1WW*(mw**2*(dcmplx(0d0,1d0)*ECw(2) + ECw(3)) + (dcmplx(0d0,1d0)*Ph(2) + Ph(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) + Pw(3))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) +&
     &   ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4))) +&
     &   dcmplx(0d0,2d0)*(a4WW*(ECw(3)*Ph(4)*Pw(1) - dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) -&
     &   ECw(1)*Ph(4)*Pw(3) + ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) +&
     &   Ph(1)*(dcmplx(0d0,1d0)*Pw(2) + Pw(3))) - ECw(3)*Ph(1)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4))))))*(Eg(2)*Pg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) + Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) -&
     &   Eg(2)*Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) - Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) -&
     &   Eg(4)*Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) - dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pw(2)&
     &   + Pw(3)))*(ECw(3)*Ph(3) + ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) + ECw(2)*(Ph(2) +&
     &   Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) -&
     &   dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4)))) + DCONJG(UBb(2))*(dcmplx(0d0,-1d0)*(a1WW*mw**2*(ECw(1)&
     &   - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) - Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) -&
     &   2*a2WW*(ECw(3)*Ph(3)*Pw(1) - ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3)&
     &   + ECw(3)*Pw(1)*Pw(3) - ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) +&
     &   ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) -&
     &   ECw(3)*Pw(3)*Pw(4) - ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2&
     &   + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1)&
     &   - ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1)&
     &   - Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(1) +&
     &   mt*Eg(1)*UBt(2) - Eg(1)*Pg(1)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(2) +&
     &   Eg(3)*Pg(3)*UBt(2) + Eg(1)*Pt(1)*UBt(2) - dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(2) -&
     &   Eg(3)*Pt(3)*UBt(2) - dcmplx(0d0,1d0)*mt*Eg(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(3) + Eg(1)*Pg(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(3) -&
     &   Eg(1)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(3) + Eg(2)*(-(Pg(4)*UBt(1)) +&
     &   Pt(4)*UBt(1) + Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Pg(3)*UBt(2) - Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(2) + mt*UBt(3) - Pg(1)*UBt(3) + Pt(1)*UBt(3)) -&
     &   Eg(1)*Pg(4)*UBt(4) + Eg(1)*Pt(4)*UBt(4) + Eg(4)*(Pg(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(1) - Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Pt(3)*UBt(1) +&
     &   Pg(4)*UBt(2) - Pt(4)*UBt(2) - mt*UBt(4) + Pg(1)*UBt(4) - Pt(1)*UBt(4))) -&
     &   dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + a1WW*(Ph(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pw(2) + Pw(3)))*(ECw(3)*Ph(3) + ECw(4)*Ph(4)&
     &   - ECw(1)*(Ph(1) + Pw(1)) + ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) + ECw(4)*Pw(4)) +&
     &   2*a4WW*(ECw(3)*Ph(4)*Pw(1) - dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   2*a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4)))))*(Eg(2)*Pg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(1) + Eg(3)*Pg(3)*UBt(1) + Eg(4)*Pg(4)*UBt(1) -&
     &   Eg(2)*Pt(2)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(1) -&
     &   dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(1) - Eg(3)*Pt(3)*UBt(1) - Eg(4)*Pt(4)*UBt(1) -&
     &   Eg(4)*Pg(2)*UBt(2) - dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(2) + Eg(2)*Pg(4)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(2) + Eg(4)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(2) - Eg(2)*Pt(4)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(2) - Eg(4)*Pg(1)*UBt(3) + Eg(4)*Pt(1)*UBt(3) -&
     &   Eg(2)*Pg(1)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(4) + Eg(2)*Pt(1)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(4) + mt*(Eg(1)*UBt(1) + Eg(4)*UBt(3) + (Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*UBt(4)) + Eg(1)*(-(Pg(1)*UBt(1)) + Pt(1)*UBt(1) + Pg(4)*UBt(3)&
     &   - Pt(4)*UBt(3) + Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))) + (a1WW*(mw**2*(dcmplx(0d0,1d0)*ECw(2) + ECw(3)) +&
     &   (dcmplx(0d0,1d0)*Ph(2) + Ph(3) + dcmplx(0d0,1d0)*Pw(2) + Pw(3))*(ECw(3)*Ph(3) +&
     &   ECw(4)*Ph(4) - ECw(1)*(Ph(1) + Pw(1)) + ECw(2)*(Ph(2) + Pw(2)) + ECw(3)*Pw(3) +&
     &   ECw(4)*Pw(4))) + dcmplx(0d0,2d0)*(a4WW*(ECw(3)*Ph(4)*Pw(1) - dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2)&
     &   - ECw(1)*Ph(4)*Pw(3) + ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) +&
     &   Ph(1)*(dcmplx(0d0,1d0)*Pw(2) + Pw(3))) - ECw(3)*Ph(1)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4) +&
     &   dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4))) +&
     &   a2WW*(dcmplx(0d0,1d0)*ECw(3)*(Ph(1)*Pw(1) + Pw(1)**2 - Ph(2)*Pw(2) +&
     &   dcmplx(0d0,1d0)*Ph(3)*Pw(2) - Pw(2)**2 + dcmplx(0d0,1d0)*Pw(2)*Pw(3) - Ph(4)*Pw(4) -&
     &   Pw(4)**2) + ECw(2)*(-(Ph(1)*Pw(1)) - Pw(1)**2 + dcmplx(0d0,1d0)*Ph(2)*Pw(3) +&
     &   Ph(3)*Pw(3) + dcmplx(0d0,1d0)*Pw(2)*Pw(3) + Pw(3)**2 + Ph(4)*Pw(4) + Pw(4)**2) +&
     &   (Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(4)*(Ph(4) +&
     &   Pw(4))))))*(-(Eg(1)*Pg(4)*UBt(1)) + Eg(1)*Pt(4)*UBt(1) + Eg(2)*Pg(1)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(1)*UBt(2) - Eg(1)*Pg(2)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(2) - Eg(2)*Pt(1)*UBt(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(2) + Eg(1)*Pt(2)*UBt(2) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(2) + Eg(1)*Pg(1)*UBt(3) - Eg(2)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(2)*Pg(3)*UBt(3) -&
     &   Eg(3)*Pg(3)*UBt(3) - Eg(1)*Pt(1)*UBt(3) + Eg(2)*Pt(2)*UBt(3) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(3) + dcmplx(0d0,1d0)*Eg(2)*Pt(3)*UBt(3) +&
     &   Eg(3)*Pt(3)*UBt(3) + mt*(Eg(4)*UBt(1) + Eg(2)*UBt(2) + dcmplx(0d0,1d0)*Eg(3)*UBt(2) +&
     &   Eg(1)*UBt(3)) - Eg(2)*Pg(4)*UBt(4) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(4) +&
     &   Eg(2)*Pt(4)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(4) + Eg(4)*(Pg(1)*UBt(1) -&
     &   Pt(1)*UBt(1) - Pg(4)*UBt(3) + Pt(4)*UBt(3) + Pg(2)*UBt(4) +&
     &   dcmplx(0d0,1d0)*Pg(3)*UBt(4) - Pt(2)*UBt(4) - dcmplx(0d0,1d0)*Pt(3)*UBt(4))) +&
     &   dcmplx(0d0,1d0)*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4)))))*(dcmplx(0d0,-1d0)*Eg(3)*Pg(1)*UBt(1) - Eg(1)*Pg(2)*UBt(1) +&
     &   dcmplx(0d0,1d0)*Eg(1)*Pg(3)*UBt(1) + dcmplx(0d0,1d0)*Eg(3)*Pt(1)*UBt(1) +&
     &   Eg(1)*Pt(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(1)*Pt(3)*UBt(1) - Eg(4)*Pg(1)*UBt(2) +&
     &   Eg(1)*Pg(4)*UBt(2) + Eg(4)*Pt(1)*UBt(2) - Eg(1)*Pt(4)*UBt(2) - Eg(4)*Pg(2)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(4)*Pg(3)*UBt(3) - dcmplx(0d0,1d0)*Eg(3)*Pg(4)*UBt(3) +&
     &   Eg(4)*Pt(2)*UBt(3) - dcmplx(0d0,1d0)*Eg(4)*Pt(3)*UBt(3) +&
     &   dcmplx(0d0,1d0)*Eg(3)*Pt(4)*UBt(3) + Eg(1)*Pg(1)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Eg(3)*Pg(2)*UBt(4) - Eg(3)*Pg(3)*UBt(4) - Eg(4)*Pg(4)*UBt(4) -&
     &   Eg(1)*Pt(1)*UBt(4) + dcmplx(0d0,1d0)*Eg(3)*Pt(2)*UBt(4) + Eg(3)*Pt(3)*UBt(4) +&
     &   Eg(4)*Pt(4)*UBt(4) + mt*(Eg(2)*UBt(1) - dcmplx(0d0,1d0)*Eg(3)*UBt(1) - Eg(4)*UBt(2) +&
     &   Eg(1)*UBt(4)) + Eg(2)*(Pg(1)*UBt(1) - Pt(1)*UBt(1) + Pg(4)*UBt(3) - Pt(4)*UBt(3) -&
     &   Pg(2)*UBt(4) + dcmplx(0d0,1d0)*Pg(3)*UBt(4) + Pt(2)*UBt(4) -&
     &   dcmplx(0d0,1d0)*Pt(3)*UBt(4))))))/(sqrt2*v**2*(mt**2 - Pg(1)**2 + Pg(2)**2 +&
     &   Pg(3)**2 + Pg(4)**2 + 2*Pg(1)*Pt(1) - Pt(1)**2 - 2*Pg(2)*Pt(2) + Pt(2)**2 -&
     &   2*Pg(3)*Pt(3) + Pt(3)**2 - 2*Pg(4)*Pt(4) + Pt(4)**2)*(mw**2 - Ph(1)**2 + Ph(2)**2 +&
     &   Ph(3)**2 + Ph(4)**2 - 2*Ph(1)*Pw(1) - Pw(1)**2 + 2*Ph(2)*Pw(2) + Pw(2)**2 +&
     &   2*Ph(3)*Pw(3) + Pw(3)**2 + 2*Ph(4)*Pw(4) + Pw(4)**2))

         amp = ampts + ampws + amptt1 + amptt2 + ampwt

      end subroutine gb_twmHamp      


      subroutine gbb_tbwpHamp(Pg,Eg,Pb,Vb,Pt,Vt,Pw,ECw,Ph,TTBHcoupl,amp)
! amplitude for production g(Pg)+b~(Pb)->t~(Pt)+wp(Pw)+h(Ph)
! allowing for scalar & pseudoscalar couplings of Higgs to top (Spinors and Pol-Vectors in Dirac representation!)
        implicit none
        real(8)    :: Pg(4),Pb(4),Pt(4),Pw(4),Ph(4)
        complex(8) :: Vb(4),Eg(4),Vt(4),ECw(4)
        complex(8) :: amp,ampts,ampws,amptt1,amptt2,ampwt
        real(8)    :: mt,mw,v, gs
        complex(8)    :: kap,kapt,a1WW,a2WW,a4WW  
        complex(8) TTBHcoupl(1:2)
        include 'includeVars.F90'            


        mw=M_W
        mt=m_Top
        v=vev
        
        gs=dsqrt(4d0*Pi*alphas);
        
        kap=TTBHcoupl(1)
        kapt=TTBHcoupl(2)
        a1WW=(2d0,0d0)
        a2WW=(0d0,0d0)
        a4WW=(0d0,0d0)
        

        ampts = -((gs*mt*mw*(dcmplx(0d0,1d0)*kapt*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(1) - kap*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(1) - kap*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) - ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(2) + dcmplx(0d0,1d0)*kapt*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) + ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(2) + kap*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(3) - dcmplx(0d0,1d0)*kapt*(((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(3) - dcmplx(0d0,1d0)*kapt*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt - Ph(1) - Pt(1)) + ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) - ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(4) + kap*(((ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(mt + Ph(1) + Pt(1)) - ((ECw(1) +&
     &   ECw(4))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2)&
     &   + Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(2) + dcmplx(0d0,1d0)*ECw(3))*((DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) +&
     &   ECw(4))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3))) + ((ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3))*((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4))) +&
     &   dcmplx(0d0,1d0)*(ECw(1) - ECw(4))*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4))) - (ECw(1) -&
     &   ECw(4))*(dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4))))*(Ph(4) + Pt(4)))*Vt(4)))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) +&
     &   Pg(2))**2 + (Pb(3) + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mt**2 - (Ph(1) + Pt(1))**2 +&
     &   (Ph(2) + Pt(2))**2 + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)))

        ampws = -((gs*mw*(-((dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) -&
     &   (-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(1))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3)))&
     &   - (DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   ((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) - (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) + ECw(4))&
     &   - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1)&
     &   + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))*Vt(1)) -&
     &   (((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) +&
     &   Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) - (-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - (DCONJG(Vb(4))*Eg(1) -&
     &   DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))*Vt(2) +&
     &   (dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) - (-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - (DCONJG(Vb(4))*Eg(1) -&
     &   DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   ((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) +&
     &   ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) +&
     &   ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) - (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(1) + ECw(4))&
     &   - a1WW*(Ph(1) + Ph(4) + Pw(1) + Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1)&
     &   + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) +&
     &   2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) + ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 +&
     &   ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) + ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) +&
     &   ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) + Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) -&
     &   Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 -&
     &   Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) + ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) +&
     &   ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) + Pw(4))) + ECw(2)*(-((Ph(1) +&
     &   Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))*Vt(3) +&
     &   (((dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) +&
     &   Eg(3)) - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(2))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) -&
     &   dcmplx(0d0,1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   (dcmplx(0d0,-1d0)*(-(DCONJG(Vb(3))*Eg(1)) + DCONJG(Vb(2))*(Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(1) + Pg(1)) +&
     &   dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) +&
     &   DCONJG(Vb(4))*Eg(4))*(Pb(2) + dcmplx(0d0,1d0)*Pb(3) + Pg(2) + dcmplx(0d0,1d0)*Pg(3)) +&
     &   (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) + DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3))&
     &   - dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(4) + Pg(4)))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   dcmplx(0d0,1d0)*((DCONJG(Vb(2))*Eg(1) - DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(4))*Eg(4))*(Pb(1) + Pg(1)) - (-(DCONJG(Vb(3))*Eg(1)) +&
     &   DCONJG(Vb(2))*(Eg(2) + dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(1))*Eg(4))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - (DCONJG(Vb(4))*Eg(1) -&
     &   DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(2))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   (dcmplx(0d0,1d0)*(DCONJG(Vb(4))*Eg(1) - DCONJG(Vb(1))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3))&
     &   + DCONJG(Vb(2))*Eg(4))*(Pb(1) + Pg(1)) + (dcmplx(0d0,1d0)*DCONJG(Vb(1))*Eg(1) +&
     &   DCONJG(Vb(4))*(dcmplx(0d0,-1d0)*Eg(2) + Eg(3)) -&
     &   dcmplx(0d0,1d0)*DCONJG(Vb(3))*Eg(4))*(Pb(2) - dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pg(2) + Pg(3))) - dcmplx(0d0,1d0)*(DCONJG(Vb(2))*Eg(1) -&
     &   DCONJG(Vb(3))*(Eg(2) - dcmplx(0d0,1d0)*Eg(3)) + DCONJG(Vb(4))*Eg(4))*(Pb(4) +&
     &   Pg(4)))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) +&
     &   Pw(4))))))*Vt(4)))/(sqrt2*v**2*(-(Pb(1) + Pg(1))**2 + (Pb(2) + Pg(2))**2 + (Pb(3)&
     &   + Pg(3))**2 + (Pb(4) + Pg(4))**2)*(mw**2 - (Ph(1) + Pw(1))**2 + (Ph(2) + Pw(2))**2 +&
     &   (Ph(3) + Pw(3))**2 + (Ph(4) + Pw(4))**2)))

        amptt1 = (gs*mt*mw*(-((dcmplx(0d0,1d0)*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((mt - Pg(1)&
     &   + Pt(1))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kapt*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))))) - Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   dcmplx(0d0,1d0)*(mt - Pg(1) + Pt(1))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + dcmplx(0d0,1d0)*Eg(1)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Pg(4) - Pt(4))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(1)) + (Eg(4)*(dcmplx(0d0,1d0)*(mt - Pg(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - (Pg(2) - dcmplx(0d0,1d0)*Pg(3) -&
     &   Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   dcmplx(0d0,1d0)*(mt - Pg(1) + Pt(1))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + dcmplx(0d0,1d0)*Eg(1)*((Pg(4) - Pt(4))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) - kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))))&
     &   + dcmplx(0d0,1d0)*(mt + Pg(1) - Pt(1))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))))&
     &   + (Pg(2) - dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(2) - (Eg(1)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(4) -&
     &   Pt(4))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   dcmplx(0d0,1d0)*(mt - Pg(1) + Pt(1))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + dcmplx(0d0,1d0)*(Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(2) - dcmplx(0d0,1d0)*Pg(3) -&
     &   Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   - dcmplx(0d0,1d0)*Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Pg(4) - Pt(4))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(3) + dcmplx(0d0,1d0)*(Eg(1)*((mt - Pg(1) +&
     &   Pt(1))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kapt*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))))) + Eg(4)*((Pg(4) -&
     &   Pt(4))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) + (Pg(2) - dcmplx(0d0,1d0)*Pg(3) -&
     &   Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))&
     &   + (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   kap*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - dcmplx(0d0,1d0)*(mt + Pg(1) -&
     &   Pt(1))*(kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   dcmplx(0d0,1d0)*kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Pg(4) - Pt(4))*(kap*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   dcmplx(0d0,1d0)*kapt*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4))))))*Vt(4)))/(sqrt2*v**2*(mt**2 - (Pg(1) - Pt(1))**2 + (Pg(2) - Pt(2))**2 +&
     &   (Pg(3) - Pt(3))**2 + (Pg(4) - Pt(4))**2)*(mt**2 - (Pb(1) - Pw(1))**2 + (Pb(2) -&
     &   Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))

        amptt2 = (gs*mt*mw*(kapt*((Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   Eg(1)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(1) + dcmplx(0d0,1d0)*kap*((mt - Ph(1) - Pt(1))*((Eg(2)&
     &   + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) + (Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))*Vt(1) - kapt*((Ph(4) +&
     &   Pt(4))*(Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) + (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) +&
     &   Pt(3)))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*(Eg(4)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4)))))*Vt(2) - dcmplx(0d0,1d0)*kap*((mt - Ph(1) - Pt(1))*(Eg(1)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(4) + Pt(4))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3)))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(2) - dcmplx(0d0,1d0)*kap*((Ph(2) +&
     &   dcmplx(0d0,1d0)*Ph(3) + Pt(2) + dcmplx(0d0,1d0)*Pt(3))*(Eg(1)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   Eg(1)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(3) - kapt*((mt - Ph(1) - Pt(1))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) + (Ph(2) + dcmplx(0d0,1d0)*Ph(3) + Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) -&
     &   (Ph(4) + Pt(4))*((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) +&
     &   Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))))*Vt(3) +&
     &   dcmplx(0d0,1d0)*kap*((Ph(4) + Pt(4))*(Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2)&
     &   - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3)))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(4)*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) - ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - Eg(1)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))) - (mt + Ph(1) + Pt(1))*(Eg(4)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) -&
     &   Pw(4)))))*Vt(4) + kapt*((mt - Ph(1) - Pt(1))*(Eg(1)*(((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) +&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) -&
     &   (Eg(2) - dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt + Pb(1) - Pw(1)) - ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) -&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(4) + Pt(4))*(Eg(4)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1)) - ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(2) -&
     &   dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) - (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) +&
     &   ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4)))) +&
     &   (Ph(2) - dcmplx(0d0,1d0)*(Ph(3) + dcmplx(0d0,1d0)*Pt(2) + Pt(3)))*((Eg(2) +&
     &   dcmplx(0d0,1d0)*Eg(3))*(((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) +&
     &   ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(2) - dcmplx(0d0,1d0)*Pb(3) - Pw(2) + dcmplx(0d0,1d0)*Pw(3)) -&
     &   ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) +&
     &   Eg(1)*(-(((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt + Pb(1) -&
     &   Pw(1))) + ((DCONJG(Vb(1)) - DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   DCONJG(Vb(2))*(ECw(1) - ECw(4)) + DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) +&
     &   dcmplx(0d0,1d0)*(Pb(3) + dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(Pb(4) - Pw(4))) + Eg(4)*(((DCONJG(Vb(2)) -&
     &   DCONJG(Vb(4)))*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) -&
     &   DCONJG(Vb(3))*(ECw(1) + ECw(4)))*(mt - Pb(1) + Pw(1)) + ((DCONJG(Vb(1)) -&
     &   DCONJG(Vb(3)))*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(2))*(ECw(1) - ECw(4)) +&
     &   DCONJG(Vb(4))*(-ECw(1) + ECw(4)))*(Pb(2) + dcmplx(0d0,1d0)*(Pb(3) +&
     &   dcmplx(0d0,1d0)*Pw(2) - Pw(3))) + ((DCONJG(Vb(2)) - DCONJG(Vb(4)))*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3)) + DCONJG(Vb(1))*(ECw(1) + ECw(4)) - DCONJG(Vb(3))*(ECw(1) +&
     &   ECw(4)))*(Pb(4) - Pw(4)))))*Vt(4)))/(sqrt2*v**2*(mt**2 - (Ph(1) + Pt(1))**2 +&
     &   (Ph(2) + Pt(2))**2 + (Ph(3) + Pt(3))**2 + (Ph(4) + Pt(4))**2)*(mt**2 - (Pb(1) -&
     &   Pw(1))**2 + (Pb(2) - Pw(2))**2 + (Pb(3) - Pw(3))**2 + (Pb(4) - Pw(4))**2))
     
        ampwt = (dcmplx(0d0,-1d0)*gs*mw*((-((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*((mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))))) +&
     &   Eg(1)*(-((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) + (mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) -&
     &   Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))))*Vt(1) -&
     &   (Eg(1)*(-((mt + Pg(1) - Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) - Eg(4)*((mt -&
     &   Pg(1) + Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))))*Vt(2) +&
     &   ((Eg(2) + dcmplx(0d0,1d0)*Eg(3))*(-((mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) +&
     &   Eg(4)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) +&
     &   Eg(1)*((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (mt - Pg(1) +&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) + (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))))*Vt(3) -&
     &   (Eg(4)*(-((mt + Pg(1) - Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) +&
     &   Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   a1WW*(Ph(2) + Pw(2) - dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) -&
     &   ECw(2)*(Ph(2) + Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) -&
     &   2*a2WW*(ECw(2) - dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) +&
     &   Pw(2)) - Pw(3)*(Ph(3) + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) - Eg(1)*((mt -&
     &   Pg(1) + Pt(1))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) - dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4)))))) + (Pg(2) -&
     &   dcmplx(0d0,1d0)*Pg(3) - Pt(2) +&
     &   dcmplx(0d0,1d0)*Pt(3))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) +&
     &   2*a2WW*(Pw(2) + dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))))) + (Eg(2) -&
     &   dcmplx(0d0,1d0)*Eg(3))*(-((Pg(2) - Pt(2) + dcmplx(0d0,1d0)*(Pg(3) -&
     &   Pt(3)))*(DCONJG(Vb(1))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(2) - dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) -&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) -&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) -&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) -&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,-1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) + dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   + dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(2))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(1) - ECw(4)) - a1WW*(Ph(1) - Ph(4) + Pw(1) -&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) - 2*a2WW*(ECw(3)*Ph(3)*Pw(1) -&
     &   ECw(1)*Ph(2)*Pw(2) - ECw(1)*Pw(2)**2 - ECw(1)*Ph(3)*Pw(3) + ECw(3)*Pw(1)*Pw(3) -&
     &   ECw(1)*Pw(3)**2 + ECw(2)*(Ph(2) + Pw(2))*(Pw(1) - Pw(4)) + ECw(1)*Ph(1)*Pw(4) -&
     &   ECw(3)*Ph(3)*Pw(4) - ECw(1)*Ph(4)*Pw(4) + ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) -&
     &   ECw(1)*Pw(4)**2 + ECw(4)*(-(Ph(1)*Pw(1)) + Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2) +&
     &   Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 + Pw(1)*Pw(4))) - 2*a4WW*(-((ECw(1) -&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(2)*((-Ph(1) + Ph(4))*Pw(3) + Ph(3)*(Pw(1) -&
     &   Pw(4))) + ECw(3)*((Ph(1) - Ph(4))*Pw(2) + Ph(2)*(-Pw(1) + Pw(4))))))) + (mt + Pg(1) -&
     &   Pt(1))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4)))))) - (Pg(4) -&
     &   Pt(4))*(DCONJG(Vb(2))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) -&
     &   DCONJG(Vb(4))*(a1WW*mw**2*(ECw(2) + dcmplx(0d0,1d0)*ECw(3)) + 2*a2WW*(Pw(2) +&
     &   dcmplx(0d0,1d0)*Pw(3))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - a1WW*(Ph(2) + Pw(2) +&
     &   dcmplx(0d0,1d0)*(Ph(3) + Pw(3)))*(ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) + Pw(2)) -&
     &   ECw(3)*(Ph(3) + Pw(3)) - ECw(4)*(Ph(4) + Pw(4))) - 2*a2WW*(ECw(2) +&
     &   dcmplx(0d0,1d0)*ECw(3))*(Pw(1)*(Ph(1) + Pw(1)) - Pw(2)*(Ph(2) + Pw(2)) - Pw(3)*(Ph(3)&
     &   + Pw(3)) - Pw(4)*(Ph(4) + Pw(4))) + 2*a4WW*(ECw(3)*Ph(4)*Pw(1) +&
     &   dcmplx(0d0,1d0)*ECw(1)*Ph(4)*Pw(2) - ECw(1)*Ph(4)*Pw(3) +&
     &   ECw(4)*(dcmplx(0d0,1d0)*Ph(2)*Pw(1) - Ph(3)*Pw(1) + Ph(1)*(dcmplx(0d0,-1d0)*Pw(2) +&
     &   Pw(3))) - ECw(3)*Ph(1)*Pw(4) - dcmplx(0d0,1d0)*ECw(1)*Ph(2)*Pw(4) + ECw(1)*Ph(3)*Pw(4)&
     &   - dcmplx(0d0,1d0)*ECw(2)*(Ph(4)*Pw(1) - Ph(1)*Pw(4)))) +&
     &   DCONJG(Vb(1))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) + Pw(4))))) -&
     &   DCONJG(Vb(3))*(a1WW*mw**2*(ECw(1) + ECw(4)) - a1WW*(Ph(1) + Ph(4) + Pw(1) +&
     &   Pw(4))*(-(ECw(3)*Ph(3)) - ECw(4)*Ph(4) + ECw(1)*(Ph(1) + Pw(1)) - ECw(2)*(Ph(2) +&
     &   Pw(2)) - ECw(3)*Pw(3) - ECw(4)*Pw(4)) + 2*a2WW*(-(ECw(3)*Ph(3)*Pw(1)) +&
     &   ECw(1)*Ph(2)*Pw(2) + ECw(1)*Pw(2)**2 + ECw(1)*Ph(3)*Pw(3) - ECw(3)*Pw(1)*Pw(3) +&
     &   ECw(1)*Pw(3)**2 + ECw(1)*Ph(1)*Pw(4) - ECw(3)*Ph(3)*Pw(4) + ECw(1)*Ph(4)*Pw(4) +&
     &   ECw(1)*Pw(1)*Pw(4) - ECw(3)*Pw(3)*Pw(4) + ECw(1)*Pw(4)**2 - ECw(2)*(Ph(2) +&
     &   Pw(2))*(Pw(1) + Pw(4)) + ECw(4)*(-(Ph(1)*Pw(1)) - Ph(4)*Pw(1) - Pw(1)**2 + Ph(2)*Pw(2)&
     &   + Pw(2)**2 + Ph(3)*Pw(3) + Pw(3)**2 - Pw(1)*Pw(4))) + 2*a4WW*(-((ECw(1) +&
     &   ECw(4))*(Ph(3)*Pw(2) - Ph(2)*Pw(3))) + ECw(3)*((Ph(1) + Ph(4))*Pw(2) - Ph(2)*(Pw(1) +&
     &   Pw(4))) + ECw(2)*(-((Ph(1) + Ph(4))*Pw(3)) + Ph(3)*(Pw(1) +&
     &   Pw(4))))))))*Vt(4)))/(sqrt2*v**2*(mt**2 - (Pg(1) - Pt(1))**2 + (Pg(2) - Pt(2))**2&
     &   + (Pg(3) - Pt(3))**2 + (Pg(4) - Pt(4))**2)*(mw**2 - (Ph(1) + Pw(1))**2 + (Ph(2) +&
     &   Pw(2))**2 + (Ph(3) + Pw(3))**2 + (Ph(4) + Pw(4))**2))

         amp = ampts + ampws + amptt1 + amptt2 + ampwt

      end subroutine gbb_tbwpHamp



    
          subroutine ubarSpi_Dirac(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8)  p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2=p0+m
          fc=cdsqrt( dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=fc
            f(2)=dcmplx(0d0,0d0)
            f(3)=-1d0*pz*fc/fc2
            f(4)=-(px-(0d0,1d0)*py)*fc/fc2
          elseif (i.eq.-1) then
            f(1)=dcmplx(0d0,0d0)
            f(2)=fc
            f(3)=-(px+(0d0,1d0)*py)*fc/fc2
            f(4)=pz*fc/fc2
          else
              print *, "wrong helicity setting in ubarSpi"
              stop
          endif

          return
          end subroutine



          subroutine vSpi_Dirac(p,m,i,f)
          implicit none
          integer i
          real(8) m
          complex(8) p(4)
          complex(8) f(4),fc
          real(8) p0,px,py,pz,fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

          fc2 = p0+m
          fc=cdsqrt(dcmplx(fc2))
!           fc=dsqrt(fc2)

          if (i.eq.1) then
            f(1)=pz*fc/fc2
            f(2)=(px+(0d0,1d0)*py)*fc/fc2
            f(3)=fc
            f(4)=dcmplx(0d0,0d0)
          elseif (i.eq.-1) then
            f(1)=(px-(0d0,1d0)*py)*fc/fc2
            f(2)=-pz*fc/fc2
            f(3)=dcmplx(0d0,0d0)
            f(4)=fc
          else
              print *, "wrong helicity setting in vSpi"
          endif

          return
          end SUBROUTINE



SUBROUTINE TopDecay(Flavor,Mom,Spinor,TopHel)
implicit none
real(8) :: Mom(1:4,1:3)
integer :: flavor
integer,optional :: TopHel
complex(8) :: Spinor(1:4)
real(8) :: TopMom(1:4),NWAFactor_Top
complex(8) :: Spi(1:4),BarSpi(1:4),BotSpi(1:4),WCurr(1:4)
real(8) :: NWAFactor_W
complex(8) :: WProp
include 'includeVars.F90'

NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
WProp = (0d0,-1d0)*NWAFactor_W



    TopMom(1:4) = Mom(1:4,1)+Mom(1:4,2)+Mom(1:4,3)

    NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top)


    if( Flavor.eq.Top_ ) then ! Top quark decay
!       assemble lepton current
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! nu or up
        WCurr(1:4)  = vbqq_Weyl(BarSpi(1:4),Spi(1:4)) * WProp * gwsq ! vbqq introduces -i/Sqrt(2)

!       connect to quark current
        BarSpi(1:4) = BotSpi(1:4)
        Spinor(1:4) = vgq_Weyl( WCurr(1:4),BarSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        Spinor(1:4) =( spb2_Weyl(Spinor(1:4),dcmplx(TopMom(1:4))) + m_Top*Spinor(1:4) ) * NWAFactor_Top
        Spinor(1:4) = WeylToDirac(Spinor(1:4))
    elseif( Flavor.eq.ATop_ ) then ! Anti-Top quark decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! nubar or up_bar
        WCurr(1:4)  = vbqq_Weyl(BarSpi(1:4),Spi(1:4)) * WProp * gwsq ! vbqq introduces -i/Sqrt(2)

!       connect to quark current:
        Spi(1:4) = BotSpi(1:4)
        Spinor(1:4) = vbqg_Weyl( Spi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
        Spinor(1:4) = ( spi2_Weyl(dcmplx(TopMom(1:4)),Spinor(1:4)) - m_Top*Spinor(1:4) ) * NWAFactor_Top
        Spinor(1:4) = WeylToDirac(Spinor(1:4))
    endif


RETURN
END SUBROUTINE






SUBROUTINE WDecay(Charge,Mom,WCurr)
implicit none
real(8) :: Mom(1:4,1:2)
integer :: Charge
real(8) :: WMom(1:4)
complex(8) :: Spi(1:4),BarSpi(1:4),BotSpi(1:4),WCurr(1:4)
real(8) :: NWAFactor_W
complex(8) :: WProp
include 'includeVars.F90'

NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
WProp = (0d0,-1d0)*dsqrt(gwsq)*NWAFactor_W



    WMom(1:4) = Mom(1:4,1)+Mom(1:4,2)


    if( Charge.eq.Wp_ ) then ! W+ decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! nu or up
        WCurr(1:4)  = vbqq_Weyl(BarSpi(1:4),Spi(1:4)) * WProp ! vbqq introduces -i/Sqrt(2)
        
    elseif( Charge.eq.Wm_ ) then ! W- decay
!       assemble lepton current
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! nubar or up_bar
        WCurr(1:4)  = vbqq_Weyl(BarSpi(1:4),Spi(1:4)) * WProp ! vbqq introduces -i/Sqrt(2)
        
    endif


RETURN
END SUBROUTINE




          SUBROUTINE ubarSpi_Weyl(p,i,ubarSpi)  ! i=+1 is ES to Chir_Weyl(.false.), i=-1 is ES to Chir_Weyl(.true.)
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: ubarSpi(4)
          complex(8) :: ephi
          real(8) :: p0,px,py,pz
          complex(8) :: fc, fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

         fc2 = p0 + pz
         fc=cdsqrt(fc2)

         if (cdabs(fc2).gt.1D-15) then

         if (i.eq.1) then
            ubarSpi(1)=(0d0,0d0)
            ubarSpi(2)=(0d0,0d0)
            ubarSpi(3)=fc
            ubarSpi(4)=(px-(0d0,1d0)*py)/fc
         elseif (i.eq.-1) then
            ubarSpi(1)=(px+(0d0,1d0)*py)/fc
            ubarSpi(2)=-fc
            ubarSpi(3)=(0d0,0d0)
            ubarSpi(4)=(0d0,0d0)
         else
          call Error("wrong helicity setting in ubarSpi_Weyl")
         endif

         else

         if (i.eq.1) then
            ubarSpi(1) = (0d0,0d0)
            ubarSpi(2) = (0d0,0d0)
            ubarSpi(3) = (0d0,0d0)
            ubarSpi(4) = dsqrt(2d0*p0)
         elseif (i.eq.-1) then
            ubarSpi(1) = dsqrt(2d0*p0)
            ubarSpi(2) = (0d0,0d0)
            ubarSpi(3) = (0d0,0d0)
            ubarSpi(4) = (0d0,0d0)
         else
            call Error("wrong helicity setting in ubarSpi_Weyl")
         endif

         endif
        return
        END SUBROUTINE





          SUBROUTINE vSpi_Weyl(p,i,vSpi)  ! i=+1 is ES to Chir_Weyl(.false.), i=-1 is ES to Chir_Weyl(.true.)
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4)
          complex(8) :: ephi
          real(8) :: p0,px,py,pz
          real(8) :: nx,ny,nz,theta,phi
          real(8) :: ct,ct2,st,st2,cphi,sphi
          complex(8) :: fc2, fc

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

         fc2 = p0 + pz
         fc=cdsqrt(fc2)

         if (cdabs(fc2).gt.1D-15) then

         if (i.eq.1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=(0d0,0d0)
            vSpi(3)=(px-(0d0,1d0)*py)/fc
            vSpi(4)=-fc
         elseif (i.eq.-1) then
            vSpi(1)=fc
            vSpi(2)=(px+(0d0,1d0)*py)/fc
            vSpi(3)=(0d0,0d0)
            vSpi(4)=(0d0,0d0)
         else
            call Error("wrong helicity setting in vSpi_Weyl")
         endif

         else

         if (i.eq.1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=(0d0,0d0)
            vSpi(3)=dsqrt(2d0*p0)
            vSpi(4)=(0d0,0d0)
         elseif (i.eq.-1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=dsqrt(2d0*p0)
            vSpi(3)=(0d0,0d0)
            vSpi(4)=(0d0,0d0)
         else
            call Error("wrong helicity setting in vSpi_Weyl")
         endif

         endif

         RETURN
         END SUBROUTINE



  function pol_mless(p,i,outgoing)
  implicit none
    complex(8), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------
    integer :: pol
    real(8) :: p0,px,py,pz
    real(8) :: pv,ct,st,cphi,sphi
    complex(8) :: pol_mless(4)
     real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
    complex(8), parameter :: czero = dcmplx(0d0,0d0), ci=dcmplx(0d0,1d0)

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),8)
    px=real(p(2),8)
    py=real(p(3),8)
    pz=real(p(4),8)
!^^^END


    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0d0-ct**2))

    if (st < 0.0000001d0) then
       cphi=1.0d0
       sphi=0.0d0
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0d0) then
       pol=i
    else
       pol=-i
    endif

    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif

    pol_mless(1)=czero
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless
  
  
   function pol_mass(p,i,outgoing)
   implicit none
   integer, intent(in) :: i
   integer :: pol
   complex(8), intent(in) :: p(4)
   logical, intent(in),optional :: outgoing
   complex(8) :: pol_mass(4)
   complex(8) :: msq,m
   real(8) :: p0,px,py,pz, pv,pvsq
   real(8) :: ct,st,cphi,sphi
   real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0   
   complex(8), parameter :: czero = dcmplx(0d0,0d0), ci=dcmplx(0d0,1d0), cone=dcmplx(1d0,0d0)   

      p0=dreal(p(1))
      px=dreal(p(2))
      py=dreal(p(3))
      pz=dreal(p(4))

      pv=px**2 + py**2 + pz**2
      m=p0**2-pv
      m=sqrt(m)
      pv=sqrt(pv)

      if(cdabs(pv/m).lt.1d-8) then
         if(i.eq.0) then
            pol_mass(1:3)=czero
            pol_mass( 4 )=cone
            return
         endif
         ct = 1d0; st=0d0
      else
         ct= pz/pv
         st= dsqrt(dabs(1.0d0-ct**2))
      endif


      if (st .lt. 1D-15) then
         cphi=1.0d0
         sphi=0d0
      else
         cphi= px/pv/st
         sphi= py/pv/st
      endif


!     i=0 is longitudinal polarization
!     the following ifstatement distinguishes between
!     positive and negative energies
      if ( p0 .gt. 0.0d0) then
         pol=i
      else
         pol=-i
      endif

      ! -- take complex conjugate for outgoing
      if (present(outgoing)) then
         if (outgoing) pol = -pol
      endif

      if(pol.eq.-1 .or. pol.eq.1) then
         pol_mass(1)=czero
         pol_mass(2)=(ct*cphi-pol*ci*sphi)/sqrt2
         pol_mass(3)=(ct*sphi+pol*ci*cphi)/sqrt2
         pol_mass(4)=-st/sqrt2
      else if(pol.eq.0) then
         pol_mass(1)= pv/m
         pol_mass(2)= p0/m/pv*px
         pol_mass(3)= p0/m/pv*py
         pol_mass(4)= p0/m/pv*pz
      else
         print *,"wrong helicity setting in pol_mass"
         stop
      endif

   end function pol_mass




        FUNCTION vbqg_Weyl(sp,e1)
        implicit none
        complex(8), intent(in) :: e1(:)
        complex(8), intent(in) :: sp(:)
        complex(8) :: vbqg_Weyl(size(sp))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            vbqg_Weyl = (0d0,-1d0)/sqrt2*spi2_Weyl(e1,sp)

        END FUNCTION




        FUNCTION vbqq_Weyl(sp1,sp2)
        implicit none
        complex(8), intent(in) :: sp1(:), sp2(:)
        integer, parameter ::  Dv=4
        integer :: i
        complex(8) :: vbqq_Weyl(Dv)
        complex(8) :: rr, va(Dv),sp1a(size(sp1))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            va=(0d0,0d0)
            vbqq_Weyl=(0d0,0d0)

            do i=1,Dv
              if (i.eq.1) then
                va(1)=(1d0,0d0)
              else
                va(i)=(-1d0,0d0)
              endif
              sp1a=spb2_Weyl(sp1,va)

              rr=(0d0,-1d0)/sqrt2*psp1_(sp1a,sp2)
              if (i.eq.1) then
                    vbqq_Weyl = vbqq_Weyl + rr*va
                else
                    vbqq_Weyl = vbqq_Weyl - rr*va
              endif
              va(i)=(0d0,0d0)
            enddo

        END FUNCTION



        function vgq_Weyl(e1,sp)
        implicit none
        complex(8), intent(in) :: e1(:)
        complex(8), intent(in) :: sp(:)
        complex(8) :: vgq_Weyl(size(sp))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            vgq_Weyl = (0d0,-1d0)/sqrt2*spb2_Weyl(sp,e1)

        end function





      function spb2_Weyl(sp,v)
      implicit none
      complex(8), intent(in) :: sp(:),v(:)
      integer, parameter ::  Dv=4, Ds=4
      complex (8) :: spb2_Weyl(size(sp))
      complex(8) :: x0(4,4),xx(4,4),xy(4,4)
      complex(8) :: xz(4,4),x5(4,4)
      complex(8) :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,imax



      imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2

           xx(1,i) = y4
           xx(2,i) = y3
           xx(3,i) = -y2
           xx(4,i) = -y1

           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=y3
           xz(2,i)=-y4
           xz(3,i)=-y1
           xz(4,i)=y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo



           do i=1,4

           spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo


           end function





         function spi2_Weyl(v,sp)
         implicit none
         complex(8), intent(in) :: sp(:),v(:)
         complex(8) :: spi2_Weyl(size(sp))
         integer, parameter ::  Dv=4,Ds=4
         complex(8) :: x0(4,4),xx(4,4),xy(4,4)
         complex(8) :: xz(4,4),x5(4,4)
         complex(8) ::  y1,y2,y3,y4,bp,bm,cp,cm
         integer :: i,i1,i2,i3,imax


         imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2


           xx(1,i) = -y4
           xx(2,i) = -y3
           xx(3,i) = y2
           xx(4,i) = y1


           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=-y3
           xz(2,i)=y4
           xz(3,i)=y1
           xz(4,i)=-y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo


           do i=1,4

           spi2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1) -v(3)*xy(i,1)-v(4)*xz(i,1)
           enddo


           end function



          function WeylToDirac(sp)   ! unitary transformation U to convert a Weyl spinor into the Dirac representation
          implicit none              ! sp can be spinor or bar-spinor, i.e. U^dagger.sp = barsp.U
          double complex :: sp(1:4)
          double complex :: WeylToDirac(1:4)
          double precision,parameter :: SqrtFac=1d0/dsqrt(2d0)

              WeylToDirac(1) = SqrtFac*(sp(1)+sp(3))
              WeylToDirac(2) = SqrtFac*(sp(2)+sp(4))
              WeylToDirac(3) = SqrtFac*(sp(1)-sp(3))
              WeylToDirac(4) = SqrtFac*(sp(2)-sp(4))
          return
          end function




          function psp1_(sp1,sp2) result(res)
          implicit none
          complex(8), intent(in) :: sp1(:)
          complex(8), intent(in) :: sp2(:)
          complex(8) :: res

            res = sum(sp1(1:)*sp2(1:))

           end function
    
    
    
END MODULE

