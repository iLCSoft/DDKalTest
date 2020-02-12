#include "DDKalTest/DDConeMeasLayer.h"
#include <UTIL/LCTrackerConf.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Vector3D.h"

#include "TVTrack.h"

#include "TMath.h"

 
//namespace{
//  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
//  inline double toBaseRange( double phi)  {
//    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
//    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
//    return phi ;
//  }
//}


DDConeMeasLayer::DDConeMeasLayer(dd4hep::rec::ISurface* surf,
				 Double_t   Bz,
				 const Char_t  *name ) :
  DDVMeasLayer(  surf, Bz, name ) ,
  Data(dynamic_cast<dd4hep::rec::ICone*>(surf)->z0()/dd4hep::mm ,
       dynamic_cast<dd4hep::rec::ICone*>(surf)->radius0()/dd4hep::mm,
       dynamic_cast<dd4hep::rec::ICone*>(surf)->z1()/dd4hep::mm,
       dynamic_cast<dd4hep::rec::ICone*>(surf)->radius1()/dd4hep::mm ),
  TCutCone(_R1*(_Z2-_Z1)/(_R2-_R1), 
   	   _R2*(_Z2-_Z1)/(_R2-_R1), 
   	   (_R2-_R1)/(_Z2-_Z1),
   	   0.,0.,(_R2*_Z1-_R1*_Z2)/(_R2-_R1)),
  fsortingPolicy( _R2 ){
  
}

DDConeMeasLayer::~DDConeMeasLayer() { }

TKalMatrix DDConeMeasLayer::XvToMv(const TVector3 &xxv) const {
  // Calculate hit coordinate information:
  //	mv(0,0) = r * phi 
  //     (1,0) = z
  
  TKalMatrix mv(kMdim,1);
  TVector3 xv = xxv - GetXc();
  Double_t r  = xv.Z()*GetTanA();
  
  mv(0,0)  = r * TMath::ATan2(xv.Y(), xv.X());
  mv(1,0)  = xv.Z();
  return mv;
}

TKalMatrix DDConeMeasLayer::XvToMv(const TVTrackHit &,
				   const TVector3   &xv) const {
  return XvToMv(xv);
}

TVector3 DDConeMeasLayer::HitToXv(const TVTrackHit &) const {
  streamlog_out( ERROR ) << "DDConeMeasLayer::HitToXv Don't use this, it's not implemented!";
  
  return TVector3(0.,0.,0.);
}

void DDConeMeasLayer::CalcDhDa(const TVTrackHit &,
			       const TVector3   &xxv,
			       const TKalMatrix &dxphiada,
			       TKalMatrix &H)  const {
  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  // where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //

  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);

  TVector3 xv = xxv - GetXc();
  Double_t x  = xv.X();
  Double_t y  = xv.Y();
  Double_t z  = xv.Z();
  Double_t xxyy = x * x + y * y;
  Double_t phi  = TMath::ATan2(y, x);
  Double_t tana = GetTanA();
  Double_t r    = z * GetTanA();

  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
   
  for (Int_t i=0; i<hdim; i++) {
    H(0,i) = - r * (y / xxyy) * dxphiada(0,i) 
      + r * (x / xxyy) * dxphiada(1,i)
      +     tana * phi * dxphiada(2,i);
    H(1,i) =  dxphiada(2,i);
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    H(1,sdim-1) = 0.;
  }
}

Bool_t DDConeMeasLayer::IsOnSurface(const TVector3 &xx) const {

  TVector3 xxc = xx - GetXc();
  Double_t r   = xxc.Perp();
  Double_t z   = xxc.Z();
  Double_t s   = (r - GetTanA()*z) * (r + GetTanA()*z);
  const Double_t tol = 1.e-8;

#if 0 
  std::cout << this->TVMeasLayer::GetName() << ":" << this->GetIndex() << ":" 
	    << "s=" << s << " xx=(" << xx.X() << "," << xx.Y() << "," << xx.Z() << ")" 
	    << "bool=" << (TMath::Abs(s) < tol && ((xx.Z()-_Z1)*(xx.Z()-_Z2) <= 0.)) 
	    << "_Z1=" << _Z1 << " _Z2=" << _Z2 << std::endl;
#endif

  return (TMath::Abs(s) < tol && ((xx.Z()-_Z1)*(xx.Z()-_Z2) <= 0.));
} 

Int_t DDConeMeasLayer::CalcXingPointWith(const TVTrack  &hel,
					     TVector3 &xx,
					     Double_t &phi,
					     Double_t  eps ) const {

  return CalcXingPointWith(hel,xx,phi,0,eps);
}


Int_t DDConeMeasLayer::CalcXingPointWith(const TVTrack  &hel,
					     TVector3 &xx,
					     Double_t &phi,
					     Int_t     mode,
					     Double_t  eps) const {

  Int_t result = TCutCone::CalcXingPointWith( hel, xx, phi, mode, eps) ;

  
  

  streamlog_out(DEBUG1) << "DDConeMeasLayer::CalcXingPointWith:on surface:" 
			 <<  IsOnSurface(xx) 
    //			<< "  (chg*phi*mode)<0: " <<  ((chg*phi*mode)<0)
			 << " x = " << xx.X()
			 << " y = " << xx.Y()
			 << " z = " << xx.Z()
			 << " r = " << xx.Perp()
			 << " phi = " << xx.Phi()
			 << " dphi = " <<  phi
			 << " " << this->TVMeasLayer::GetName()
			 << std::endl;
  
  streamlog_out(DEBUG)    << " Drho       : " << hel.GetDrho     () 
			  << " Phi0       : " << hel.GetPhi0     () 
			  << " Kappa      : " << hel.GetKappa    () 
			  << " Dz         : " << hel.GetDz       () 
			  << " TanLambda  : " << hel.GetTanLambda() 
			  << " Pivot      : " << hel.GetPivot().X() 
			  << "', " << hel.GetPivot().Y() 
			  << "', " << hel.GetPivot().Z()
			  << " Rho        : " << hel.GetRho      () 
			  << " PtoR       : " << hel.GetPtoR     ()  
			  << " MagField   : " << hel.GetMagField ()  
			  << std::endl;



  return result ;
}
