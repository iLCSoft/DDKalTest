
#include <iostream>

#include "DDKalTest/DDDiscMeasLayer.h"
#include "DDKalTest/DDPlanarHit.h"

#include "TVTrack.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

#include <EVENT/TrackerHitPlane.h>

#include "DDRec/Vector3D.h"

#include "streamlog/streamlog.h"


TKalMatrix DDDiscMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  // Calculate measurement vector (hit coordinates) from global coordinates:
  
  TKalMatrix mv(kMdim,1);
  
  mv(0,0)  = xv.X() ;
  
  
  mv(1,0)  = xv.Y() ;
  return mv;
  
}


TVector3 DDDiscMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  //  const DDPlanarHit &mv = dynamic_cast<const DDPlanarHit &>(vht);
  
  double x =   vht(0,0) ;
  double y =   vht(1,0) ;
  
  double z = GetXc().Z() ;
  
  return TVector3(x,y,z);
}

void DDDiscMeasLayer::CalcDhDa(const TVTrackHit &vht,
                                const TVector3   &xxv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const
{
  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  // where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
  
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) = dxphiada(0,i);
    H(1,i) = dxphiada(1,i) ;
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    H(1,sdim-1) = 0.;
  }
  
}

Int_t DDDiscMeasLayer::CalcXingPointWith(const TVTrack  &hel,
                                          TVector3 &xx,
                                          Double_t &phi,
                                          Int_t     mode,
                                          Double_t  eps) const{
    
  phi = 0.0;
  
  xx.SetX(0.0);
  xx.SetY(0.0);
  xx.SetZ(0.0);

  
  // check that direction has one of the correct values
  if( !( mode == 0 || mode == 1 || mode == -1) ) return -1 ;
  
  // get helix parameters
  Double_t dr     = hel.GetDrho();
  Double_t phi0   = hel.GetPhi0();  //
  Double_t kappa  = hel.GetKappa();
  Double_t rho    = hel.GetRho();
  Double_t omega  = 1.0 / rho;
  Double_t z0     = hel.GetDz();
  Double_t tanl   = hel.GetTanLambda();
  
  TVector3 ref_point = hel.GetPivot();
  
  
  //
  // Check if charge is nonzero.
  //
  
  Int_t    chg = (Int_t)TMath::Sign(1.1,kappa);
  if (!chg) {
    streamlog_out(ERROR) << ">>>> Error >>>> DDDiscMeasLayer::CalcXingPointWith" << std::endl
    << "      Kappa = 0 is invalid for a helix "          << std::endl;
    return -1;
  }
  
  const double sin_phi0 = sin(phi0); 
  const double cos_phi0 = cos(phi0); 
  
  const double x_pca = ref_point.x() + dr * cos_phi0 ; 
  const double y_pca = ref_point.y() + dr * sin_phi0 ; 
  const double z_pca = ref_point.z() + z0 ;
  
  const double z = this->GetXc().Z() ;
  // get path length to crossing point 
  
  const double s = ( z - z_pca ) / tanl ;
  
//  streamlog_out(DEBUG0) << "DDDiscMeasLayer::CalcXingPointWith "
//  << " ref_point.z()  = " << ref_point.z()
//  << " z = " << z
//  << " z0  = " << z0
//  << " z_pca  = " << z_pca
//  << " tanl  = " << tanl
//  << " z - z_pca  = " << z - z_pca
//  << std::endl;
  
//  TVector3 xx_n;
//  int cuts = TVSurface::CalcXingPointWith(hel, xx_n, phi, 0, eps);
//  streamlog_out(DEBUG0) << "DDDiscMeasLayer::CalcXingPointWith from Newton: cuts = " << cuts << " x = " << xx_n.x() << " y = "<< xx_n.y() << " z = " << xx_n.z() << " r = " << xx_n.Perp() << " phi = " << xx_n.Phi() << " dphi = " <<  phi << std::endl;

  
  phi = -omega * s;
  
  const double delta_phi_half = -phi/2.0 ;
  
  
  double x;
  double y;
  
  if( fabs(s) > FLT_MIN ){ // protect against starting on the plane

    x = x_pca - s * ( sin(delta_phi_half) / delta_phi_half ) *  sin( phi0 - delta_phi_half ) ;
    
    y = y_pca + s * ( sin(delta_phi_half) / delta_phi_half ) *  cos( phi0 - delta_phi_half ) ;

  }
  else{
    streamlog_out(DEBUG0) << "DDDiscMeasLayer::CalcXingPointWith Using PCA values " << std::endl;
    x = x_pca;
    y = y_pca;
    phi = 0;
  }
  
  
  // check if intersection with plane is within boundaries
  
  xx.SetXYZ(x, y, z);
  
  
  streamlog_out(DEBUG0) << "DDDiscMeasLayer::CalcXingPointWith            : cuts = " << (IsOnSurface(xx) ? 1 : 0) << " x = " << xx.x() << " y = "<< xx.y() << " z = " << xx.z() << " r = " << xx.Perp() << " phi = " << xx.Phi() << " dphi = " <<  phi << " s = " << s << " " << this->TVMeasLayer::GetName() << std::endl;  

  if( mode!=0 && fabs(phi)>1.e-10 ){ // (+1,-1) = (fwd,bwd)
    if( chg*phi*mode > 0){
      return 0;
    }
  }
  
  return (IsOnSurface(xx) ? 1 : 0);  
  
}


Bool_t DDDiscMeasLayer::IsOnSurface(const TVector3 &xx) const
{
    
  bool onSurface = false ;
  
  TKalMatrix mv = XvToMv(xx);
  
  // check whether the hit lies in the same plane as the surface
  if( TMath::Abs((xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() + (xx.Z()-GetXc().Z())*GetNormal().Z()) < 1e-4){
    // check whether the hit lies within the boundary of the surface 
    
    double r2 = mv(0,0) * mv(0,0) + mv(1,0) * mv(1,0) ;
    
    if(  r2 <= _rMax*_rMax && r2 >= _rMin*_rMin )
        { 
          onSurface = true ;
        }    
  }
  
  return onSurface;
  
}


DDVTrackHit* DDDiscMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  
  if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  
  dd4hep::rec::Vector3D U(1.0,plane_hit->getU()[1],plane_hit->getU()[0],dd4hep::rec::Vector3D::spherical);
  dd4hep::rec::Vector3D V(1.0,plane_hit->getV()[1],plane_hit->getV()[0],dd4hep::rec::Vector3D::spherical);
  dd4hep::rec::Vector3D X(1.0,0.0,0.0);
  dd4hep::rec::Vector3D Y(0.0,1.0,0.0);
  
  const float eps = 1.0e-07;
  // U must be the global X axis 
  if( fabs(1.0 - U.dot(X)) > eps ) {
    streamlog_out(ERROR) << "DDDiscMeasLayer: TrackerHitPlane measurment vectors U is not equal to the global X axis. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
    exit(1);
  }
  
  // V must be the global X axis 
  if( fabs(1.0 - V.dot(Y)) > eps ) {
    streamlog_out(ERROR) << "DDDiscMeasLayer: TrackerHitPlane measurment vectors V is not equal to the global Y axis. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
    exit(1);
  }
  
  
  const TVector3 hit( plane_hit->getPosition()[0], plane_hit->getPosition()[1], plane_hit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  double  x[2] ;
  double dx[2] ;
  
  x[0] = h(0, 0);
  x[1] = h(1, 0);
  
  dx[0] = plane_hit->getdU() ;
  dx[1] = plane_hit->getdV() ;
  
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG1) << "DDDiscMeasLayer::ConvertLCIOTrkHit DDPlanarHit created" 
  << " u = "  <<  x[0]
  << " v = "  <<  x[1]
  << " du = " << dx[0]
  << " dv = " << dx[1]
  << " x = " << plane_hit->getPosition()[0]
  << " y = " << plane_hit->getPosition()[1]
  << " z = " << plane_hit->getPosition()[2]
  << " onSurface = " << hit_on_surface
  << std::endl ;
  
  return hit_on_surface ? new DDPlanarHit( *this , x, dx, this->GetBz(), trkhit) : NULL; 
  
}
