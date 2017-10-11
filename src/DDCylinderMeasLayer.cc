#include "TKalTrack.h" 

#include "DDKalTest/DDCylinderMeasLayer.h"
#include "DDKalTest/DDCylinderHit.h"
#include <UTIL/LCTrackerConf.h>

#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <UTIL/Operators.h>

#include "DD4hep/DD4hepUnits.h"
#include "DDSurfaces/Vector3D.h"

#include "aidaTT/trajectory.hh"

#include "streamlog/streamlog.h"

#include "TMath.h"
#include <cmath>

namespace{
  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi)  {
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }
}

using namespace UTIL ;

DDCylinderMeasLayer::DDCylinderMeasLayer(DDSurfaces::ISurface* surf,
					 Double_t   Bz,
					 const Char_t  *name ) :
  DDVMeasLayer(  surf, Bz, name ) ,
  
  TCylinder(  dynamic_cast<DDSurfaces::ICylinder*>(surf)->radius()/dd4hep::mm , 
	      surf->length_along_v()/dd4hep::mm / 2. , 
	      dynamic_cast<DDSurfaces::ICylinder*>(surf)->center().x()/dd4hep::mm, 
	      dynamic_cast<DDSurfaces::ICylinder*>(surf)->center().y()/dd4hep::mm , 
	      dynamic_cast<DDSurfaces::ICylinder*>(surf)->center().z()/dd4hep::mm ),
  
  fSortingPolicy(0.),
  
  fMDim( surf->type().isMeasurement1D() ?  1 :  2 )  {
  
  static double epsilon=1e-4 ;

  UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
  encoder.setValue( surf->id() );

  int side = encoder[ UTIL::LCTrackerCellID::side() ] ;

  // if we have an unbounded surface and the side is set to one, we also add the cellID for the other side (multilayer)
  if( side == 1 && surf->type().isUnbounded() ){
    encoder[ UTIL::LCTrackerCellID::side() ] = -1 ;
    _cellIDs.push_back( encoder.lowWord() ) ;
  }

  fSortingPolicy = dynamic_cast<DDSurfaces::ICylinder*>(surf)->radius()/dd4hep::mm + side * epsilon ;

  // assumptions made here: the cylinder runs parallel to z and v ...
  
  streamlog_out(DEBUG1) << "DDCylinderMeasLayer created" 
			<< " Layer x0 = " << this->GetXc().X() 
			<< " y0 = " << this->GetXc().Y() 
			<< " z0 = " << this->GetXc().Z() 
			<< " R = " << this->GetR() 
			<< " phi = " << this->GetXc().Phi() 
    			<< " sorting policy : " << GetSortingPolicy()
			<< " is_active = " << surf->type().isSensitive()  
			<< " CellID = " << UTIL::LCTrackerCellID::valueString( surf->id() ) 
			<< " name = " << this->DDVMeasLayer::GetName()  
			<< std::endl ;

  // for a cylindrical layer we also set the side to 0 in the layerId 
  // ( in the case the cylinder is split between forward and backward )
  encoder[ UTIL::LCTrackerCellID::side() ] = 0;
  encoder[ UTIL::LCTrackerCellID::module() ] = 0;
  encoder[ UTIL::LCTrackerCellID::sensor() ] = 0;
  
  _layerID = encoder.lowWord();

}


/** Global to Local coordinates */

TKalMatrix DDCylinderMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  // // Calculate hit coordinate information:
  // //   mv(0, 0) = r * phi
  // //     (1, 0) = drift distance
  // // account for cylinder not centered at x=0.0, y=0.0
  // TVector3 xxv = xv - GetXc();
  // Double_t phi = TMath::ATan2(xxv.Y(), xxv.X());
  // // bring phi back into +/- Pi range
  // static Double_t kPi    = TMath::Pi();
  // static Double_t kTwoPi = 2 * kPi;
  // while (phi < -kPi) phi += kTwoPi;
  // while (phi >  kPi) phi -= kTwoPi;
  // TKalMatrix mv(kMdim, 1);
  // mv(0, 0) = GetR() * phi;
  // mv(1, 0) = xxv.Z();

  TKalMatrix mv( fMDim , 1 );
  
  DDSurfaces::Vector2D lv = _surf->globalToLocal( DDSurfaces::Vector3D( xv.X()*dd4hep::mm ,  
									xv.Y()*dd4hep::mm ,  
									xv.Z()*dd4hep::mm ) ) ;
  
  mv(0,0) = lv[0] / dd4hep::mm ; 
  
  if( fMDim == 2 )
    mv(1,0) = lv[1] / dd4hep::mm ;
  
  streamlog_out(DEBUG0) << "\t DDCylinderMeasLayer::XvToMv: "
			<< " x = " << xv.X() 
			<< " y = " << xv.Y() 
			<< " z = " << xv.Z() 
			<< " mv(0,0) = " << mv(0,0) 
			<< " mv(1,0) = " << ( fMDim==2 ?  mv(1,0) : 0.0 ) 
			// << " old code : mv(0,0) = " << GetR() * phi
			// << " mv(1,0) = " <<  xxv.Z()
			<< std::endl;

  // streamlog_out(DEBUG0) <<"\t surface : " << *_surf  << std::endl ;
  
  return mv;
}


/** Local to Global coordinates */

TVector3 DDCylinderMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  
  // Double_t phi = vht(0, 0) / GetR() ;
  // Double_t z0   = vht(1, 0);
  // // account for cylinder not centered at x=0.0, y=0.0
  // Double_t x0   = GetR() * TMath::Cos(phi) + GetXc().X();
  // Double_t y0   = GetR() * TMath::Sin(phi) + GetXc().Y();
  // return TVector3(x, y, z);

  DDSurfaces::Vector3D v = ( fMDim == 2 ? 
			     _surf->localToGlobal( DDSurfaces::Vector2D (  vht(0,0)*dd4hep::mm, vht(1,0) *dd4hep::mm ) )  :
			     _surf->localToGlobal( DDSurfaces::Vector2D (  vht(0,0)*dd4hep::mm,    0.  ) ) 
			     ) ;
  
  double x = v[0] / dd4hep::mm ;
  double y = v[1] / dd4hep::mm ;
  double z = v[2] / dd4hep::mm ;
  
  
  streamlog_out(DEBUG0) << "\t DDCylinderMeasLayer::HitToXv: "
			<< " vht(0,0) = " << vht(0,0)
			<< " vht(1,0) = " << ( fMDim==2 ? vht(1,0) : 0. )
			<< " x = " << x 
			<< " y = " << y 
			<< " z = " << z 
			// << " x_old = " << x0 
			// << " y_old = " << y0 
			// << " z_old = " << z0 
			<< std::endl;
  
  
  return TVector3(x,y,z);



}


/** Calculate Projector Matrix */

void DDCylinderMeasLayer::CalcDhDa(const TVTrackHit &vht, // tracker hit not used here
                                    const TVector3   &xxv,
                                    const TKalMatrix &dxphiada,
                                    TKalMatrix &H) const
{
  
  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  //  where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5, sdim - 1);
  

#if 0 // original code from KalDet --------
  // account for cylinder not centered at x=0.0, y=0.0
  TVector3 xxvc = xxv - GetXc();
  
  Double_t xv   = xxvc.X();
  Double_t yv   = xxvc.Y();
  Double_t xxyy = xv * xv + yv * yv;
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
  
  for (Int_t i = 0; i < hdim; i++) {
    H(0, i)  = - (yv / xxyy) * dxphiada(0, i)
    + (xv / xxyy) * dxphiada(1, i);
    H(0, i) *= GetR();
    
    H(1, i)  = dxphiada(2, i);
  }
  
  if (sdim == 6) {
    H(0, sdim - 1) = 0.;
  }
  
#else // ---- new code using surfaces 

  DDSurfaces::Vector3D u = _surf->u(  DDSurfaces::Vector3D( xxv.x(),xxv.Y(),xxv.Z() ) ) ;
  DDSurfaces::Vector3D v = _surf->v(  DDSurfaces::Vector3D( xxv.x(),xxv.Y(),xxv.Z() ) ) ;
  
  double uv = u * v ;
  DDSurfaces::Vector3D uprime = ( u - uv * v ).unit() ; 
  DDSurfaces::Vector3D vprime = ( v - uv * u ).unit() ; 
  double uup = u * uprime ;
  double vvp = v * vprime ;
  
  DDSurfaces::Vector3D dudx =  1./uup * uprime  ;
  DDSurfaces::Vector3D dvdx =  1./vvp * vprime;
  
  //------------------------------------------------------
  
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) =  dudx[0] * dxphiada(0,i) + dudx[1] * dxphiada(1,i) + dudx[2] * dxphiada(2,i) ;   
    
    if( fMDim == 2 )  
      H(1,i) =  dvdx[0] * dxphiada(0,i) + dvdx[1] * dxphiada(1,i) + dvdx[2] * dxphiada(2,i) ;   
    
  }
  
  if (sdim == 6) {
    
    H(0,sdim-1) = 0.;
    
    if( fMDim == 2 ) H(1,sdim-1) = 0.;
  }
  
#endif
}


/** Convert LCIO Tracker Hit to an DDCylinderHit  */

DDVTrackHit* DDCylinderMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  if ( ! trkhit) {
    streamlog_out(ERROR) << "DDCylinderMeasLayer::ConvertLCIOTrkHit trkhit pointer is NULL" << std::endl;
    return NULL;
  }

  const TVector3 hit( trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  Double_t  x[2] ;
  Double_t dx[2] ;
  
  x[0] = h(0, 0);
  x[1] = h(1, 0);
  
  
  EVENT::TrackerHitZCylinder* cylinder_hit = dynamic_cast<EVENT::TrackerHitZCylinder*>( trkhit ) ;
  
  if(cylinder_hit){
    // convert errors
    dx[0] = cylinder_hit->getdRPhi();
    dx[1] = cylinder_hit->getdZ();
  }
  else {
    // convert errors
    dx[0] = sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]) ;
    dx[1] = sqrt(trkhit->getCovMatrix()[5]) ; 
  }
  
    
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG1) << "DDCylinderMeasLayer::ConvertLCIOTrkHit DDCylinderHit created" 
  << " R = " << hit.Perp()
  << " Layer R = " << this->GetR() 
  << " RPhi = "  <<  x[0]
  << " Z = "     <<  x[1]
  << " dRPhi = " << dx[0]
  << " dZ = "    << dx[1]
  << " x = " << trkhit->getPosition()[0]
  << " y = " << trkhit->getPosition()[1]
  << " z = " << trkhit->getPosition()[2]
  << " onSurface = " << hit_on_surface
  << std::endl ;  
  
  // we sometimes get hits which are not on the surface when running with ddsim and the 'old' digitizers
  // for now just return the hit anyways ...
  if( ! hit_on_surface ){
    streamlog_out(DEBUG9) << " DDCylinderMeasLayer::ConvertLCIOTrkHit: hit is not on surface: " << *trkhit << std::endl ;
  }
  return  new DDCylinderHit( *this , x, dx, this->GetBz(), trkhit );
 
  //  return hit_on_surface ? new DDCylinderHit( *this , x, dx, this->GetBz(), trkhit) : NULL; 
  
}

Int_t DDCylinderMeasLayer::CalcXingPointWith(const TVTrack  &hel,
					     TVector3 &xx,
					     Double_t &phi,
					     Double_t  eps ) const {

  return CalcXingPointWith(hel,xx,phi,0,eps);
}


Int_t DDCylinderMeasLayer::CalcXingPointWith(const TVTrack  &hel,
					     TVector3 &xx,
					     Double_t &phi,
					     Int_t     mode,
					     Double_t  eps) const {
  
  // check that direction has one of the correct values
  if( !( mode == 0 || mode == 1 || mode == -1) ) return -1 ;

 //fixme: currently the mode is ignored  
  //       this might cause issues in track state propagation 
  //       for curlers -> to be investigated ...   


  // This assumes nonzero B field.
  //
  // Copy helix parameters to local variables.
  //
  
  Double_t dr     = hel.GetDrho();
  Double_t phi0   = hel.GetPhi0(); //
  Double_t kappa  = hel.GetKappa();
  Double_t rho    = hel.GetRho();
  Double_t omega  = 1.0 / rho;
  Double_t r      = TMath::Abs(rho);
  Double_t z0     = hel.GetDz();
  Double_t tanl   = hel.GetTanLambda();
  
  TVector3 ref_point = hel.GetPivot();
  

  // ---  Check if charge is nonzero.
  Int_t    chg = (Int_t)TMath::Sign(1.1,kappa);

  if (!chg) {
    
    streamlog_out(ERROR) << ">>>> Error >>>> DDCylinderMeasLayer::CalcXingPointWith" << std::endl
			 << "      Kappa = 0 is invalid for a helix "          << std::endl;
    return -1;
  }
  
   
  // =============== use code from aidaTT for computing the intersection with the plane =======

  double d0 = - dr ;
  double phi0_lcio =  toBaseRange( phi0 + M_PI/2. );

  // // aidaTT::trackParameters trkParam  ;

  // // trkParam.setTrackParameters( aidaTT::Vector5( omega/dd4hep::mm , tanl, phi0_lcio , d0*dd4hep::mm ,  z0 *dd4hep::mm )  ) ; 

  // // //order defined in ./helpers/utilities.cc
  // // //  L3 type: [ Omega, tan(lambda), phi_0, d_0, z_0 ]
  // // trkParam.setReferencePoint( aidaTT::Vector3D( ref_point.X()*dd4hep::mm , 
  // // 						ref_point.Y()*dd4hep::mm, 
  // // 						ref_point.Z()*dd4hep::mm ) ) ;

  // // aidaTT::trajectory traj( trkParam , 0 ) ;
  
  // // double s = 0. ;
  // // aidaTT::Vector3D xxV3 ;
  // // bool foundIntersect = traj._calculateIntersectionWithSurface( _surf , s , ( aidaTT::Vector2D*) 0 , &xxV3 );
  
  double s = 0. ;
  aidaTT::Vector3D xxV3 ;

  aidaTT::Vector5 hp( omega/dd4hep::mm , tanl, phi0_lcio , d0*dd4hep::mm ,  z0 *dd4hep::mm ) ;
  
  aidaTT::Vector3D rp( ref_point.X()*dd4hep::mm , ref_point.Y()*dd4hep::mm, ref_point.Z()*dd4hep::mm ) ;

  bool foundIntersect = aidaTT::intersectWithZCylinder( _surf, hp, rp, s, xxV3, mode , true ) ;
  

  if( foundIntersect ){
    
    s /= dd4hep::mm ; 

    xx.SetXYZ( xxV3[0]/dd4hep::mm , xxV3[1]/dd4hep::mm,  xxV3[2]/dd4hep::mm) ;   

    streamlog_out( DEBUG3 ) << " ++++  intersection found for surface : " << UTIL::LCTrackerCellID::valueString(_surf->id()) << std::endl
			    << "       at s = " << s
			    << "       xx   = ( " << xx.X() << ", " << xx.Y() << ", " << xx.Z() << ") " << std::endl
			    << " track parameters: " <<  aidaTT::trackParameters( hp, rp )
			    << " mode: " << mode
			    <<  std::endl ;
    
    phi = -s * omega ; 
    
  } else {

    streamlog_out( DEBUG3 ) << " ++++ no intersection found for surface : " << UTIL::LCTrackerCellID::valueString(_surf->id()) << std::endl
			    << " track parameters: " <<  aidaTT::trackParameters( hp, rp ) 
  			    << " mode : " << mode
  			    << std::endl ;

    return 0 ;
  }
 

  //=============================================================================================

  streamlog_out(DEBUG1) << "DDCylinderMeasLayer::CalcXingPointWith:on surface:" <<  IsOnSurface(xx) 
  			<< "  (chg*phi*mode)<0: " <<  ((chg*phi*mode)<0)
  			<< " x = " << xx.X()
  			<< " y = " << xx.Y()
  			<< " z = " << xx.Z()
  			<< " r = " << xx.Perp()
  			<< " phi = " << xx.Phi()
  			<< " dphi = " <<  phi
  			<< " " << this->TVMeasLayer::GetName() 
  			<< std::endl;
  
    
#if 0  // DEBUG compare w/ KalTest code for crossing points

  TVector3 xxDeb ;

  // this code seems to be quite a bit faster  - need to check  !!!!!
  TCylinder::CalcXingPointWith( hel,xxDeb,phi,mode,eps) ;

   streamlog_out(DEBUG1) << "TCylinder::CalcXingPointWith:on surface:          " <<  IsOnSurface(xxDeb) 
   			<< "  (chg*phi*mode)<0: " <<  ((chg*phi*mode)<0)
   			<< " x = " << xxDeb.X()
   			<< " y = " << xxDeb.Y()
   			<< " z = " << xxDeb.Z()
   			<< " r = " << xxDeb.Perp()
   			<< " phi = " << xxDeb.Phi()
   			<< " dphi = " <<  phi
   			<< " " << this->TVMeasLayer::GetName() 
   			<< std::endl;

   streamlog_out( DEBUG ) << "DDCylinderMeasLayer::CalcXingPointWith: point on trajectory at given s (aidaTT ): " << traj.pointAt( s * dd4hep::mm ) << std::endl ;
   streamlog_out( DEBUG ) << "DDCylinderMeasLayer::CalcXingPointWith: point on trajectory at given s (KalTest): " << traj.pointAt( - phi / omega * dd4hep::mm ) << std::endl ;

   streamlog_out( DEBUG ) << "DDCylinderMeasLayer::CalcXingPointWith: distance to surface (aidaTT ): " <<  _surf->distance( xxV3 ) << std::endl ;

   streamlog_out( DEBUG ) << "DDCylinderMeasLayer::CalcXingPointWith: distance to surface (KalTest): " <<  _surf->distance( aidaTT::Vector3D( xxDeb.X() *dd4hep::mm,  xxDeb.Y()*dd4hep::mm ,   xxDeb.Z()*dd4hep::mm  )  ) << std::endl ;

//  xx = xxDeb ;
#endif
 

  if( mode!=0 && fabs(phi)>1.e-10){ // (+1,-1) = (fwd,bwd)
    if( chg*phi*mode > 0){
      return 0;
    }
  }



  //  return TCylinder::CalcXingPointWith( hel,xx,phi,mode,eps) ;

  return (IsOnSurface(xx) ? 1 : 0);
}
