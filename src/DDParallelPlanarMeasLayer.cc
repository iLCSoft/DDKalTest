#include "DDKalTest/DDParallelPlanarMeasLayer.h"
#include "DDKalTest/DDKalTestConf.h"

#include "TVTrack.h"

#include "aidaTT/trajectory.hh"

#include "DD4hep/DD4hepUnits.h"

#include "streamlog/streamlog.h"

namespace{
  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi)  {
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }
}



Int_t DDParallelPlanarMeasLayer::CalcXingPointWith(const TVTrack  &hel,
						   TVector3 &xx,
						   Double_t &phi,
						   Double_t  eps ) const {
  return CalcXingPointWith(hel,xx,phi,0,eps);
}


Int_t DDParallelPlanarMeasLayer::CalcXingPointWith(const TVTrack  &hel,
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
    
    streamlog_out(ERROR) << ">>>> Error >>>> DDParallelPlanarMeasLayer::CalcXingPointWith" << std::endl
			 << "      Kappa = 0 is invalid for a helix "          << std::endl;
    return -1;
  }
  
   
  // =============== use code from aidaTT for computing the intersection with the plane =======

  aidaTT::trackParameters trkParam  ;

  double d0 = - dr ;
  double phi0_lcio =  toBaseRange( phi0 + M_PI/2. );

  trkParam.setTrackParameters( aidaTT::Vector5( omega/dd4hep::mm , tanl, phi0_lcio , d0*dd4hep::mm ,  z0 *dd4hep::mm )  ) ; 

  //order defined in ./helpers/utilities.cc
  //  L3 type: [ Omega, tan(lambda), phi_0, d_0, z_0 ]
  trkParam.setReferencePoint( aidaTT::Vector3D( ref_point.X()*dd4hep::mm , 
						ref_point.Y()*dd4hep::mm, 
						ref_point.Z()*dd4hep::mm ) ) ;

  aidaTT::trajectory traj( trkParam , 0 ) ;
  
  double s = 0. ;
  aidaTT::Vector3D xxV3 ;
  bool foundIntersect = traj._calculateIntersectionWithSurface( _surf , s , ( aidaTT::Vector2D*) 0 , &xxV3 );
  
  
  if( foundIntersect ){
    
    s /= dd4hep::mm ; 

    xx.SetXYZ( xxV3[0]/dd4hep::mm , xxV3[1]/dd4hep::mm,  xxV3[2]/dd4hep::mm) ;   
    
    streamlog_out( DEBUG0 ) << " ++++  intersection found for surface : " << DDKalTest::CellIDEncoding::valueString(_surf->id()) << std::endl 
     			    << "       at s = " << s 
     			    << "       xx   = ( " << xx.X() << ", " << xx.Y() << ", " << xx.Z() << ") " << std::endl 
        		    << " track parameters: " << trkParam 
			    << " mode: " << mode
     			    <<  std::endl ;
    
    phi = -s * omega ; 
    
  } else {

    streamlog_out( DEBUG0 ) << " ++++ no intersection found for surface : " << DDKalTest::CellIDEncoding::valueString(_surf->id()) << std::endl
			    << " track parameters: " << trkParam 
			    << " mode : " << mode
			    << std::endl ;

    return 0 ;
  }
 

  //=============================================================================================

  streamlog_out(DEBUG1) << "DDParallelPlanarMeasLayer::CalcXingPointWith:on surface:" <<  IsOnSurface(xx) 
			<< "  (chg*phi*mode)<0: " <<  ((chg*phi*mode)<0)
			<< " x = " << xx.X()
			<< " y = " << xx.Y()
			<< " z = " << xx.Z()
			<< " r = " << xx.Perp()
			<< " phi = " << xx.Phi()
			<< " dphi = " <<  phi
			<< " " << this->TVMeasLayer::GetName() 
			<< std::endl;
  
  if( mode!=0 && fabs(phi)>1.e-10){ // (+1,-1) = (fwd,bwd)
    if( chg*phi*mode > 0){
      return 0;
    }
  }
    
  return (IsOnSurface(xx) ? 1 : 0);
}




