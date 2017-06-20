#include <iostream>
#include <cmath>

#include "DDKalTest/DDPlanarMeasLayer.h"
#include "DDKalTest/DDPlanarHit.h"
//include "DDKalTest/MaterialMap.h"
#include <UTIL/LCTrackerConf.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Volumes.h"

#include "DDRec/Surface.h"
#include "DDSurfaces/Vector3D.h"

#include <EVENT/TrackerHitPlane.h>
#include <UTIL/Operators.h>

#include "streamlog/streamlog.h"


using namespace UTIL ;

DDPlanarMeasLayer::DDPlanarMeasLayer(dd4hep::rec::ISurface* surf, Double_t   Bz, const Char_t *name) :
  
  DDVMeasLayer(  surf, Bz, name ),
  // DDVMeasLayer(  surf,
  // 		 MaterialMap::get( surf->innerMaterial() ) , 
  // 		 MaterialMap::get( surf->outerMaterial() ) ,
  // 		 Bz, 
  // 		 surf->type().isSensitive() ,
  // 		 surf->id(), 
  // 		 name ),
  
  TPlane( TVector3( surf->origin()[0]/dd4hep::mm, 
		    surf->origin()[1]/dd4hep::mm , 
		    surf->origin()[2]/dd4hep::mm ) , TVector3( surf->normal() ) ),
  
  fMDim( surf->type().isMeasurement1D() ?  1 :  2 ) ,
  fSortingPolicy( 0.) {

  static int count=0 ;
  static const double epsilon=1e-6 ;
  static const double epsZ=1e-8 ;
  
  //fg: the sorting policy is used to find the measurment layers that are going to be hit by a track ...
  //    simply add an epslion to the radius in order to have a loop over all sensors in a layer
  //    NB:  this does not deal with the overlap region, e.g. in the VTX
  //         so in this region we might loose one of two hits on a layer
  //         -> to be done ...

  if( surf->type().isParallelToZ() ){
    
    fSortingPolicy = surf->origin().rho()/dd4hep::mm +  epsilon * count++ ;
    
  } else if( surf->type().isZDisk()  ){
    
    // need to compute the extend in radius of the surface 
    // do this along the direction of the origin vector 
    dd4hep::Volume vol = ((DD4hep::DDRec::Surface*)surf)->volume() ;

    // get global/local origin
    const dd4hep::rec::Vector3D& o = surf->origin() ;
    const dd4hep::rec::Vector3D& oL = ((DD4hep::DDRec::Surface*)surf)->volSurface().origin() ;

    // dd4hep::rec::Vector3D oR( o[0] , o[1] , 0 ) ; // radial direction of origin in global coordinates
    // if ( oR.trans2() < epsilon ){ // if origin has zero length use x-axis 
    //   oR.fill( 0., 1., 0 ) ;
    // }
    // this direction would have to be rotated into the local system of the volume
    // but we don't have access to the worlTransfrom matrix here, so
    // for now we use the y-axis ( works for Traps and Tubs )
    dd4hep::rec::Vector3D oR( 0. , 1. , 0 );

    double dist_r = 0. ;

    if( vol->GetShape()->Contains( oL.const_array() ) ){

      dist_r = vol->GetShape()->DistFromInside( const_cast<double*> ( oL.const_array() ) ,
						const_cast<double*> ( oR.const_array() ) ) ;

    } else {
      
      streamlog_out( WARNING ) << " ZDisk that does not contain its origin : " 
			       << UTIL::LCTrackerCellID::valueString( surf->id() ) 
			       << *surf  << std::endl ;
    }
    
    // streamlog_out( MESSAGE ) << " surface:    " << UTIL::LCTrackerCellID::valueString( surf->id() )  << *surf 
    // 			        << " vol shape : " << std::endl ; 
    // vol->GetShape()->Dump() ;
    
    double rMax = o.rho() + std::abs( dist_r ) ;


    fSortingPolicy =  rMax/dd4hep::mm + epsilon * count++ ;


  } else{

    fSortingPolicy =  surf->distance( dd4hep::rec::Vector3D( 0.,0.,0.) )/dd4hep::mm   +  epsilon * count++ ;
  }


  streamlog_out(DEBUG1) << "DDPlanarMeasLayer created" 
			<< " Layer x0 = " << this->GetXc().X() 
			<< " y0 = " << this->GetXc().Y() 
			<< " z0 = " << this->GetXc().Z() 
			<< " R = " << this->GetXc().Perp() 
			<< " phi = " << this->GetXc().Phi() 
			<< " soting policy : " << fSortingPolicy
			<< " is_active = " << surf->type().isSensitive()  
			<< " CellID = " << UTIL::LCTrackerCellID::valueString( surf->id() ) 
			<< " name = " << this->DDVMeasLayer::GetName()  
			<< std::endl ;
  }


TKalMatrix DDPlanarMeasLayer::XvToMv(const TVector3 &xv) const {
  
  TKalMatrix mv( fMDim , 1 );
  
  dd4hep::rec::Vector2D lv = _surf->globalToLocal( dd4hep::rec::Vector3D( xv.X()*dd4hep::mm ,  xv.Y()*dd4hep::mm ,  xv.Z()*dd4hep::mm ) ) ;
  
  mv(0,0) = lv[0] / dd4hep::mm ; 
  
  if( fMDim == 2 )
    mv(1,0) = lv[1] / dd4hep::mm ;
  
  streamlog_out(DEBUG0) << "\t DDPlanarMeasLayer::XvToMv: "
			<< " x = " << xv.X() 
			<< " y = " << xv.Y() 
			<< " z = " << xv.Z() 
			<< " mv(0,0) = " << mv(0,0) 
			<< " mv(1,0) = " << ( fMDim==2 ?  mv(1,0) : 0.0 ) 
			<< std::endl;

  // streamlog_out(DEBUG0) <<"\t surface : " << *_surf  << std::endl ;

  return mv;
}

TKalMatrix DDPlanarMeasLayer::XvToMv(const TVTrackHit &,
				     const TVector3   &xv) const {

  return XvToMv(xv);
}

TVector3 DDPlanarMeasLayer::HitToXv(const TVTrackHit &vht) const {
  
  dd4hep::rec::Vector3D v = ( fMDim == 2 ?
			     _surf->localToGlobal( dd4hep::rec::Vector2D (  vht(0,0)*dd4hep::mm, vht(1,0) *dd4hep::mm ) )  :
			     _surf->localToGlobal( dd4hep::rec::Vector2D (  vht(0,0)*dd4hep::mm,    0.  ) )
			     ) ;
  
  double x = v[0] / dd4hep::mm ;
  double y = v[1] / dd4hep::mm ;
  double z = v[2] / dd4hep::mm ;
  
  
  streamlog_out(DEBUG0) << "\t DDPlanarMeasLayer::HitToXv: "
			<< " vht(0,0) = " << vht(0,0)
			<< " vht(1,0) = " << ( fMDim==2 ? vht(1,0) : 0. )
			<< " x = " << x 
			<< " y = " << y 
			<< " z = " << z 
			<< std::endl;
  
  
  return TVector3(x,y,z);
}


void DDPlanarMeasLayer::CalcDhDa(const TVTrackHit &vht,
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
  

  //----------------------------------------------------
  //fixme: the derivatives should be stored either in the surface class or here
  dd4hep::rec::Vector3D u = _surf->u() ;
  dd4hep::rec::Vector3D v = _surf->v() ;
  
  double uv = u * v ;
  dd4hep::rec::Vector3D uprime = ( u - uv * v ).unit() ;
  dd4hep::rec::Vector3D vprime = ( v - uv * u ).unit() ;
  double uup = u * uprime ;
  double vvp = v * vprime ;
  
  
  dd4hep::rec::Vector3D dudx =  1./uup * uprime  ;
  dd4hep::rec::Vector3D dvdx =  1./vvp * vprime;
  
  //------------------------------------------------------
  
  streamlog_out(DEBUG0) << "\t DDPlanarMeasLayer::CalcDhDa: "
			<< " dudx = " << dudx[0]
			<< " dudy = " << dudx[1]
			<< " dudz = " << dudx[2]
			<< " dvdx = " << dvdx[0] 
			<< " dvdy = " << dvdx[1] 
			<< " dvdz = " << dvdx[2] 
			<< " xxv: " << xxv.X() <<", "  << xxv.Y() <<", "  << xxv.Z() 
			<< std::endl;


  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) =  dudx[0] * dxphiada(0,i) + dudx[1] * dxphiada(1,i) + dudx[2] * dxphiada(2,i) ;   
    
    if( fMDim == 2 )  
      H(1,i) =  dvdx[0] * dxphiada(0,i) + dvdx[1] * dxphiada(1,i) + dvdx[2] * dxphiada(2,i) ;   
    
  }

  if (sdim == 6) {
    
    H(0,sdim-1) = 0.;
    
    if( fMDim == 2 ) H(1,sdim-1) = 0.;
  }

}

Bool_t DDPlanarMeasLayer::IsOnSurface(const TVector3 &xx) const {
  
  //fg: here we ask the surface implementation for the bounds (using the volume) 
  //    this could be optimized by a faster algorithm (for simple bounds) 
  return _surf->insideBounds( dd4hep::rec::Vector3D( xx.x()*dd4hep::mm , xx.Y()*dd4hep::mm,  xx.Z()*dd4hep::mm ) )  ;
}


int DDPlanarMeasLayer::getIntersectionAndCellID(const TVTrack  &hel,
						TVector3 &xx,
						Double_t &phi,
						Int_t    &CellID,
						Int_t     mode,
						Double_t  eps ) const {
  
  CellID = this->getCellIDs()[0]; // not multilayer
  return this->CalcXingPointWith(hel,xx,phi,0,eps);
}
  
 
DDVTrackHit* DDPlanarMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  
  if( plane_hit == NULL ){
    streamlog_out( ERROR ) << " DDPlanarMeasLayer::ConvertLCIOTrkHit called with something that isn't an EVENT::TrackerHitPlane : " << *trkhit << std::endl ;
    return NULL; 
  }
  
  
  // get the measurment durections from the surface ( should of course be identical to the ones in the lcio hit ...)
  dd4hep::rec::Vector3D u = _surf->u() ;
  dd4hep::rec::Vector3D v = _surf->v() ;
  
  // dd4hep::rec::Vector3D U(1.0,plane_hit->getU()[1],plane_hit->getU()[0],dd4hep::rec::Vector3D::spherical);
  // dd4hep::rec::Vector3D V(1.0,plane_hit->getV()[1],plane_hit->getV()[0],dd4hep::rec::Vector3D::spherical);
  // dd4hep::rec::Vector3D Z(0.0,0.0,1.0);
  // streamlog_out(DEBUG1) << "DDPlanarMeasLayer::ConvertLCIOTrkHit : " 
  // 			// << "\n U : " << U  
  // 			<< "\n u : " << u 
  // 			<< "\n V : " << V 
  // 			<< "\n v : " << v
  // 			<< std::endl ; 
  
  // remember here the "position" of the hit in fact defines the origin of the plane it defines so u and v are per definition 0. 
  const TVector3 hit( plane_hit->getPosition()[0], plane_hit->getPosition()[1], plane_hit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  Double_t  x[ 2 ] ;
  Double_t dx[ 2 ] ;
  
  x[0] = h(0, 0);
  dx[0] = plane_hit->getdU() ;
  
  if( fMDim ==2 ) {
    x[1] = h(1, 0);
    dx[1] = plane_hit->getdV() ;
  }
  
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG0) << "DDPlanarMeasLayer::ConvertLCIOTrkHit creating DDPlanarHit " 
			<< *plane_hit 
			<< *_surf 
			<< " Layer R = " << this->GetXc().Perp() 
			<< " Layer phi = " << this->GetXc().Phi() 
			<< " Layer z0 = " << this->GetXc().Z() 
			<< " u = "  <<  x[0]
			<< " v = "  <<  x[1]
			<< " du = " << dx[0]
			<< " dv = " << dx[1]
			<< " x = " << plane_hit->getPosition()[0]
			<< " y = " << plane_hit->getPosition()[1]
			<< " z = " << plane_hit->getPosition()[2]
			<< " r = " << sqrt( plane_hit->getPosition()[0]*plane_hit->getPosition()[0] + plane_hit->getPosition()[1]*plane_hit->getPosition()[1])
			<< std::endl ;
  
  if( ! hit_on_surface )   
    streamlog_out( WARNING )  << "DDPlanarMeasLayer::ConvertLCIOTrkHit: hit is not on surface : " 
			      <<  *plane_hit 
			      << " cellID: " << UTIL::LCTrackerCellID::valueString( plane_hit->getCellID0() ) 
			      << std::endl ;
  
  return hit_on_surface ? new DDPlanarHit( *this , x, dx, this->GetBz(), trkhit, fMDim ) : NULL; 

}
