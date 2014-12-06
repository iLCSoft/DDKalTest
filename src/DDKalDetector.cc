
#include "DDKalTest/DDKalDetector.h"
#include "DDKalTest/DDCylinderMeasLayer.h"
#include "DDKalTest/DDParallelPlanarMeasLayer.h"
#include "DDKalTest/DDParallelPlanarStripMeasLayer.h"

#include "DD4Hep/LCDD.h"
#include "DD4Hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"

#include "streamlog/streamlog.h"


DDKalDetector::DDKalDetector( DD4hep::Geometry::DetElement det ){
  
  
  streamlog_out( DEBUG4 ) << " DDKalDetector::DDKalDetector() initialize surfaces for detector " 
			  << det.name() << " type: " << det.type() << std::endl ;
  
  // --- get B field :
  double origin[3] = { 0., 0., 0. } , bfield[3] ;

  DD4hep::Geometry::OverlayedField ovField = DD4hep::Geometry::LCDD::getInstance().field() ;

  ovField.magneticField( origin , bfield ) ;

  double Bz = bfield[2] / dd4hep::tesla ;

  streamlog_out( DEBUG4 ) << " - use Bz = " << Bz << " Tesla " << std::endl ;
  //-------------

  DD4hep::DDRec::SurfaceManager ds( det ) ;
  
  const DD4hep::DDRec::SurfaceList& detSL = ds.surfaceList() ;
  
  for( DD4hep::DDRec::SurfaceList::const_iterator it = detSL.begin() ; it != detSL.end() ; ++it ){
    
    DD4hep::DDRec::Surface* surf =  *it ;
    
    // streamlog_out( DEBUG ) << " ------------------------- "
    // 			   << " surface: "  << *surf         << std::endl
    // 			   << " ------------------------- "  << std::endl ;


    if( surf->type().isCylinder() ) {
      
      // double  lhalf, x0, y0, z0, bz, dummy ;

      // Add( new DDCylinderMeasLayer( MaterialMap::get( surf->innerMaterial() ) ,
      // 				    MaterialMap::get( surf->outerMaterial()  ) ,
      // 				    surf->origin().rho() , // radius of cylinder 
      // 				    surf-> ,
      // 				    x0, y0, z0, bz, dummy,-1,"TPCInnerFCInr" ) );
    
    }
    
    if( surf->type().isPlane() && surf->type().isParallelToZ() ) {
      
      // if( surf->v().rho() > 0.001 ){

      // 	Add( new DDParallelPlanarStripMeasLayer( surf , Bz ) ) ;

      // } else {

 	Add( new DDParallelPlanarMeasLayer( surf , Bz ) ) ;
     // }


    }



  }


}
