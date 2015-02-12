
#include "DDKalTest/DDKalDetector.h"
#include "DDKalTest/DDCylinderMeasLayer.h"
#include "DDKalTest/DDParallelPlanarMeasLayer.h"
#include "DDKalTest/DDDiscMeasLayer.h"
//#include "DDKalTest/DDParallelPlanarStripMeasLayer.h"

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
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

      Add( new DDCylinderMeasLayer( surf , Bz ) ) ;
    }
    
    else if( surf->type().isPlane() ){ 

      if( surf->type().isParallelToZ() ) {
      
 	Add( new DDParallelPlanarMeasLayer( surf , Bz ) ) ;

      } else if(  surf->type().isOrthogonalToZ() ){
	
	Add( new DDDiscMeasLayer( surf , Bz ) ) ;

      } else{ 

 	Add( new DDPlanarMeasLayer( surf , Bz ) ) ;
      }

    }



  }


}
