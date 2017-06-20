// $Id:$
//====================================================================
//  AIDA Detector description implementation for LCD
//--------------------------------------------------------------------
//
//  Simple program to dump the complete DetElement hierarchy
// 
//  Author     : F.Gaede, CERN/DESY
//  Date       : 07 Nov 2014
//====================================================================
#include "DD4hep/DD4hepUnits.h"

#include "DDRec/Surface.h"
#include "DDRec/DetectorSurfaces.h"
#include "DDRec/SurfaceHelper.h"

#include <list>


using namespace std ;

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::SurfaceHelper;
using dd4hep::rec::SurfaceList;
using dd4hep::rec::ISurface;

//=============================================================================

int main(int argc, char** argv ){
    
  if( argc != 2 ) {
    std::cout << " usage: dumpdetector compact.xml " << std::endl ;
    exit(1) ;
  }
  
  std::string inFile =  argv[1] ;

  Detector& theDetector = Detector::getInstance();

  theDetector.fromCompact( inFile );

  DetElement world = theDetector.world() ;

  std::cout << "############################################################################### "  << std::endl  << std::endl  ;
	

  //------------------ breadth first tree traversal ---------
  std::list< DetElement > dets ;
  std::list< DetElement > daugs ; 
  std::list< DetElement > gdaugs ; 
  daugs.push_back( world ) ;
  while( ! daugs.empty() ) {
    //   for( std::list< DetElement >::iterator li=daugs.begin() ; li != daugs.end() ; ++li ){
    //     DetElement dau = *li ;
    //     DetElement::Children chMap = dau.children() ;
    //     for ( DetElement::Children::const_iterator it=chMap.begin() ; it != chMap.end() ; ++it ){
    //     	DetElement de = (*it).second ;
    //     	gdaugs.push_back( de ) ;
    //     }  
    //   }
    dets.splice( dets.end() , daugs ) ;
    daugs.splice( daugs.end() , gdaugs ) ;
  }
  //------------------ end tree traversal ---------
  

  for (  std::list< DetElement >::const_iterator it=dets.begin() ; it != dets.end() ; ++it ){

    DetElement de = (*it) ;
    
    SurfaceHelper surfMan( de ) ;

    const SurfaceList& sL = surfMan.surfaceList() ;

    std::cout << "DetElement: " << de.name() <<  " surfaces : " <<  ( sL.empty() ? 0 : sL.size()  ) << std::endl ;

    for( SurfaceList::const_iterator sit = sL.begin() ; sit != sL.end() ; ++sit ){
    
      ISurface* surf =  *sit ;
    
      std::cout << " ------------------------- " 
		<< " surface: "  << *surf         << std::endl
		<< " ------------------------- "  << std::endl ;

    }
  }

  std::cout << "############################################################################### "  << std::endl  << std::endl  ;


  return 0;
}

//=============================================================================
