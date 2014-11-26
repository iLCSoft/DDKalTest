#ifndef __DDParallelPlanarStripMeasLayer__
#define __DDParallelPlanarStripMeasLayer__

/** DDParallelPlanarStripMeasLayer: User defined KalTest measurement layer class 
 *
 * @author S.Aplin DESY
 */

//#include "TKalMatrix.h"
//#include "TVector3.h"
//#include "TVTrackHit.h"


#include "DDParallelPlanarMeasLayer.h"

#include "streamlog/streamlog.h"

class DDParallelPlanarStripMeasLayer : public DDParallelPlanarMeasLayer {
  
public:
  
  /// c'tor using a DDRec::Surface
  DDParallelPlanarStripMeasLayer( DD4hep::DDRec::Surface* surf, Double_t   Bz) :  
    DDParallelPlanarMeasLayer( surf , Bz )
  { 
    
    double angle = surf->v().theta() * 2. ; // FIXME: had 3.5 in SIT_SimplePlanar_geo !!!
    
    // check sign of rotation via helper vector cerated from cross product w/ z-axis
    DDSurfaces::Vector3D n  =  surf->v().cross(  DDSurfaces::Vector3D( 0, 0, 1.)  ) ;
    
    // does this point into the same hemisphere as the normal ?
    double ndot = n.dot( surf->normal()  ) ; 
    
    if( ndot < 0 )
      angle = - angle ;



    streamlog_out( DEBUG1 ) << "  DDParallelPlanarStripMeasLayer : set strip angle to " << angle*180./M_PI 
			    << " ndot : " << ndot
			    << std::endl ; 

    _stripAngle = angle ;
    
  }



  /** Constructor Taking inner and outer materials, distance and phi of plane pca to origin, B-Field, Sorting policy, plane transverse witdth and offset of centre, longitudinal width, whether the layer is sensitive, Cell ID, and an optional name */
  DDParallelPlanarStripMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t   r,
                             Double_t   phi,
                             Double_t   Bz,
                             Double_t   SortingPolicy,
                             Double_t   xiwidth,
                             Double_t   zetawidth,
                             Double_t   xioffset,
                             Double_t   zoffset,
                             Double_t   UOrigin,
                             Double_t   stripAngle,
                             Int_t      CellID = -1,
                             const Char_t    *name = "DDParallelPlanarStripMeasLayer")
  :
  DDParallelPlanarMeasLayer(min,mout,r,phi,Bz,SortingPolicy,xiwidth,zetawidth,xioffset,zoffset,UOrigin,true,CellID,name), _stripAngle(stripAngle)
  
  { /* no op */ }
  
  
  // Parent's pure virtuals that must be implemented

  TKalMatrix XvToMv(const TVector3 &xv) const;

  TKalMatrix XvToMv(const TVTrackHit &, const TVector3   &xv) const {
    return XvToMv(xv);
  }

  TVector3 HitToXv(const TVTrackHit &vht) const ;
  
  void CalcDhDa(const TVTrackHit &vht, const TVector3   &xxv, const TKalMatrix &dxphiada, TKalMatrix &H)  const;
    
  DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const;
  
private:
  
  double _stripAngle;
  
};

#endif

