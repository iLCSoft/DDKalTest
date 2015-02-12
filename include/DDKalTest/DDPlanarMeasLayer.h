#ifndef DDPlanarMeasLayer_h
#define DDPlanarMeasLayer_h

#include "DDVMeasLayer.h"
#include "TPlane.h"

#include "DDRec/Surface.h"

class TVTrackHit;
class TVector3 ;

/** DDPlanarMeasLayer provides a generic planar measurment for 1-dim and 2-dim hits
 *  using a DD4hep::DDRec::Surface.
 *  
 *  @author F.Gaede CERN/DESY
 *  @date Dec 2014
 *  @version $Id:$
 */
class DDPlanarMeasLayer : public DDVMeasLayer, public TPlane {
  
public:

  /// Ctor: initialize with Surface
  DDPlanarMeasLayer( DD4hep::DDRec::Surface* surf,
		     Double_t   Bz,		    
		     const Char_t  *name = "DDPlanarMeasL");
  
  
  virtual ~DDPlanarMeasLayer() {} ;
  
  // Parrent's pure virtuals that must be implemented
  
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const;
  
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  /** Get the intersection and the CellID, needed for multilayers 
   *  ( default implementation returns this cellID - overwrite for multilayers )
   */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                               TVector3 &xx,
                               Double_t &phi,
                               Int_t    &CellID,
                               Int_t     mode,
                               Double_t  eps = 1.e-8) const ; 


  virtual DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  virtual Bool_t IsOnSurface (const TVector3 &xx) const;
  
  Double_t GetSortingPolicy() const { return fSortingPolicy; }
  
protected:
  unsigned fMDim ;
  Double_t fSortingPolicy;

};

#endif
