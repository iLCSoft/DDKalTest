#ifndef __DDVMEASLAYER__
#define __DDVMEASLAYER__

#include "TVector3.h"
#include "kaltest/TKalMatrix.h"
#include "kaltest/TCylinder.h"
#include "kaltest/TVMeasLayer.h"
#include "kaltest/TAttDrawable.h"
#include "kaltest/KalTrackDim.h"
#include "TString.h"

#include "DDSurfaces/ISurface.h"

#include <vector>

class TVTrackHit;
class TNode;
class DDVTrackHit;

namespace EVENT{
  class TrackerHit;
}

/** DDVMeasLayer: Virtual measurement layer class used by DD[X]MeasLayer Classes.
 *  Main methods are GetEnergyLoss() and CalcQms() which both use the DDSurfaces::ISurface
 *  and compute energy loss and material from the inner and outer materials assigned
 *  to this surface.
 *  The main difference to the implementation in KalTest::TVMeasLayer is that only
 *  one surface is needed for a wafer, as we project the path length of the particle
 *  to the normal of the surface using its thickness.
 * 
 * @author F. Gaede, CERN/DESY, S.Aplin DESY
 * @date Dec 2014
 * @version $Id:$
 */
class DDVMeasLayer : public TVMeasLayer {
public:
  
  static Bool_t kActive;
  static Bool_t kDummy;
  
  /** Get the layer ID */
  inline int getLayerID() const { return _layerID ; } 
  
  /** Get the Cell ID associated with this measurement layer */
  inline const std::vector<int>& getCellIDs() const { return _cellIDs ; }
  
  /** Get the number of Cell ID associated with this measurement layer */
  inline unsigned int getNCellIDs() const { return _cellIDs.size() ; }
    
  /** Get the Magnetic field at the measurement surface */
  inline Double_t GetBz() const { return _Bz; }
  
  /** Convert LCIO Tracker Hit to an DDPLanarTrackHit  */
  virtual DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const = 0 ;
  
  /** Check whether the measurement layer represents a series of detector elements */
  bool isMultilayer() const { return _isMultiLayer; } 
  
  /** Get the intersection and the CellID, needed for multilayers 
   */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                               TVector3 &xx,
                               Double_t &phi,
                               Int_t    &CellID,
                               Int_t     mode,
                               Double_t  eps = 1.e-8) const = 0 ; 
  

  /** Overwrite of energy loss method using the path length from DDRec::Surface
   *  and the inner and outer material assigned to the surface.
   */
  virtual Double_t GetEnergyLoss( Bool_t isoutgoing,
				  const TVTrack  &hel,
				  Double_t  df) const ;
  


  /** Overwrite of multiple scattering method using the path length from DDRec::Surface
   *  and the inner and outer material assigned to the surface.
   */
  virtual void CalcQms( Bool_t isoutgoing,
			const TVTrack &hel,
			Double_t  df,
			TKalMatrix  &Qms) const ;
    
  const DDSurfaces::ISurface* surface() const { return _surf ; }

protected:
  
  ///Simple c'tor that takes a surface, the b-field and a name - should be the only one !?
  DDVMeasLayer( DDSurfaces::ISurface* surf,
		Double_t   Bz,
		const Char_t    *name) ;

  DDVMeasLayer( DDSurfaces::ISurface* surf,
		TMaterial &min,
                TMaterial &mout,
                Double_t  Bz,
                Bool_t    is_active = DDVMeasLayer::kActive,
                int CellID = -1 , 
                const Char_t    *name = "DDMeasL");
  
  DDVMeasLayer( DDSurfaces::ISurface* surf,
		TMaterial &min,
	        TMaterial &mout,
                Double_t  Bz,
                const std::vector<int>& cellIDs,
                Bool_t    is_active = DDVMeasLayer::kActive,
                const Char_t    *name = "DDMeasL");
  
  
  
  Double_t _Bz ;       // Magnitude of B-Field in Z
  int _layerID ;
  std::vector<int> _cellIDs ;

  bool _isMultiLayer;

  DDSurfaces::ISurface* _surf ;
  
};

#endif
