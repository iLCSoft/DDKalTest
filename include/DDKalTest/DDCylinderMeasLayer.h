#ifndef DDCylinderMeasLayer_H
#define DDCylinderMeasLayer_H

#include "DDVMeasLayer.h"
#include <iostream>
#include <cmath>
//#include "streamlog/streamlog.h"

#include "DDSurfaces/ISurface.h"

/** DDCylinderMeasLayer provides a generic planar measurement for 1-dim and 2-dim hits
 *  using a DDSurfaces::ISurface.
 *  
 *  @author F.Gaede CERN/DESY
 *  @date Feb 2015
 *  @version $Id:$
 */
class DDCylinderMeasLayer : public DDVMeasLayer, public TCylinder {
  
public:
  
  /// Constructor: initialize with Surface and B-field
  DDCylinderMeasLayer(DDSurfaces::ISurface* surf,
		      Double_t   Bz,
		      const Char_t    *name = "DDCylinderMeasL") ; 
  
  
  Bool_t IsOnSurface(const TVector3 &xx) const {

    //fg: leave this code for now - we are restricted to cylinders around the z-axis
    bool z = ( _surf->type().isUnbounded() ?  true : (xx.Z() >= GetZmin() && xx.Z() <= GetZmax()) );
    bool r = std::fabs( (xx-this->GetXc()).Perp() - this->GetR() ) < 1.e-3; // for very short, very stiff tracks this can be poorly defined, so we relax this here a bit to 1 micron

    //    streamlog_out(DEBUG0) << "DDCylinderMeasLayer IsOnSurface for " << this->TVMeasLayer::GetName() << " R =  " << this->GetR() 
    //    << "  GetZmin() = " << GetZmin() << " GetZmax() = " << GetZmax()
    //    << " dr = " << std::fabs( (xx-this->GetXc()).Perp() - this->GetR() ) << " r = " << r << " z = " << z 
    //    << std::endl;

    return r && z;
  }

  /** overloaded version of CalcXingPointWith using closed solution */
  virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                     TVector3 &xx,
                                     Double_t &phi,
                                     Int_t     mode,
                                     Double_t  eps = 1.e-8) const;  

  /** overloaded version of CalcXingPointWith using closed solution */
  virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
                                     TVector3 &xx,
                                     Double_t &phi,
                                     Double_t  eps = 1.e-8) const ; 

  // Parent's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv)   const;
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv)   const 
  
  { return this->XvToMv(xv); }  


  /** Local to Global coordinates */
  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
  
  
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)    const;
  
  /** Convert LCIO Tracker Hit to an DDCylinderHit  */
  virtual DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                                       TVector3 &xx,
                                       Double_t &phi,
                                       Int_t    &CellID,
                                       Int_t     mode,
                                       Double_t  eps = 1.e-8) const {
    
    int ret = this->CalcXingPointWith(hel,xx,phi,0,eps) ;

    CellID = ( ( this->getCellIDs().size() > 1 && xx.Z() < 0. ) ?
	       this->getCellIDs()[1] :  this->getCellIDs()[0]  ) ;    // multilayer cylinder ?

    return ret ;
  }

  Double_t GetSortingPolicy() const { return fSortingPolicy; }
 
protected:
  unsigned fMDim ;
  Double_t fSortingPolicy;
  
private:
  
};
#endif
