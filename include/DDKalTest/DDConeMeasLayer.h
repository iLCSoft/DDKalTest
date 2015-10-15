#ifndef DDConeMeasLayer_H
#define DDConeMeasLayer_H

#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCutCone.h"
#include "KalTrackDim.h"
#include "DDVMeasLayer.h"
#include "iostream"
#include "streamlog/streamlog.h"

//helper struct serving as private base
namespace{
  struct Data{
    Double_t _Z1;      // z of front face
    Double_t _R1;      // r of front face
    Double_t _Z2;      // z of back end
    Double_t _R2;      // r of back end
    Data(double z1,double r1, double z2,double r2 ) : _Z1(z1),_R1(r1),_Z2(z2),_R2(r2) {}
  } ;
}

/** DDConeMeasLayer provides a conical helper layer (no measurements )
 *  using a DDSurfaces::ISurface.
 *  
 *  @author F.Gaede DESY
 *  @date Nov 2015
 *  @version $Id:$
 */
class DDConeMeasLayer : public DDVMeasLayer, private Data, public TCutCone {
public:
  // Ctors and Dtor
  /// Constructor: initialize with Surface and B-field
  DDConeMeasLayer(DDSurfaces::ISurface* surf,
		  Double_t   Bz,
		  const Char_t    *name = "DDConeMeasL") ; 

  virtual ~DDConeMeasLayer();

  // Parrent's pure virtuals that must be implemented
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
				const TVector3   &xv) const;
   
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
   
  /** Local to Global coordinates */
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
   
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
				const TVector3   &xv,
				const TKalMatrix &dxphiada,
				TKalMatrix &H)  const;

  virtual Bool_t IsOnSurface(const TVector3 &xx) const;

  /** Convert LCIO Tracker Hit to an DDConeHit  */
  virtual DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* ) const {
      
    streamlog_out( ERROR ) << "DDConeMeasLayer::ConvertLCIOTrkHit Don't use this, it's not implemented!";
    return NULL;
  }
   
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
				       TVector3 &xx,
				       Double_t &phi,
				       Int_t    &CellID,
				       Int_t     /*mode*/,
				       Double_t  eps = 1.e-8) const {
                                           
    CellID = this->getCellIDs()[0]; // not multilayer
    return this->CalcXingPointWith(hel,xx,phi,0,eps);
                                           
                                           
  }

  /** Get sorting policy for this plane  */
  virtual double GetSortingPolicy() const { return fsortingPolicy; }


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


   
private:
  // Double_t fZ1;      // z of front face
  // Double_t fR1;      // r of front face
  // Double_t fZ2;      // z of back end
  // Double_t fR2;      // r of back end
  Double_t fsortingPolicy; // used for sorting the layers in to out

  
};

#endif
