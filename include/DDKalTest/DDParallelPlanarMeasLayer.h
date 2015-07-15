#ifndef __DDParallelPlanarMeasLayer__
#define __DDParallelPlanarMeasLayer__

#include "DDPlanarMeasLayer.h"
#include "DDSurfaces/ISurface.h"


/** DDParallelPlanarMeasLayer: specialization of DDPlanarMeasruement layer 
 *  for planes that are parallel to the z-axis, where the crossing point
 *  with a helix can be comuted analytically. 
 *  We use the DDSurfaces::ISurface and aidaTT::trajectory for the implementation 
 *  of CalcXingPointWith(). 
 * 
 * @author F. Gaede CERN/DESY, S.Aplin DESY
 * @date Dec 2014
 * @version $Id:$
 */
class DDParallelPlanarMeasLayer : public DDPlanarMeasLayer {
  
public:
  
  /// c'tor using a DDRec::Surface
  DDParallelPlanarMeasLayer( DDSurfaces::ISurface* surf, Double_t   Bz) :  DDPlanarMeasLayer( surf , Bz )  {}

  
  // Parent's pure virtuals that must be implemented
  
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
  
  //fg: now implemented in DDVMeasLayer ...
  // // /** Get the intersection and the CellID, needed for multilayers */
  // // virtual int getIntersectionAndCellID(const TVTrack  &hel,
  // //                                      TVector3 &xx,
  // //                                      Double_t &phi,
  // //                                      Int_t    &CellID,
  // //                                      Int_t     mode,
  // //                                      Double_t  eps = 1.e-8) const {
  
  // //   CellID = this->getCellIDs()[0]; // not multilayer
  // //   return CalcXingPointWith(hel,xx,phi,0,eps);
  // // }
  
  
// protected:
  
//   Double_t _r;
//   Double_t _phi;
//   Double_t _cos_phi;
//   Double_t _sin_phi;

};

#endif

