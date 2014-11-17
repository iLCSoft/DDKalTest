#ifndef DDPLANARHIT_H
#define DDPLANARHIT_H

/** DDPlanarHit: User defined KalTest hit class using u and v coordinates, which provides coordinate vector as defined by the MeasLayer 
 *
 * @author S.Aplin DESY
 */

#include "KalTrackDim.h"

#include "DDVTrackHit.h"

#define DDPlanarHit_DIM 2

class DDPlanarHit : public DDVTrackHit {
  
public:
  
  /** Constructor Taking u and v coordinates and associated measurement layer, with bfield */
  DDPlanarHit(const TVMeasLayer  &ms,
               Double_t           *x,
               Double_t           *dx,
               Double_t           bfield,
               EVENT::TrackerHit* trkhit) 
  : DDVTrackHit(ms, x, dx, bfield, DDPlanarHit_DIM,trkhit)
  { /* no op */ } 
  
  // TVTrackHit's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
  
  /** Print Debug information */
  virtual void       DebugPrint(Option_t *opt = "")           const;
  
  
private:
  
  
};
#endif
