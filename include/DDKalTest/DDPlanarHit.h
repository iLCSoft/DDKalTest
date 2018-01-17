#ifndef DDPlanarHit_H
#define DDPlanarHit_H

//#include "KalTrackDim.h"
#include "DDVTrackHit.h"

/** DDPlanarHit: generic KalTest hit class for planar measurement surfaces using u and v coordinates.
 *  Can be used for 2-dim and 1-dim hits. 
 *
 * @author F.Gaede, CERN/DESY, S.Aplin DESY
 * @version $Id:$
 */
class DDPlanarHit : public DDVTrackHit {
  
public:
  
  /** Constructor Taking u and v coordinates and associated measurement layer, with bfield  and the dimension of the hit (1 for strip hits) */
  DDPlanarHit(const TVMeasLayer  &ms,
	      Double_t           *x,
	      Double_t           *dx,
	      Double_t           bfield,
	      const EVENT::TrackerHit* trkhit,
	      unsigned dimension ) 
    : DDVTrackHit(ms, x, dx, bfield, dimension ,trkhit)
  { /* no op */ } 
  
  // TVTrackHit's pure virtuals that must be implemented

  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
  
  /** Print Debug information */
  virtual void DebugPrint(Option_t *opt = "") const ;
  
  
private:
  
  
};
#endif
