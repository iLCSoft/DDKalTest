#ifndef DDCYLINDERHIT_H
#define DDCYLINDERHIT_H

/** DDCylinderHit: User defined KalTest hit class using R and Rphi coordinates,
 *  which provides coordinate vector as defined by the MeasLayer
 *
 * @author  F.Gaede, S.Aplin DESY
 */

#include "kaltest/KalTrackDim.h"
#include "DDVTrackHit.h"

#include <iosfwd>


class DDCylinderHit : public DDVTrackHit {
  
public:
  
  
  /** Constructor Taking R and Rphi coordinates and associated measurement layer, with bfield */
 DDCylinderHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx,
	       Double_t bfield, EVENT::TrackerHit* trkhit )
   : DDVTrackHit(ms, x, dx, bfield, 2, trkhit) {}
    
  
  // TVTrackHit's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  TKalMatrix XvToMv(const TVector3 &xv, Double_t t0) const override;
  
  /** Print Debug information */
  void DebugPrint(Option_t *opt = "") const;
  void DebugPrint(Option_t *opt, Int_t nc) const override;
  void DebugPrint(std::ostream &os, Option_t *opt = "", Int_t nc = 5) const override;
  
};
#endif
