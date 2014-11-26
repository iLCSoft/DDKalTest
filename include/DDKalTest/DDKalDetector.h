#ifndef DDKalDetector_h
#define DDKalDetector_h

#include "kaltest/TVKalDetector.h"
#include "DD4hep/LCDD.h"


/** Generic TVKalDetector that uses surfaces defined
 *  in a DD4hep::Geometry::DetElement.
 *
 * @author F.Gaede CERN/DESY
 * @date  20 Nov 2014
 */

class DDKalDetector : public TVKalDetector {

public:
  
  /** Initialize the detector from a DD4hep::GeometryDetElement */
  DDKalDetector( DD4hep::Geometry::DetElement det );
  
  
private:
  
};

#endif
