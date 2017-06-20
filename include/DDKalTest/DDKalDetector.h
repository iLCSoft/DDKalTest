#ifndef DDKalDetector_h
#define DDKalDetector_h

#include "kaltest/TVKalDetector.h"
#include <DD4hep/DetElement.h>


/** Generic TVKalDetector that uses surfaces defined
 *  in a dd4hep::DetElement.
 *
 * @author F.Gaede CERN/DESY
 * @date  20 Nov 2014
 */

class DDKalDetector : public TVKalDetector {

public:
  
  /** Initialize the detector from a dd4hep::DetElement */
  DDKalDetector( dd4hep::DetElement det );
  
  
private:
  
};

#endif
