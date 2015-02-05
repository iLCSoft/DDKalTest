#ifndef __DDDiscMeasLayer__
#define __DDDiscMeasLayer__

#include "DDParallelPlanarMeasLayer.h"


/** DDDiscMeasLayer: specialization of DDPlanarMeasruement layer 
 *  As we use the DD4hep::DDRec::Surface and aidaTT::trajectory for the implementation 
 *  of CalcXingPointWith() and the disctinction between parallel and orthogonal
 *  planes is done there, this class is currently identical to  DDParallelPlanarMeasLayer
 *  and we simply create a typedef for now.
 * 
 * @author F. Gaede CERN/DESY
 * @date Feb 2015
 * @version $Id:$
 */

typedef DDParallelPlanarMeasLayer DDDiscMeasLayer ;


#endif

