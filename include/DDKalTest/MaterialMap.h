#ifndef MaterialMap_h
#define MaterialMap_h

#include <map>
#include <string>

#include "DDRec/IMaterial.h"

#include "streamlog/streamlog.h"

#include "TMaterial.h"

/** Singleton map for holding TMaterial objects, created 
 *  from dd4hep::rec::IMaterial objects - assumes that material names
 *  are unique.
 *
 * @author F.Gaede CERN/DESY
 * @date  20 Nov 2014
 */

class MaterialMap {

  typedef std::map<std::string,TMaterial*> MAP ;

public:
  
  static TMaterial& get( const dd4hep::rec::IMaterial& mat ){

    static MAP _map ;

    MAP::iterator it = _map.find( mat.name() ) ;

    if( it == _map.end() ) {

      TMaterial* m = new TMaterial( mat.name().c_str() , "" , 
				    mat.A(), mat.Z() , mat.density() , 
				    mat.radiationLength() , mat.interactionLength() ) ;

      _map[ mat.name()  ] = m ;

      streamlog_out( DEBUG3 ) << " -- MaterialMap add new material " << mat.name() << std::endl ; 

      return *m;
    }
    
    return *it->second ;
  }
  
protected:
  /// default c'tor
  MaterialMap();
  
};

#endif
