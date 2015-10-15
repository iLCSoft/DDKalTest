#ifndef DDKalTestConf_h
#define DDKalTestConf_h 1

#include <string>
#include "UTIL/BitField64.h"

namespace DDKalTest{


  /** Singleton Helper class for dealing with the encoding of the cellIDs.
   *  Might be replaced by a more generic approach for encoding/decoding
   *  the cellID encoding.
   *  
   *  @author F.Gaede, CERN/DESY
   *  @date 10 Nov, 2014
   */
  class CellIDEncoding{
    
  public:

    static CellIDEncoding& instance() { 
      static CellIDEncoding _me ;
      return _me ;
    } 

    /// c'tor initialize the encoding string with the 'canonical' encoding 
    CellIDEncoding() : _encoding("subdet:5,side:-2,layer:9,module:8,sensor:8" ),
		       _subdet (0),
		       _side   (1),
		       _layer  (2),
		       _module (3),
		       _sensor (4){
    }

    // get the current encoding string
    const std::string& encoding_string() { return _encoding ; }

    /// return a string with the details of the given id:
    static std::string valueString( unsigned val ){
      UTIL::BitField64 encoder( instance().encoding_string() ) ;
      encoder.setValue( val ) ;
      return encoder.valueString() ;
    }

    // void set_encoding_string( const std::string& enc_str )  { _encoding = enc_str ; }

    /// index of subdet in cellID
    int subdet() const { return  _subdet ; } 

    /// index of side in cellID
    int side  () const { return  _side   ; }

    /// index of layer in cellID
    int layer () const { return  _layer  ; }

    /// index of module in cellID
    int module() const { return  _module ; }

    /// index of sensor in cellID
    int sensor() const { return  _sensor ; }

  protected:
    std::string _encoding ;
    const int _subdet ;
    const int _side   ;
    const int _layer ;
    const int _module ;
    const int _sensor ;
  } ;

}

#endif
