#ifndef DDKalTestConf_h
#define DDKalTestConf_h 1

#include <string>

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

    // void set_encoding_string( const std::string& enc_str )  { _encoding = enc_str ; }

    /// index of subdet in cellID
    int subdet() { return  _subdet ; } 

    /// index of side in cellID
    int side  () { return  _side   ; }

    /// index of layer in cellID
    int layer () { return  _layer  ; }

    /// index of module in cellID
    int module() { return  _module ; }

    /// index of sensor in cellID
    int sensor() { return  _sensor ; }

  protected:
    int _subdet ;
    int _side   ;
    int _layer ;
    int _module ;
    int _sensor ;
    std::string _encoding ;
  } ;


}

#endif
