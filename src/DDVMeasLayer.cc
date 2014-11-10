
#include "DDVMeasLayer.h"

#include <UTIL/BitField64.h>
#include <DDKalTestConf.h>

#include "streamlog/streamlog.h"

Bool_t   DDVMeasLayer::kActive = kTRUE;
Bool_t   DDVMeasLayer::kDummy = kFALSE;

//ClassImp(DDVMeasLayer)

DDVMeasLayer::DDVMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t   Bz,
                             Bool_t     is_active,
                             int        cellID ,
                             const Char_t    *name)  
: TVMeasLayer(min, mout, is_active, name),
_Bz(Bz),
_isMultiLayer(false)
{
  _cellIDs.push_back(cellID);
  
  UTIL::BitField64 encoder( DDKalTest::CellIDEncoding::instance().encoding_string() ) ; 
  encoder.setValue(cellID);
  encoder[ DDKalTest::CellIDEncoding::instance().module() ] = 0;
  encoder[ DDKalTest::CellIDEncoding::instance().sensor() ] = 0;
  
  _layerID = encoder.lowWord();
  
}


DDVMeasLayer::DDVMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t  Bz,
                             const std::vector<int>& cellIDs,
                             Bool_t    is_active,
                             const Char_t    *name)
: TVMeasLayer(min, mout, is_active, name),
_Bz(Bz),
_cellIDs(cellIDs),
_isMultiLayer(true)
{
  
  if (cellIDs.size() == 0 ) {
    streamlog_out(ERROR) << __FILE__ << " line " << __LINE__ << " size of cellIDs == 0" << std::endl;
  }

 UTIL::BitField64 encoder( DDKalTest::CellIDEncoding::instance().encoding_string() ) ;
 encoder.setValue(cellIDs.at(0));
 encoder[DDKalTest::CellIDEncoding::instance().module() ] = 0;
 encoder[DDKalTest::CellIDEncoding::instance().sensor() ] = 0;

 _layerID = encoder.lowWord();
  
}



