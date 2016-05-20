
#include "DDKalTest/DDVMeasLayer.h"
#include "DDKalTest/MaterialMap.h"

#include "kaltest/TVMeasLayer.h"  // from KalTrackLib
#include "kaltest/TKalTrack.h"    // from KalTrackLib
#include "kaltest/TVTrack.h"      // from KalTrackLib

#include "aidaTT/materialUtils.hh"

#include <UTIL/BitField64.h>
#include <DDKalTest/DDKalTestConf.h>

#include <DD4hep/DD4hepUnits.h>

#include "streamlog/streamlog.h"

Bool_t   DDVMeasLayer::kActive = kTRUE;
Bool_t   DDVMeasLayer::kDummy = kFALSE;

//ClassImp(DDVMeasLayer)

namespace{
  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi)  {
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }
}



DDVMeasLayer::DDVMeasLayer( DDSurfaces::ISurface* surf,
			    Double_t   Bz,
			    const Char_t    *name)  
  : 
  TVMeasLayer(MaterialMap::get( surf->innerMaterial() ), 
	      MaterialMap::get( surf->outerMaterial() ) ,
	      surf->type().isSensitive(), 
	      name),
  _surf( surf) ,
  _Bz(Bz),
  _isMultiLayer(false) {
  
  unsigned cellID = surf->id() ;
  _cellIDs.push_back( cellID );
  
  UTIL::BitField64 encoder( DDKalTest::CellIDEncoding::instance().encoding_string() ) ; 
  encoder.setValue(cellID);
  encoder[ DDKalTest::CellIDEncoding::instance().module() ] = 0;
  encoder[ DDKalTest::CellIDEncoding::instance().sensor() ] = 0;
  
  _layerID = encoder.lowWord();
  
}


DDVMeasLayer::DDVMeasLayer( DDSurfaces::ISurface* surf,
			    TMaterial &min,
			    TMaterial &mout,
			    Double_t   Bz,
			    Bool_t     is_active,
			    int        cellID ,
			    const Char_t    *name)  
  : TVMeasLayer(min, mout, is_active, name),
    _surf( surf) ,
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


DDVMeasLayer::DDVMeasLayer(  DDSurfaces::ISurface* surf,
			     TMaterial &min,
                             TMaterial &mout,
                             Double_t  Bz,
                             const std::vector<int>& cellIDs,
                             Bool_t    is_active,
                             const Char_t    *name)
  : TVMeasLayer(min, mout, is_active, name),
    _surf( surf) ,
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



/** Function to compute the energy loss per path length and density.
 *  Could make this a global static function that could be modified 
 *  if needed ...
 */
double computeDEdx( const TMaterial &mat, double mass, double mom2 ){
  // -----------------------------------------
  // Bethe-Bloch eq. (Physical Review D P195.)
  // -----------------------------------------
  static const Double_t kK   = 0.307075e-3;     // [GeV*cm^2]
  static const Double_t kMe  = 0.510998902e-3;  // electron mass [GeV]
  
  Double_t dnsty = mat.GetDensity();		 // density
  Double_t A     = mat.GetA();                   // atomic mass
  Double_t Z     = mat.GetZ();                   // atomic number
  //Double_t I    = Z * 1.e-8;			 // mean excitation energy [GeV]
  //Double_t I    = (2.4 +Z) * 1.e-8;		 // mean excitation energy [GeV]
  Double_t I    = (9.76 * Z + 58.8 * TMath::Power(Z, -0.19)) * 1.e-9;
  Double_t hwp  = 28.816 * TMath::Sqrt(dnsty * Z/A) * 1.e-9;
  Double_t bg2  = mom2 / (mass * mass);
  Double_t gm2  = 1. + bg2;
  Double_t meM  = kMe / mass;
  Double_t x    = log10(TMath::Sqrt(bg2));
  Double_t C0   = - (2. * log(I/hwp) + 1.);
  Double_t a    = -C0/27.;
  Double_t del;
  if (x >= 3.)            del = 4.606 * x + C0;
  else if (0.<=x && x<3.) del = 4.606 * x + C0 + a * TMath::Power(3.-x, 3.);
  else                    del = 0.;
  Double_t tmax = 2.*kMe*bg2 / (1. + meM*(2.*TMath::Sqrt(gm2) + meM)); 
  Double_t dedx = kK * Z/A * gm2/bg2 * (0.5*log(2.*kMe*bg2*tmax / (I*I))
   					- bg2/gm2 - del);
  
  return dedx ;

}

//_________________________________________________________________________
// -----------------
//  GetEnergyLoss: code copied from KalTes::TVMeasLayer.cc - modified for usage of 
//                 a single surface, using inner and outer material and corresponding 
//                 thicknesses
// -----------------

Double_t DDVMeasLayer::GetEnergyLoss( Bool_t    isoutgoing,
				      const TVTrack  &hel,
				      Double_t  df) const {

  Double_t cpa    = hel.GetKappa();
  Double_t tnl    = hel.GetTanLambda(); 
  Double_t tnl2   = tnl * tnl;
  Double_t tnl21  = 1. + tnl2;
  Double_t cslinv = TMath::Sqrt(tnl21);
  Double_t mom2   = tnl21 / (cpa * cpa);
  Double_t phi0   = hel.GetPhi0();

  // For straight track, cpa is 0.
  if(!hel.IsInB()) { mom2 = hel.GetMomentum(); mom2 *= mom2; }

  static const Double_t kMpi = 0.13957018;      // pion mass [GeV]

  TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSystem::GetCurInstancePtr());
  Double_t   mass = ktp ? ktp->GetMass() : kMpi;


  //fixme: hack : use hardcoded muon mass:
  mass = 0.105658 ;

  double edep = 0 ;

  bool use_aidaTT = false ;
  if( use_aidaTT ){  // ----------------------------------------------
    Double_t dr     = hel.GetDrho();
    //Double_t kappa  = hel.GetKappa();
    Double_t rho    = hel.GetRho();
    Double_t omega  = 1.0 / rho;
    //  Double_t r      = TMath::Abs(rho);
    Double_t z0     = hel.GetDz();
    Double_t tanl   = hel.GetTanLambda();
    double d0 = - dr ;
    double phi0_lcio =  toBaseRange( phi0 + M_PI/2. );
    TVector3 ref_point = hel.GetPivot(); 
    aidaTT::trackParameters trkParam  ;
    trkParam.setTrackParameters( aidaTT::Vector5( omega/dd4hep::mm , tnl, phi0_lcio , d0*dd4hep::mm ,  z0 *dd4hep::mm )  ) ;
    trkParam.setReferencePoint( aidaTT::Vector3D( ref_point.X()*dd4hep::mm ,
						  ref_point.Y()*dd4hep::mm,
						  ref_point.Z()*dd4hep::mm ) ) ;
    double energy(0), beta(0) ;
    edep = aidaTT::computeEnergyLoss( _surf, trkParam , energy, beta, mass ) ;
    
  }else{
    
    //----- add energy loss from both sides of the surface 
    
    //--- compute projection cosine for momentum unit vector and surface normal
    //    this assumes that the track can be considered as being straight
    //    across the thickness of the surface
    //    also assumes that the reference point is near by 
    //    -fixme: check these conditions !
    DDSurfaces::Vector3D p( - std::sin( phi0 ), std::cos( phi0 ) , tnl ) ;
    DDSurfaces::Vector3D up = p.unit() ;
    
    // need to get the normal at crossing point ( should be the current helix' reference point) 
    const TVector3& piv = hel.GetPivot() ;
    DDSurfaces::Vector3D xx( piv.X()*dd4hep::mm,piv.Y()*dd4hep::mm,piv.Z()*dd4hep::mm) ;
    const DDSurfaces::Vector3D& n = _surf->normal(xx) ;
    
    Double_t cosTrk = std::fabs( up * n )  ;
    
    
    bool outerMat = true ;
    bool innerMat = ! outerMat ;
    
    const TMaterial &mat_o = GetMaterial( outerMat ) ;
    Double_t dnsty         = mat_o.GetDensity();		 // density
    Double_t dedx          = computeDEdx( mat_o, mass , mom2 ) ;
    
    Double_t projectedPath = _surf->outerThickness() ;
    
    //note: projectedPath is already in dd4hep(TGeo) units, i.e. cm !
    projectedPath /= cosTrk ; 
    
    //----  path from last step:
    Double_t path = hel.IsInB()
      ? TMath::Abs(hel.GetRho()*df)*cslinv
      : TMath::Abs(df)*cslinv;
    path /= 10. ; 
    
    // take the smaller of the complete step and the one projected to the surface
    // not sure if this is really needed (or even correct) as the surface thicknesses
    // should be small compared to the distance from the previous mesurement ...
    //-------------
    // path = ( projectedPath < path  ?  projectedPath  : path ) ; 
    
    Double_t edep = dedx * dnsty * projectedPath ;
    
    streamlog_out( DEBUG1) << "\n ** in  DDVMeasLayer::GetEnergyLoss: " 
			   << "\n outer material: " << mat_o.GetName()  
			   << "\n dedx: " << dedx 
			   << "\n path: " << path
			   << "\n projectedPath: " << projectedPath 
			   << "\n edep: " << edep
			   << "\n isoutgoing: " << isoutgoing
			   << "\n up : " << up
			   << "\n normal: " << n
			   << "\n cosTrk: " << cosTrk
			   << std::endl ;
    
    const TMaterial &mat_i    = GetMaterial( innerMat ) ;
    dnsty  = mat_i.GetDensity();
    dedx   = computeDEdx( mat_i, mass , mom2 ) ;
    
    projectedPath = _surf->innerThickness() ;
    
    //note: projectedPath is already in dd4hep(TGeo) units, i.e. cm !
    projectedPath /= cosTrk ; 
    
    // take the smaller of the complete step and the one projected to the surface
    //  path = ( projectedPath < path  ?  projectedPath  : path ) ; 
    
    edep += dedx * dnsty * projectedPath ;
    
    streamlog_out( DEBUG1) << "\n ** in  DDVMeasLayer::GetEnergyLoss: " 
			   << "\n inner material: " << mat_i.GetName()  
			   << "\n dedx: " << dedx 
			   << "\n path: " << path
			   << "\n projectedPath: " << projectedPath 
			   << "\n edep: " << edep
			   << "\n isoutgoing: " << isoutgoing
			   << "\n surface: " << *_surf
			   << "\n up : " << up
			   << "\n normal: " << n
			   << "\n cosTrk: " << cosTrk
			   << std::endl ;
  }


  // streamlog_out(DEBUG7) << " @@@@ eloss aidaTT: " << edepAidaTT << " eloss DDKaltest: " << edep 
  // 			<< " surface: " << *_surf 
  // 			<< std::endl ;
  
  //-----------------------
  // FIXME: Debug hack:
  //  edep *= 1.2 ;
  //----------------------
  
  if(!hel.IsInB()) return edep;
  
  Double_t cpaa = TMath::Sqrt(tnl21 / (mom2 + edep
				       * (edep + 2. * TMath::Sqrt(mom2 + mass * mass))));
  Double_t dcpa = TMath::Abs(cpa) - cpaa;

  static const Bool_t kForward  = kTRUE;
  static const Bool_t kBackward = kFALSE;
  Bool_t isfwd = ((cpa > 0 && df < 0) || (cpa <= 0 && df > 0)) ? kForward : kBackward;
  return isfwd ? (cpa > 0 ? dcpa : -dcpa) : (cpa > 0 ? -dcpa : dcpa);
}

//_________________________________________________________________________
// -----------------
//  CalQms
// -----------------
//    calculates process noise matrix for multiple scattering with
//    thin layer approximation.
//
void DDVMeasLayer::CalcQms( Bool_t        isoutgoing,
			    const TVTrack &hel,
			    Double_t      df,
			    TKalMatrix    &Qms) const
{
  Double_t cpa    = hel.GetKappa();
  Double_t tnl    = hel.GetTanLambda(); 
  Double_t tnl2   = tnl * tnl;
  Double_t tnl21  = 1. + tnl2;
  Double_t cpatnl = cpa * tnl;
  Double_t cslinv = TMath::Sqrt(tnl21);
  Double_t mom    = TMath::Abs(1. / cpa) * cslinv;
  if(!hel.IsInB()) mom = hel.GetMomentum();

  static const Double_t kMpi = 0.13957018; // pion mass [GeV]
  TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSystem::GetCurInstancePtr());
  Double_t   mass = ktp ? ktp->GetMass() : kMpi;
  Double_t   beta = mom / TMath::Sqrt(mom * mom + mass * mass);

  // *Calculate sigma_ms0 =============================================
  static const Double_t kMS1  = 0.0136;
  static const Double_t kMS12 = kMS1 * kMS1;
  static const Double_t kMS2  = 0.038;


  // average the X0 of the inner and outer materials:
  // could also move to c'tor and cache value ...
  const DDSurfaces::IMaterial& mat_i = _surf->innerMaterial() ;
  const DDSurfaces::IMaterial& mat_o = _surf->outerMaterial() ;
  double x_i = mat_i.radiationLength() ;
  double x_o = mat_o.radiationLength() ;
  double l_i = _surf->innerThickness() ;
  double l_o = _surf->outerThickness() ;

  Double_t x0inv = ( l_i/x_i + l_o/x_o ) / ( l_i + l_o ) ; 

  //compute path as projection of (straight) track to surface normal:
  Double_t phi0   = hel.GetPhi0();
  DDSurfaces::Vector3D p( - std::sin( phi0 ), std::cos( phi0 ) , tnl ) ;
  DDSurfaces::Vector3D up = p.unit() ;

  // need to get the normal at crossing point ( should be the current helix' reference point) 
  const TVector3& piv = hel.GetPivot() ;
  DDSurfaces::Vector3D xx( piv.X()*dd4hep::mm,piv.Y()*dd4hep::mm,piv.Z()*dd4hep::mm) ;
  const DDSurfaces::Vector3D& n = _surf->normal( xx ) ;

  Double_t cosTrk = std::fabs( up * n )  ;
  
  double path = l_i + l_o ;

  //note: projectedPath is already in dd4hep(TGeo) units, i.e. cm !
  path /= cosTrk ; 

  Double_t xl   = path * x0inv;

  // ------------------------------------------------------------------
  // Very Crude Treatment!!
  //  Double_t tmp = 1. + kMS2 * TMath::Log(TMath::Max(1.e-4, xl));
  //fg: we should not limit the path length here, as we have less surfaces...
  Double_t tmp = 1. + kMS2 * TMath::Log( xl );
  tmp /= (mom * beta);
  Double_t sgms2 = kMS12 * xl * tmp * tmp;
  // ------------------------------------------------------------------


  streamlog_out( DEBUG1 ) << " ** in  DDVMeasLayer::CalcQms: "
			  << " inner material: " << mat_i.name()  
			  << " outer material: " << mat_o.name()  << std::scientific
			  << " path: " << path
			  << " x0inv: " << x0inv
			  << " sgms2: " << sgms2
			  << " cosTrk: " << cosTrk
			  << std::endl ;


  Qms(1,1) = sgms2 * tnl21;
  Qms(2,2) = sgms2 * cpatnl * cpatnl;
  Qms(2,4) = sgms2 * cpatnl * tnl21;
  Qms(4,2) = sgms2 * cpatnl * tnl21;
  Qms(4,4) = sgms2 * tnl21  * tnl21;
}
