
#include "DDKalTest/DDVMeasLayer.h"
#include "DDKalTest/MaterialMap.h"

#include "kaltest/TVMeasLayer.h"  // from KalTrackLib
#include "kaltest/TKalTrack.h"    // from KalTrackLib
#include "kaltest/TVTrack.h"      // from KalTrackLib

#include "aidaTT/materialUtils.hh"

#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>

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



DDVMeasLayer::DDVMeasLayer( dd4hep::rec::ISurface* surf,
			    Double_t   Bz,
			    const Char_t    *name)  
  : 
  TVMeasLayer(MaterialMap::get( surf->innerMaterial() ), 
	      MaterialMap::get( surf->outerMaterial() ) ,
	      surf->type().isSensitive(), 
	      name),
  _Bz(Bz),
  _surf( surf) {
  
  unsigned cellID = surf->id() ;
  _cellIDs.push_back( cellID );
  
  UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ; 
  encoder.setValue(cellID);
  encoder[ UTIL::LCTrackerCellID::module() ] = 0;
  encoder[ UTIL::LCTrackerCellID::sensor() ] = 0;
  
  _layerID = encoder.lowWord();
  
}


DDVMeasLayer::DDVMeasLayer( dd4hep::rec::ISurface* surf,
			    TMaterial &min,
			    TMaterial &mout,
			    Double_t   Bz,
			    Bool_t     is_active,
			    int        cellID ,
			    const Char_t    *name)  
  : TVMeasLayer(min, mout, is_active, name),
    _Bz(Bz),
    _surf( surf) {

  _cellIDs.push_back(cellID);
  
  UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ; 
  encoder.setValue(cellID);
  encoder[ UTIL::LCTrackerCellID::module() ] = 0;
  encoder[ UTIL::LCTrackerCellID::sensor() ] = 0;
  
  _layerID = encoder.lowWord();
  
}


DDVMeasLayer::DDVMeasLayer(  dd4hep::rec::ISurface* surf,
			     TMaterial &min,
                             TMaterial &mout,
                             Double_t  Bz,
                             const std::vector<int>& cellIDs,
                             Bool_t    is_active,
                             const Char_t    *name)
  : TVMeasLayer(min, mout, is_active, name),
    _Bz(Bz),
    _cellIDs(cellIDs),
    _isMultiLayer(true),
    _surf( surf) {
  
  if (cellIDs.size() == 0 ) {
    streamlog_out(ERROR) << __FILE__ << " line " << __LINE__ << " size of cellIDs == 0" << std::endl;
  }
  
  UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
  encoder.setValue(cellIDs.at(0));
  encoder[ UTIL::LCTrackerCellID::module() ] = 0;
  encoder[ UTIL::LCTrackerCellID::sensor() ] = 0;
  
  _layerID = encoder.lowWord();

}



/** Function to compute the energy loss per path length and density.
 *  Could make this a global static function that could be modified 
 *  if needed ...
 */


// Ziegler param of proton stopping power on various materials
//   adapted from Geant4 class G4hZiegler1985p.cc

const double G4hZiegler1985p_a[92][8] = {
  {0.0091827, 0.0053496, 0.69741, 0.48493, 316.07 , 1.0143, 9329.3, 0.053989},
  {0.11393,   0.0051984, 1.0822,  0.39252, 1081.0 , 1.0645, 4068.5, 0.017699},
  {0.85837,   0.0050147, 1.6044,  0.38844, 1337.3 , 1.047,  2659.2, 0.01898},
  {0.8781,    0.0051049, 5.4232,  0.2032 , 1200.6 , 1.0211, 1401.8, 0.038529},
  {1.4608,    0.0048836, 2.338,   0.44249, 1801.3 , 1.0352, 1784.1, 0.02024},
  {3.2579,    0.0049148, 2.7156,  0.36473, 2092.2 , 1.0291, 2643.6, 0.018237},
  {0.59674,   0.0050837, 4.2073,  0.30612, 2394.2 , 1.0255, 4892.1, 0.016006},
  {0.75253,   0.0050314, 4.0824,  0.30067, 2455.8 , 1.0181, 5069.7, 0.017426},
  {1.226,     0.0051385, 3.2246,  0.32703, 2525.0 , 1.0142, 7563.6, 0.019469},
  {1.0332,    0.0051645, 3.004,   0.33889, 2338.6 ,0.99997, 6991.2, 0.021799},

  {6.0972,    0.0044292, 3.1929,  0.45763, 1363.3 , 0.95182, 2380.6, 0.081835},
  {14.013,    0.0043646, 2.2641,  0.36326, 2187.4 , 0.99098, 6264.8, 0.0462},
  {0.039001,  0.0045415, 5.5463,  0.39562, 1589.2 , 0.95316, 816.16, 0.047484},
  {2.072,     0.0044516, 3.5585,  0.53933, 1515.2 , 0.93161, 1790.3, 0.035189},
  {17.575,    0.0038346, 0.078694,1.2388,  2806.0 , 0.97284, 1037.6, 0.012879},
  {16.126,    0.0038315, 0.054164,1.3104,  2813.3 , 0.96587, 1251.4, 0.011847},
  {3.217,     0.0044579, 3.6696,  0.5091,  2734.6 , 0.96253, 2187.5, 0.016907},
  {2.0379,    0.0044775, 3.0743,  0.54773, 3505.0 , 0.96575, 1714.0, 0.011701},
  {0.74171,   0.0043051, 1.1515,  0.95083, 917.21 , 0.8782,  389.93, 0.18926},
  {9.1316,    0.0043809, 5.4611,  0.31327, 3891.8 , 0.97933, 6267.9, 0.015196},
  
  {7.2247,    0.0043718, 6.1017,  0.37511, 2829.2 , 0.95218, 6376.1, 0.020398},
  {0.147,     0.0048456, 6.3485,  0.41057, 2164.1 , 0.94028, 5292.6, 0.050263},
  {5.0611,    0.0039867, 2.6174,  0.57957, 2218.9 , 0.92361, 6323.0, 0.025669},
  {0.53267,   0.0042968, 0.39005, 1.2725,  1872.7 , 0.90776, 64.166, 0.030107},
  {0.47697,   0.0043038, 0.31452, 1.3289,  1920.5 , 0.90649, 45.576, 0.027469},
  {0.027426,  0.0035443, 0.031563,2.1755,  1919.5 , 0.90099, 23.902, 0.025363},
  {0.16383,   0.0043042, 0.073454,1.8592,  1918.4 , 0.89678, 27.61,  0.023184},
  {4.2562,    0.0043737, 1.5606,  0.72067, 1546.8 , 0.87958, 302.02, 0.040944},
  {2.3508,    0.0043237, 2.882,   0.50113, 1837.7 , 0.89992, 2377.0, 0.04965},
  {3.1095,    0.0038455, 0.11477, 1.5037,  2184.7 , 0.89309, 67.306, 0.016588},

  {15.322,     0.0040306, 0.65391, 0.67668, 3001.7 , 0.92484, 3344.2, 0.016366},
  {3.6932,    0.0044813, 8.608,   0.27638, 2982.7 , 0.9276,  3166.6, 0.030874},
  {7.1373,    0.0043134, 9.4247,  0.27937, 2725.8 , 0.91597, 3166.1, 0.025008},
  {4.8979,    0.0042937, 3.7793,  0.50004, 2824.5 , 0.91028, 1282.4, 0.017061},
  {1.3683,    0.0043024, 2.5679,  0.60822, 6907.8 , 0.9817,  628.01, 0.0068055},
  {1.8301,    0.0042983, 2.9057,  0.6038,  4744.6 , 0.94722, 936.64, 0.0092242},
  {0.42056,   0.0041169, 0.01695, 2.3616,  2252.7 , 0.89192, 39.752, 0.027757},
  {30.78,     0.0037736, 0.55813, 0.76816, 7113.2 , 0.97697, 1604.4, 0.0065268},
  {11.576,    0.0042119, 7.0244,  0.37764, 4713.5 , 0.94264, 2493.2, 0.01127},
  {6.2406,    0.0041916, 5.2701,  0.49453, 4234.6 , 0.93232, 2063.9, 0.011844},
  {0.33073,   0.0041243, 1.7246,  1.1062,  1930.2 , 0.86907, 27.416, 0.038208},
  {0.017747,  0.0041715, 0.14586, 1.7305,  1803.6 , 0.86315, 29.669, 0.032123},
  {3.7229,    0.0041768, 4.6286,  0.56769, 1678.0 , 0.86202, 3094.0, 0.06244},
  {0.13998,   0.0041329, 0.25573, 1.4241,  1919.3 , 0.86326, 72.797, 0.032235},
  {0.2859,    0.0041386, 0.31301, 1.3424,  1954.8 , 0.86175, 115.18, 0.029342},
  {0.76002,   0.0042179, 3.386,   0.76285, 1867.4 , 0.85805, 69.994, 0.036448},
  {6.3957,    0.0041935, 5.4689,  0.41378, 1712.6 , 0.85397, 18493., 0.056471},
  {3.4717,    0.0041344, 3.2337,  0.63788, 1116.4 , 0.81959, 4766.0, 0.1179},
  {2.5265,    0.0042282, 4.532,   0.53562, 1030.8 , 0.81652, 16252., 0.19722},
  {7.3683,    0.0041007, 4.6791,  0.51428, 1160.0 , 0.82454, 17956., 0.13316},

  {7.7197,    0.004388,  3.242,   0.68434, 1428.1 , 0.83389, 1786.7, 0.066512},
  {16.78,     0.0041918, 9.3198,  0.29569, 3370.9 , 0.90298, 7431.7, 0.02616},
  {4.2132,    0.0042098, 4.6753,  0.57945, 3503.9 , 0.89261, 1468.9, 0.014359},
  {4.0818,    0.004214,  4.4425,  0.58393, 3945.3 , 0.90281, 1340.5, 0.013414},
  {0.18517,   0.0036215,0.00058788,3.5315, 2931.3 , 0.88936,  26.18, 0.026393},
  {4.8248,    0.0041458, 6.0934,  0.57026, 2300.1 , 0.86359, 2980.7, 0.038679},
  {0.49857,   0.0041054, 1.9775,  0.95877, 786.55 , 0.78509,  806.6, 0.40882},
  {3.2754,    0.0042177, 5.768,   0.54054, 6631.3 , 0.94282, 744.07, 0.0083026},
  {2.9978,    0.0040901, 4.5299,  0.62025, 2161.2 , 0.85669, 1268.6, 0.043031},
  {2.8701,    0.004096,  4.2568,  0.6138,  2130.4 , 0.85235, 1704.1, 0.039385},
  
  {10.853,    0.0041149, 5.8906,  0.46834, 2857.2 , 0.87550, 3654.2, 0.029955},
  {3.6407,    0.0041782, 4.8742,  0.57861, 1267.7 , 0.82211, 3508.2, 0.24174},
  {17.645,    0.0040992, 6.5855,  0.32734, 3931.3 , 0.90754, 5156.7, 0.036278},
  {7.5309,    0.0040814, 4.9389,  0.50679, 2519.7 , 0.85819, 3314.6, 0.030514},
  {5.4742,    0.0040829, 4.897,   0.51113, 2340.1 , 0.85296, 2342.7, 0.035662},
  {4.2661,    0.0040667, 4.5032,  0.55257, 2076.4 , 0.84151, 1666.6, 0.040801},
  {6.8313,    0.0040486, 4.3987,  0.51675, 2003.0 , 0.83437, 1410.4, 0.03478},
  {1.2707,    0.0040553, 4.6295,  0.57428, 1626.3 , 0.81858, 995.68, 0.055319},
  {5.7561,    0.0040491, 4.357,   0.52496, 2207.3 , 0.83796, 1579.5, 0.027165},
  {14.127,    0.0040596, 5.8304,  0.37755, 3645.9 , 0.87823, 3411.8, 0.016392},

  {6.6948,    0.0040603, 4.9361,  0.47961, 2719.0 , 0.85249, 1885.8, 0.019713},
  {3.0619,    0.0040511, 3.5803,  0.59082, 2346.1 , 0.83713, 1222.0, 0.020072},
  {10.811,    0.0033008, 1.3767,  0.76512, 2003.7 , 0.82269, 1110.6, 0.024958},
  {2.7101,    0.0040961, 1.2289,  0.98598, 1232.4 , 0.79066, 155.42, 0.047294},
  {0.52345,   0.0040244, 1.4038,  0.8551,  1461.4 , 0.79677, 503.34, 0.036789},
  {0.4616,    0.0040203, 1.3014,  0.87043, 1473.5 , 0.79687, 443.09, 0.036301},
  {0.97814,   0.0040374, 2.0127,  0.7225,  1890.8 , 0.81747, 930.7,  0.02769},
  {3.2086,    0.0040510, 3.6658,  0.53618, 3091.2 , 0.85602, 1508.1, 0.015401},
  {2.0035,    0.0040431, 7.4882,  0.3561,  4464.3 , 0.88836, 3966.5, 0.012839},
  {15.43,     0.0039432, 1.1237,  0.70703, 4595.7 , 0.88437, 1576.5, 0.0088534},

  {3.1512,    0.0040524, 4.0996,  0.5425,  3246.3 , 0.85772, 1691.8, 0.015058},
  {7.1896,    0.0040588, 8.6927,  0.35842, 4760.6 , 0.88833, 2888.3, 0.011029},
  {9.3209,    0.0040540, 11.543,  0.32027, 4866.2 , 0.89124, 3213.4, 0.011935},
  {29.242,    0.0036195, 0.16864, 1.1226,  5688.0 , 0.89812, 1033.3, 0.0071303},
  {1.8522,    0.0039973, 3.1556,  0.65096, 3755.0 , 0.86383, 1602.0, 0.012042},
  {3.222,     0.0040041, 5.9024,  0.52678, 4040.2 , 0.86804, 1658.4, 0.011747},
  {9.3412,    0.0039661, 7.921,   0.42977, 5180.9 , 0.88773, 2173.2, 0.0092007},
  {36.183,    0.0036003, 0.58341, 0.86747, 6990.2 , 0.91082, 1417.1, 0.0062187},
  {5.9284,    0.0039695, 6.4082,  0.52122, 4619.5 , 0.88083, 2323.5, 0.011627},
  {5.2454,    0.0039744, 6.7969,  0.48542, 4586.3 , 0.87794, 2481.5, 0.011282},
  
  {33.702,    0.0036901, 0.47257, 0.89235, 5295.7 , 0.8893,  2053.3, 0.0091908},
  {2.7589,    0.0039806, 3.2092,  0.66122, 2505.4 , 0.82863, 2065.1, 0.022816}
};
  
// Ziegler param of proton stopping power on various materials
//   adapted from Geant4 class G4hZiegler1985p.cc

Double_t G4hZiegler1985p_ElectronicStoppingPower(Double_t A, Double_t Z, double mass, double mom2 ) { //double kineticEnergy) {
  double ionloss ;
  int i = int(Z) - 1 ;  // index of atom
  if(i < 0)  i = 0 ;
  if(i > 91) i = 91 ;

  Double_t bg2  = mom2 / (mass * mass);
  const float protonMassKeV = 938272.;
  double T = (sqrt( bg2 + 1. ) - 1.)*protonMassKeV; // proton KE in keV at this betagamma

  // The data and the fit from: 
  // J.F.Ziegler, J.P.Biersack, U.Littmark The Stoping and
  // Range of Ions in Solids, Vol.1, Pergamon Press, 1985
  // Proton kinetic energy for parametrisation in Ziegler's units (keV/amu)  
  double e = T ;
  if ( T < 25.0 ) e = 25.0 ;
  // universal approximation  
  double slow  = G4hZiegler1985p_a[i][0] * std::pow(e, G4hZiegler1985p_a[i][1]) + G4hZiegler1985p_a[i][2] * std::pow(e, G4hZiegler1985p_a[i][3])  ;
  double shigh = std::log( G4hZiegler1985p_a[i][6]/e + G4hZiegler1985p_a[i][7]*e ) * G4hZiegler1985p_a[i][4] / std::pow(e, G4hZiegler1985p_a[i][5]) ;
  ionloss = slow*shigh / (slow + shigh) ; 
  // low energy region
  if ( T < 25.0 ) {
    double  sLocal = 0.45 ;
    // light elements
    if(6.5 > Z) sLocal = 0.25 ;
    // semiconductors
    if(5 == i || 13 == i || 31 == i) sLocal = 0.375 ;
    ionloss *= std::pow(T/25.0, sLocal) ;
  }
  if ( ionloss < 0.0) ionloss = 0.0 ;
  // ionloss is now in units of eV / 10^15 atoms / cm2
  ionloss *= (0.6022/A);  // convert to MeV/g/cm2 (?)
  return ionloss;
}


double dedx_bethe_bloch( Double_t dnsty, Double_t A, Double_t Z , double mass, double mom2 ) {
  // -----------------------------------------
  // Bethe-Bloch eq. (Physical Review D P195.)
  // -----------------------------------------
  static const Double_t kK   = 0.307075e-3;     // [GeV*cm^2]
  static const Double_t kMe  = 0.510998902e-3;  // electron mass [GeV]
  
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
					- bg2/gm2 - 0.5*del);
  
  return dedx ;  // this is in units of MeV / g / cm2
}

double computeDEdx( const TMaterial &mat, double mass, double mom2 ){
  // get material properties
  // computes dedx using bethe-bloch or Ziegler param, depending on value of beta-gamma

  Double_t dnsty = mat.GetDensity();		 // density
  Double_t A     = mat.GetA();                   // atomic mass
  Double_t Z     = mat.GetZ();                   // atomic number
  Double_t bg2  = mom2 / (mass * mass);

  double dedx=0;
  if ( bg2 > 0.01 ) { // beta-gamma>0.1 -> use Bethe-Bloch
    dedx = dedx_bethe_bloch( dnsty, A, Z, mass, mom2 );
  } else { // low beta-gamma, use the Ziegler param.
    dedx = G4hZiegler1985p_ElectronicStoppingPower(A, Z, mass, mom2);
  }
  if (dedx<0) dedx=0;

  return dedx;
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

  double edep = 0 ;

  bool use_aidaTT = false ;
  if( use_aidaTT ){  // ----------------------------------------------
    Double_t dr     = hel.GetDrho();
    Double_t rho    = hel.GetRho();
    Double_t omega  = 1.0 / rho;
    Double_t z0     = hel.GetDz();
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
    dd4hep::rec::Vector3D p( - std::sin( phi0 ), std::cos( phi0 ) , tnl ) ;
    dd4hep::rec::Vector3D up = p.unit() ;
    
    // need to get the normal at crossing point ( should be the current helix' reference point) 
    const TVector3& piv = hel.GetPivot() ;
    dd4hep::rec::Vector3D xx( piv.X()*dd4hep::mm,piv.Y()*dd4hep::mm,piv.Z()*dd4hep::mm) ;
    const dd4hep::rec::Vector3D& n = _surf->normal(xx) ;
    
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
    
    edep = dedx * dnsty * projectedPath ;

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


  streamlog_out(DEBUG2) << " eloss LCTracker: " << edep << ", " << UTIL::LCTrackerCellID::valueString( _surf->id() )
  			<< " surface: " << *_surf 
  			<< std::endl ;
  
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
void DDVMeasLayer::CalcQms( Bool_t        /*isoutgoing*/,
			    const TVTrack &hel,
			    Double_t      /*df*/,
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
  const dd4hep::rec::IMaterial& mat_i = _surf->innerMaterial() ;
  const dd4hep::rec::IMaterial& mat_o = _surf->outerMaterial() ;
  double x_i = mat_i.radiationLength() ;
  double x_o = mat_o.radiationLength() ;
  double l_i = _surf->innerThickness() ;
  double l_o = _surf->outerThickness() ;

  Double_t x0inv = ( l_i/x_i + l_o/x_o ) / ( l_i + l_o ) ; 

  //compute path as projection of (straight) track to surface normal:
  Double_t phi0   = hel.GetPhi0();
  dd4hep::rec::Vector3D p( - std::sin( phi0 ), std::cos( phi0 ) , tnl ) ;
  dd4hep::rec::Vector3D up = p.unit() ;

  // need to get the normal at crossing point ( should be the current helix' reference point) 
  const TVector3& piv = hel.GetPivot() ;
  dd4hep::rec::Vector3D xx( piv.X()*dd4hep::mm,piv.Y()*dd4hep::mm,piv.Z()*dd4hep::mm) ;
  const dd4hep::rec::Vector3D& n = _surf->normal( xx ) ;

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
