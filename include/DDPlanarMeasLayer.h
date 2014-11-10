#ifndef __DDPLANARMEASLAYER__
#define __DDPLANARMEASLAYER__
//*************************************************************************
//* ===================
//*  DDPlanarMeasLayer Class
//* ===================
//*
//* (Description)
//*   Planar measurement layer class used with DDPLanarTrackHit.
//* (Requires)
//*   DDVMeasLayer
//* (Provides)
//*     class DDPlanarMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*   2011/06/17  D.Kamai           Modified to handle ladder structure.
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TPlane.h"
#include "DDVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"
#include <sstream>

class TVTrackHit;

class DDPlanarMeasLayer : public DDVMeasLayer, public TPlane {
public:
  // Ctors and Dtor
  
  DDPlanarMeasLayer(TMaterial &min,
                     TMaterial &mout,
                     const TVector3  &center,
                     const TVector3  &normal,
                     Double_t   Bz,
                     Double_t   SortingPolicy,
                     Double_t   xiwidth,
                     Double_t   zetawidth,
                     Double_t   xioffset,
                     Double_t   fUOrigin,
                     Bool_t     is_active,
                     Int_t      CellID = -1,
                     const Char_t    *name = "DDPlanarMeasL");
  
  virtual ~DDPlanarMeasLayer();
  
  // Parrent's pure virtuals that must be implemented
  
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const;
  
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  virtual DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  Double_t GetSortingPolicy() const { return fSortingPolicy; }
  Double_t GetXiwidth() const { return fXiwidth; }
  Double_t GetZetawidth() const { return fZetawidth; }
  Double_t GetXioffset() const { return fXioffset; }
  
protected:
  Double_t fSortingPolicy;
  Double_t fXiwidth;
  Double_t fZetawidth;
  Double_t fXioffset; //determines how far the centre of the plane is translated in the direction positive rotation
  Double_t fUOrigin;  //determines origin of the transverse coordinate
  
};

#endif
