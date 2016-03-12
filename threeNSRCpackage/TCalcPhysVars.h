//  TCalcPhysVars.h
//  Created by Erez Cohen on 8/10/15.
#include "TDataMining.h"


#define MyHomePath "/Users/erezcohen/Desktop"
#define MyDMPath Form("%s/DataMining",MyHomePath)
#define MyDMAnaPath Form("%s/AnaFiles",MyDMPath)
#define MyDMSchemedPath Form("%s/3NSRC/SchemedData",MyDMPath)
#define RAD2DEG TMath::RadToDeg()


using namespace std;
# define Me 0.000511
# define Mp 0.938
# define Mp2 Mp*Mp
# define Mn 0.939


#ifdef __MAKECINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#endif


#ifndef __TCalcPhysVars_H__
#define __TCalcPhysVars_H__

class TCalcPhysVars{
    
private:
    
    // globals
    TString     InFileName;
    TTree       * InTree    ,   * OutTree;
    int         Nentries    ,   A   ,   counter ,   Entry;
    bool        PrintedDataThisEvent    ,   RawData;         // input types....
    TCalculations * calc;
    TPlots      * plot;
    TStopwatch  timer;
    TRandom2    * random;
    
    
  
    
    
    // Integer variables
    int             targ_type;                        // atomic number
    int             N_t     , Np   , N_PiP , N_PiN;             // number of particles
    int             AllRecoilsAbove350;
    int             N_n     , N_id[20]  , NumberOfNeutrons  , NumberOfNeutronsAndNeutronsSC;
    Int_t           unsortedpPID[20]    , unsortedpCut[20]  , pPID[20]  , pCut[20];
    Int_t           NpBack;
    
    
    // TLorentzVector variables
    TLorentzVector  q       , TargetAtRest , p4momentum , A_4Vector , A_Np_1_4Vector , Wtilde;
    std::vector<TLorentzVector>     P   , pBack;      // all protons, and backward going protons....
    
    
    // TVector3 variables
    TVector3        q3Vector, UnsortedPp[20];
    TVector3        *P_e    , eVertex   ;
    TVector3        pLead   , Pcm       , Pmiss     ;
    std::vector<TVector3>     Pp        , pVertex   ; // momenta in q system (q.Vect()||z)
    
    
    
    
    // Float variables
    Float_t         Px_e        , Py_e      , Pz_e      , X_e     ,   Y_e     ,   Z_e   ;
    Float_t         N_Px[20]    , N_Py[20]  , N_Pz[20]   ;                          //for raw data the electrons are in [0]
    Float_t         N_X[20]     , N_Y[20]   , N_Z[20]   ;                           //for raw data the electrons are in [0]
    Float_t         Xb          , Q2        , W         , Nu;
    Float_t         PpX[20]     , PpY[20]   , PpZ[20]   , pMag[20];                 //proton momentum - up to 20 protons
    Float_t         unsortedCTOF[20]        , CTOF[20];                             //corrected TOF and after proton sorting
    Float_t         unsortedpEdep[20]       , pEdep[20];                           //proton energy deposition in TOF scintillators
    Float_t         unsortedX_proton[20]    , unsortedY_proton[20] , unsortedZ_proton[20]; //reaction vertex position
    Float_t         pZ[20]      , Max_Z_vertex  ;
    Float_t         mA          , CoulombDeltaE;                                     // atomic mass
    Float_t         mA_3        , CoulombDeltaE_3;                                  // for 3NSRC - the A-3 mass
    Float_t         Emiss       , ppHe3cmEmiss, Mmiss[3]  , ppHe3cmMmiss , ppHe3Mmiss;   // Mmiss of 1, 2 and 3 protons, 3He-hit Mmiss
    Float_t         PoverQ      , ThetaPQ;
    Float_t         Total_Pt    , alpha_q       , alpha[20]     , SumAlpha;
    Float_t         N_EC_X[20]  , N_EC_Y[20]    , N_EC_Z[20]    , N_TOF[20] ;
    Float_t         N_PathSC[20], N_TimeSC[20]  , ePathSC       , eTimeSC;          // negative particle - we use the [0] - leading electron
    Float_t         pPathSC[20] , pTimeSC[20] ;                                     // positive particles - we are interested in protons...
    Float_t         q_phi       , q_theta       , Pmiss_Phi;
    Float_t         W2tilde     , Xbtilde;
    Float_t         A2mA        , theta_pq[20]  , pBackAlpha[20], SumpBackAlpha;

    
    
    
    
    
    
public:
    
    
    ~TCalcPhysVars             (){   }; 
    TCalcPhysVars              ( TTree * fInTree, TTree * fOutTree, int fA = 12, bool fRawData = true );
    void DifferentTargets       ();
    void InitializeInputTree    ();
    void InitializeOutputTree   ();
    void AcquireAndSortNucleons ( int , bool DoPrint = false);
    void ComputePhysVariables   ();
    void RotateVector           ( TVector3* );
    void PrintData              ( TString Stage );

    TTree * GetOutTree          (){return OutTree;};

    void CloseCalculation       ();

};




#endif




