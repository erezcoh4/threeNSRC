//  TSchemeEG2.h
//  Created by Erez Cohen on 8/10/15.
#include "TAnalysisEG2.h"



using namespace std;

#ifndef __TSchemeEG2_H__
#define __TSchemeEG2_H__

class TSchemeEG2{
private:
    
    // globals
    TString InFileName;
    TCut    InnerCut;
    TFile * InFile;
    TTree * InTree;
    Long64_t Nentries;
    TFile * OutFile;
    TTree * OutTree;
    TCalculations calc;
    bool DoRawData;
    
    // integer variables
    Int_t P_nmb     , P_cut[20], P_PID[20]  ;   //positive particles
    Int_t N_SCStat  , Status;
    Int_t N_n       , N_id[20] , NumberOfNeutrons , NumberOfNeutronsAndNeutronsSC;
    Int_t targ_type;
    
    // TVector3 variables
    TVector3 NeutronECHit[20]   , NeutronMomentum[20] , PnStruck;

    // Float variables
    Float_t Xb          , Q2;
    Float_t PpX[20]     , PpY[20]   , PpZ[20];  //proton momentum - up to 20 protons
    Float_t Px_e        , Py_e      , Pz_e  ,   X_e , Y_e , Z_e , STT;
    Float_t N_Px[20]    , N_Py[20]  , N_Pz[20]; // for raw data
    Float_t N_Z[20]     ;                       // electron vertex
    Float_t N_EC_X[20]  , N_EC_Y[20], N_EC_Z[20] , N_TOF[20] , N_Beta[20]   ;
    
public:
    
    ~TSchemeEG2(){};
    TSchemeEG2                  (bool fDoRawData=false){ DoRawData = fDoRawData; };
    TSchemeEG2                  ( TString , double , bool fDoRawData=false, TString type = "SRCPmissXb");
    TSchemeEG2                  ( Float_t , Float_t , TString , double );
    TSchemeEG2                  ( TString , Float_t , Float_t , double );
    
    
    
    void LoadInTree             (TString );
    void CreateOutTree          (TString, TString);
    void WriteOutFile           ();
    
    void pBack                  ( TString , Float_t);
    void SRCPmissXb             ( TString , Float_t);
    void SetAdressesInputTree   ();
    
};


#endif


