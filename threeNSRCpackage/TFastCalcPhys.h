//  TFastCalcPhys.h
//  Created by Erez Cohen on 8/10/15.
#include "TAnalysisEG2.h"

#define MyHomePath "/Users/erezcohen/Desktop"
#define MyDMPath Form("%s/DataMining",MyHomePath)
#define MyDMAnaPath Form("%s/AnaFiles",MyDMPath)

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

#ifndef __TFastCalcPhys_H__
#define __TFastCalcPhys_H__

class TFastCalcPhys{
    
private:
    
    // globals
    TString     InFileName;
    TTree       * InTree    ,   * OutTree;
    int         Nentries    ,   A   ;
    TCalculations * calculations;
    TStopwatch  timer;
    
    
    
    int             targ_type   , Ntot          , Np        ;
    Float_t         Px_e        , Py_e          , Pz_e      , ThetaPQ   ,   PoverQ  ,   Q2  ,   Xb , mA , Ebeam ,   CoulombDeltaE;
    Float_t         PpX[20]     , PpY[20]       , PpZ[20]   , pMag[20];
    
    
    TVector3                        Plead       , Pcm       , Pmiss     ,   Prec;
    TLorentzVector                  e           , ePrime    , q         ,   Proton;
    std::vector<TLorentzVector>     Uprotons    , protons   ;
    
    
    
    
public:
    
    
    ~TFastCalcPhys              (){   };
    TFastCalcPhys               ( TTree * fInTree, TTree * fOutTree, int fA = 12 );
    void InitializeInputTree    ();
    void InitializeOutputTree   ();
    void ComputePhysVariables   ( int , bool DoPrint = false);
    void PrintData              ( TString Stage );

    TTree * GetOutTree          (){return OutTree;};

};




#endif




