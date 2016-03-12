//  TAnalysisEG2.h
//  Created by Erez Cohen on 1/31/15.

#include "/Users/erezcohen/larlite/UserDev/mySoftware/mySoftwarePackage/TPlots.h"
#include "/Users/erezcohen/larlite/UserDev/mySoftware/mySoftwarePackage/TAnalysis.h"
#include "/Users/erezcohen/larlite/UserDev/mySoftware/mySoftwarePackage/TCalculations.h"


#ifndef __TAnalysisEG2_H__
#define __TAnalysisEG2_H__

class TAnalysisEG2 : public TPlots  , public TAnalysis
{
    
    private:
    
    TString FileName;
    TFile * file;
    TTree * Tree;
    int     Nentries;
    bool    Mixed;
    
    
    
    public:
    
    // TVector3
    std::vector<TVectcor3>   protons ;
    TVector3              Pcm     ,   Pmiss     , * Prec    , * PrecDiff;
    
    // TLorentzVector
    std::vector<TLorentzVector> * p ;
    TLorentzVector              * q ;
    
    // Int_t
    Int_t       A           , Np;
    Int_t       pPID[10]    , pCut[10]  ,   pCTOFCut[10];

    
    // Float_t
    Float_t     alpha[10]   , alpha_q   , SumAlpha  , pEdep[10]     , pCTOF[10];
    Float_t     Xb          , ThetaPQ   ,   PoverQ  , ipAngle       , oopAngle;

    
    
    
    // constructors
    ~TAnalysisEG2               (){};
    TAnalysisEG2                (){};
    TAnalysisEG2                ( TTree *, TString );
    TAnalysisEG2                ( TString FileName , bool mixed = false );
    
 

    
    //--------methods-------
    void SetTrioInitialState    ( bool DoMix = false );
    void Plot_3p_Momemta        ( TCut cut = "" );
    void Plot3pMomenta          ( int );
    void Print_3p_Momemta       ( TCut cut = "" );
    void PrintMomentaSet        ( int , int );
    
    void ShemeTreeToCut         ( TCut );
    void Write2File             ( TTree * , TCut cut = ""  );
    void Mix2File               ( TTree * , TCut cut = ""  );
    void SetBranches            ( bool DoMix = false );
    void SetOutputTreeBranches  ( TTree * );

    
    TH2F * GetCTOFSignalToNoise (TCut, float XMin = -5 , float XMax = 5);

    
    
  
};


#endif



