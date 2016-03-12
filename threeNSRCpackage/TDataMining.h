//  TDataMining.h
//  Created by Erez Cohen on 5/15/15.


#include "TAnalysisEG2.h"
#include "TSchemeEG2.h"
#include "TCalcPhysVars.h"

#ifndef __TDataMining_H__
#define __TDataMining_H__


class TDataMining{
    
    
private:
    
    TCalculations calc;
    
public:
    
    // constructor
    TDataMining                                 (){printf("Constructing TDataMining\n");};
    ~TDataMining                                (){};
    
    TAnalysisEG2 * CombineRawDataChain          (int, int Nthreads = 50);
    TAnalysisEG2 * GoThroughSchemingChain       (TString,int,TString Type = "RawData");
    TAnalysisEG2 * GoThroughSchemingChain       (int);
    TAnalysisEG2 * GoThroughSchemingNStruck     (int);
    TAnalysisEG2 * GoThroughSchemingMyNStruck   (int);
    
    TAnalysisEG2 * GetTriosRawData              (int);
    TAnalysisEG2 * GetTriosData                 (int);
    TAnalysisEG2 * GetMixTriosData              (int);

    
    TAnalysisEG2 * GetpBackAnalysisReady        (int,float evntsfrac = 1.);
   
    
    TCutG *        ProtonEdepVsMomentumCut      (int);
    TCutG *        ProtonEdepVsMomentumCut      (TString);


};


#endif
