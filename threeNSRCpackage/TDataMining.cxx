//  TDataMining.cpp
//  Created by Erez Cohen on 8/10/15.
#include "TDataMining.h"





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::CombineRawDataChain(int A, int Nthreads){ // Jan-14, 2016
    TString Name        = "anaTree";
    TString FileName    = "Ana_RawData_C12_SRCPmiss";
    TString Path        = "/Volumes/Lacie/DataMining/AnaFiles/MTAnaFiles";
    TChain * chain      = new TChain(Name);
    for (int i_thread = 0 ; i_thread < Nthreads ; i_thread ++ )
        chain -> Add(Form("%s/%s_%d.root",Path.Data(),FileName.Data(),i_thread));
    Printf("Added %d trees to a total of %lld events",Nthreads,chain->GetEntries());
    TAnalysisEG2 * ana = new TAnalysisEG2( (TTree*)chain , "Combined RawData trees" );
    return ana;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::GoThroughSchemingChain(int A){
    TString target = calc.TargetAsString(A);
//    TFile SchFile( Form("/Users/erezcohen/Desktop/DataMining/3NSRC/SchemedData/DATA_%s_SRCPmissXb.root",target.Data() ) );
//    if ( SchFile.IsZombie() )
//        TSchemeEG2 * scheme     = new TSchemeEG2( Form("DATA_%s",target.Data() ) , 1 );
//    else
//        printf("File DATA_%s already schemed to SRC kinematics...\n",target.Data());
//    
//    TFile AnaFile( Form("~/Desktop/DataMining/3NSRC/Analysis/Analysis_DATA_%s_SRCPmissXb.root",target.Data() ) );
//    if ( AnaFile.IsZombie() )
//        TCalcPhysVars * calc    = new TCalcPhysVars( Form("DATA_%s_SRCPmissXb.root",target.Data() ), A , 1);
//    else
//        Printf("physical variables already calculated...");
//    
    TAnalysisEG2 * ana          = new TAnalysisEG2( Form("Analysis_DATA_%s_SRCPmissXb.root",target.Data() ) );
    return ana;
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::GoThroughSchemingChain(TString filename, int A, TString Type ){
//    TFile SchFile( Form("/Users/erezcohen/Desktop/DataMining/3NSRC/SchemedData/%s_SRCPmissXb.root",filename.Data() ) );
//    if ( SchFile.IsZombie() ){
//        printf("scheming %s to SRC kinematics...\n",filename.Data());
//        TSchemeEG2 * scheme     = new TSchemeEG2( Form("%s",filename.Data() ) , 1 , true);
//    } else
//        printf("File %s already schemed to SRC kinematics...\n",filename.Data());
//
//    Printf("/Users/erezcohen/Desktop/DataMining/AnaFiles/Ana_%s_SRCPmissXb.root",filename.Data() );
//    TFile AnaFile( Form("/Users/erezcohen/Desktop/DataMining/AnaFiles/Ana_%s_SRCPmissXb.root",filename.Data() ) );
//    if ( AnaFile.IsZombie() )
//        TCalcPhysVars * calc    = new TCalcPhysVars( Form("%s_SRCPmissXb.root",filename.Data() ), A, 1 , Type);
//    else
//        Printf("physical variables already calculated...");
    TAnalysisEG2 * ana          = new TAnalysisEG2( Form("Ana_%s_SRCPmissXb.root",filename.Data() ) );
    return ana;
}









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::GetTriosRawData(int A){ // get the schemed e,e'ppp data tree
    return (new TAnalysisEG2( Form("Raw%s_pppSRCTrios.root",calc.TargetAsString(A).Data())  ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::GetTriosData(int A){ // get the schemed e,e'ppp data tree
    return (new TAnalysisEG2( Form("%s_pppSRCTrios.root",calc.TargetAsString(A).Data())  ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::GetMixTriosData(int A){ // get the schemed e,e'ppp data tree
    return (new TAnalysisEG2( Form("%s_MixedTrios.root",calc.TargetAsString(A).Data()) , true ));
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2 * TDataMining::GetpBackAnalysisReady(int A, float EventsFraction){
    TString target = calc.TargetAsString(A);
//    TFile SchFile( Form("/Users/erezcohen/Desktop/DataMining/3NSRC/SchemedData/RawData_%s_pBack.root",target.Data() ) );
//    if ( SchFile.IsZombie() )
//        TSchemeEG2 * scheme     = new TSchemeEG2( Form("RawData_%s",target.Data() ) , EventsFraction , false , "pBack");
//    else
//        printf("RawData_%s already schemed to p Backwards...\n",target.Data());
//    
//    TFile AnaFile( Form("~/Desktop/DataMining/AnaFiles/Ana_RawData_%s_pBack.root",target.Data() ) );
//    if ( AnaFile.IsZombie() )
//        TCalcPhysVars * calc    = new TCalcPhysVars( Form("RawData_%s_pBack.root",target.Data() ), A , EventsFraction);
//    else
//        Printf("physical variables already calculated...");
//    
    TAnalysisEG2 * ana          = new TAnalysisEG2( Form("Ana_RawData_%s_pBack.root",target.Data() ) );
    return ana;
}







//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TCutG * TDataMining::ProtonEdepVsMomentumCut(TString name){
    if (name == "p1")
        return ProtonEdepVsMomentumCut(0);
    else if (name == "p2")
        return ProtonEdepVsMomentumCut(1);
    else // "p3"
        return ProtonEdepVsMomentumCut(2);
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TCutG * TDataMining::ProtonEdepVsMomentumCut(int p){
    TCutG *cutg = new TCutG(Form("EdepCut_p%d",p),30);
    cutg->SetVarX(Form("P[%d].P()",p));
    cutg->SetVarY(Form("pEdep[%d]",p));
    cutg->SetPoint(0,0.239224,18.4722);
    cutg->SetPoint(1,0.256857,27.4167);
    cutg->SetPoint(2,0.280368,35.9722);
    cutg->SetPoint(3,0.30094,42.8426);
    cutg->SetPoint(4,0.347962,53.213);
    cutg->SetPoint(5,0.368534,59.4352);
    cutg->SetPoint(6,0.392045,65.1389);
    cutg->SetPoint(7,0.474334,61.8981);
    cutg->SetPoint(8,0.489028,52.1759);
    cutg->SetPoint(9,0.533111,43.75);
    cutg->SetPoint(10,0.603644,34.1574);
    cutg->SetPoint(11,0.79761,22.8796);
    cutg->SetPoint(12,1.05623,18.7315);
    cutg->SetPoint(13,1.39714,16.787);
    cutg->SetPoint(14,1.50294,14.3241);
    cutg->SetPoint(15,1.78801,12.5093);
    cutg->SetPoint(16,1.54114,10.4352);
    cutg->SetPoint(16,1.54114,10.4352);
    cutg->SetPoint(17,1.15027,11.7315);
    cutg->SetPoint(18,1.01215,13.5463);
    cutg->SetPoint(19,0.876959,14.0648);
    cutg->SetPoint(20,0.765282,16.1389);
    cutg->SetPoint(21,0.677116,17.8241);
    cutg->SetPoint(22,0.556622,22.8796);
    cutg->SetPoint(23,0.456701,29.6204);
    cutg->SetPoint(24,0.427312,33.7685);
    cutg->SetPoint(25,0.389107,32.8611);
    cutg->SetPoint(26,0.374412,20.287);
    cutg->SetPoint(27,0.32739,16.787);
    cutg->SetPoint(28,0.315635,14.713);
    cutg->SetPoint(29,0.236285,18.7315);
    cutg->SetPoint(30,0.239224,18.4722);
    return cutg;
}