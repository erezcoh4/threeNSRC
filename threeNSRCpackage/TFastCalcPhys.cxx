//  TFastCalcPhys.C
//  Created by Erez Cohen 01/11/15.
#include "TFastCalcPhys.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TFastCalcPhys::TFastCalcPhys( TTree * fInTree , TTree * fOutTree , int fA ){
    InTree  = fInTree;
    OutTree = fOutTree;
    A       = fA;
    Ebeam   = 5.009;
    calc    = new TCalculations();
    calc -> TargetMassAndDeltaE( A      , &mA   , &CoulombDeltaE);
    InitializeInputTree();
    InitializeOutputTree();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TFastCalcPhys::ComputePhysVariables(int entry , bool DoPrint){
    // everything done in the lab frame for rapid processing
    
    InTree -> GetEntry(entry);

    // electron
    e.SetVectM( TVector3(0    ,   0   ,  Ebeam) , Me);
    ePrime.SetVectM( TVector3(Px_e ,  Py_e ,  Pz_e) , Me );
    q = e - ePrime;
    
    
    // get protons - energy loss correction and Coulomb corrections
    if (!Uprotons.empty())    Uprotons.clear();
    if (!protons.empty())     protons.clear();
    
    for (int p = 0 ; p < Np ; p++ ){
        Proton.SetVectM(TVector3(PpX[p],PpY[p],PpZ[p]) , Mp);
        Proton.SetVect(calculations -> EnergyLossCorrrection(Proton.Vect()));
        Proton.SetVect(calculations -> CoulombCorrection( Proton.Vect() , CoulombDeltaE ));
        Uprotons.push_back( Proton );
        pMag[p] = Proton.P();
    }
    
    // sort the protons....
    
    int InDeX[20];
    TMath::Sort(Np , pMag , InDeX);      // sort the protons according to their magnitude
    for (int p = 0 ; p < Np ; p++ )
        protons.push_back(Uprotons.at(InDeX[p]));

    Plead       = protons.at(0).Vect();
    Pmiss       = Plead - q.Vect();
    Prec        = (Np == 2) ? protons.at(1).Vect() : ( (Np == 3) ? protons.at(0).Vect() + protons.at(1).Vect() : TVector3() ) ;
    Pcm         = (Np == 1 || Np == 2 || Np == 3) ? Pmiss + Prec : TVector3()  ;
    //     p/q , 𝜃(p,q)
    PoverQ      = Plead.Mag() / q.P();
    ThetaPQ     = RAD2DEG*( Plead.Angle(q.Vect()) );
    OutTree -> Fill();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TFastCalcPhys::PrintData(TString Stage){
    timer.Stop();
    std::cout << "--------------------------" << std::endl ;
    timer.Start();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TFastCalcPhys::InitializeInputTree(){
    InTree -> SetBranchAddress("targ_type"      , &targ_type);
    InTree -> SetBranchAddress("Xb"             , &Xb);
    InTree -> SetBranchAddress("P_nmb"          , &Np);
    InTree -> SetBranchAddress("T_nmb"          , &Ntot);
    InTree -> SetBranchAddress("Q2"             , &Q2);
    InTree -> SetBranchAddress("Px_e"           , &Px_e);
    InTree -> SetBranchAddress("Py_e"           , &Py_e);
    InTree -> SetBranchAddress("Pz_e"           , &Pz_e);
    InTree -> SetBranchAddress("Px"             , &PpX);
    InTree -> SetBranchAddress("Py"             , &PpY);
    InTree -> SetBranchAddress("Pz"             , &PpZ);
    Nentries    = InTree -> GetEntries();
    std::cout << "Initialized Input InTree TFastCalcPhys, Nentries = " <<  Nentries << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TFastCalcPhys::InitializeOutputTree(){
    OutTree -> Branch("Np"                  ,&Np                    , "Np/I");
    OutTree -> Branch("Ntot"                ,&Ntot                  , "Ntot/I");

    OutTree -> Branch("targ_type"           ,&targ_type             , "targ_type/I");
    OutTree -> Branch("Xb"                  ,&Xb                    , "Xb/F");
    OutTree -> Branch("Q2"                  ,&Q2                    , "Q2/F");
    OutTree -> Branch("ThetaPQ"             ,&ThetaPQ               , "ThetaPQ/F");
    OutTree -> Branch("PoverQ"              ,&PoverQ                , "PoverQ/F");

    OutTree -> Branch("Plead"               ,"TVector3"             ,&Plead);
    OutTree -> Branch("Pmiss"               ,"TVector3"             ,&Pmiss);
    OutTree -> Branch("Prec"                ,"TVector3"             ,&Prec);
    OutTree -> Branch("Pcm"                 ,"TVector3"             ,&Pcm);

    OutTree -> Branch("e"                   ,"TLorentzVector"       ,&e);
    OutTree -> Branch("ePrime"              ,"TLorentzVector"       ,&ePrime);
    OutTree -> Branch("q"                   ,"TLorentzVector"       ,&q);

    OutTree -> Branch("protons"             ,&protons );    // final state

    std::cout << "Initialized Output Tree TFastCalcPhys" << std::endl;
}







