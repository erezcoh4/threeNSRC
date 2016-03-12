//  TSchemeEG2.cpp
//  Created by Erez Cohen on 08/10/15.
//  Scheme EG2 data tree 

#include "TSchemeEG2.h"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  pppSRC : one with Np > 0 ,Xb > 1.2 and 0.3 < P(miss) < 1.0 GeV/c
TSchemeEG2::TSchemeEG2( TString file_name , double NeventsFraction, bool fDoRawData, TString type ){
    Printf("Scheming %s for %s...",file_name.Data(),type.Data());
    DoRawData = fDoRawData;
    LoadInTree(file_name);
    CreateOutTree(file_name,type);
    std::cout << "\n";
    for (Long64_t i = 0; i < Nentries*NeventsFraction ; i++) {
        if (i%100==0)
            printf("\rScheming...[%.0f%%]",(100.0*(float)i)/Nentries);
        InTree -> GetEntry(i);
        if(DoRawData){
            Px_e = N_Px[0]; Py_e = N_Py[0]; Pz_e = N_Pz[0];
        }
        TVector3 * q        = new TVector3( - Px_e , - Py_e , 5.009 - Pz_e );
        TVector3 * StruckPp = new TVector3();
        for (int p = 0 ; p < P_nmb ; p++){
            if( P_cut[p] == 1 && P_PID[p] == 1 ){    // this is a proton with momentum |p|<2.8 and 'good' CTOF
                TVector3 * Pp = new TVector3(PpX[p],PpY[p],PpZ[p]);
                if (Pp->Mag() > StruckPp->Mag())    // this is a faster proton
                    StruckPp = Pp;
            }
        }
        TVector3 Pmiss = *StruckPp - *q;
        if( (targ_type==2) && (0.3 < Pmiss.Mag()) && (Pmiss.Mag() < 1.0) )// && (Xb > 1.)
            OutTree -> Fill();
    }
    WriteOutFile();
}











//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  pppSRC : one with Np > 0 ,Xb > 1.2 and 0.3 < P(miss) < 1.0 GeV/c
void TSchemeEG2::SRCPmissXb( TString file_name , Float_t NeventsFraction ){
    LoadInTree(file_name);
    CreateOutTree(file_name,"SRCPmissXb");
    std::cout << "\n";
    for (Long64_t i = 0; i < Nentries*NeventsFraction ; i++) {
        if (i%(Nentries/10000)==0) printf("\rscheming...[%.1lf%%] (wrote %d events)    ",(100.0*(float)i)/Nentries, (int)OutTree->GetEntries());
        InTree -> GetEntry(i);
        
        if(DoRawData){
            Px_e = N_Px[0];
            Py_e = N_Py[0];
            Pz_e = N_Pz[0];
        }
        
        TVector3 * q        = new TVector3( - Px_e , - Py_e , 5.009 - Pz_e );
        TVector3 * StruckPp = new TVector3();
        for (int p = 0 ; p < P_nmb ; p++){
            if( P_cut[p] == 1 && P_PID[p] == 1 ){    // this is a proton with momentum |p|<2.8 and 'good' CTOF
                TVector3 * Pp = new TVector3(PpX[p],PpY[p],PpZ[p]);
                if (Pp->Mag() > StruckPp->Mag())    // this is a faster proton
                    StruckPp = Pp;
            }
        }
        TVector3 Pmiss = *StruckPp - *q;
        if( (P_nmb>0) && (targ_type==2) && (0.3 < Pmiss.Mag()) && (Pmiss.Mag() < 1.0) && (Xb > 1.) )
            OutTree -> Fill();
    }
    WriteOutFile();
}










//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// scheme for proton backwards (xB > 0.5 and all positive are protons....)
void TSchemeEG2::pBack( TString file_name , Float_t NeventsFraction ){
    
    float BeamEnergy            = 5.009;// GeV
    float MinimalTheta_pBack    = 110;  // deg.
    float MinimalProtonMomemtum = 0.3;  // GeV/c
    float MinimalXb             = 0.0;
    int counter                 = 0;
    double rad2deg              = TMath::RadToDeg();
    Float_t             theta_pq[20];
    TVector3                * Pp[20];
   
    LoadInTree(file_name);
    CreateOutTree(file_name,"pBack");
    std::cout << "\n";
    for (Long64_t i = 0; i < Nentries*NeventsFraction ; i++) {
        int     NpBack          = 0;
        int     NGoodBackProtons= 0;
        if (i%(Nentries/5000)==0) printf("\rscheming...[%.1lf%%]  (wrote %d events)        "
                                         ,(100.0*(float)i)/Nentries,counter);
        InTree -> GetEntry(i);
        if ( targ_type==2  &&  Xb>MinimalXb ){
            
            if(DoRawData){
                Px_e = N_Px[0];
                Py_e = N_Py[0];
                Pz_e = N_Pz[0];
            }
            
            TVector3 q3Vector( - Px_e , - Py_e , BeamEnergy - Pz_e );
            
            for (int p = 0 ; p < P_nmb ; p++){
                
                Pp[p] = new TVector3(PpX[p],PpY[p],PpZ[p]);
                
                theta_pq[p] = rad2deg*( Pp[p] -> Angle( q3Vector ) );
                
                if ( theta_pq[p] > MinimalTheta_pBack ){
                    NpBack ++ ;
                    if ( Pp[p]->Mag() > MinimalProtonMomemtum )
                        NGoodBackProtons ++ ;
                }
            }
            if ( NGoodBackProtons==NpBack && NpBack>0 ){
                counter++;
                OutTree -> Fill();
            }
        }
    }
    WriteOutFile();
}









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TSchemeEG2::LoadInTree(TString file_name){
    InFileName  = file_name;
    InFile      = new TFile(Form("/Users/erezcohen/Desktop/DataMining/3NSRC/EG2_DATA/%s.root",file_name.Data())); // my mac
//    InFile      = new TFile(Form("/home/erez/DataMining/3NSRC/EG2_DATA/%s.root",file_name.Data()));
    InTree      = (TTree*)InFile -> Get("T");
    Nentries    = InTree -> GetEntries();
    InTree -> SetBranchAddress("P_nmb"              , &P_nmb);  // number of positive particles
    InTree -> SetBranchAddress("Xb"                 , &Xb);
    InTree -> SetBranchAddress("targ_type"          , &targ_type);
    
    if(DoRawData) { // For the raw data we have the following branches
        
        // positive particles momenta
        InTree -> SetBranchAddress("P_Px"           , &PpX);    // positive particles momenta
        InTree -> SetBranchAddress("P_Py"           , &PpY);
        InTree -> SetBranchAddress("P_Pz"           , &PpZ);
        InTree -> SetBranchAddress("P_PID"          , &P_PID);    // positive particles momenta
        InTree -> SetBranchAddress("P_cut"          , &P_cut);    // positive particles momenta
        // negative particles momenta
        InTree -> SetBranchAddress("N_Px"           , &N_Px);    // negative particles momenta (electron is the first)
        InTree -> SetBranchAddress("N_Py"           , &N_Py);    // negative particles momenta (electron is the first)
        InTree -> SetBranchAddress("N_Pz"           , &N_Pz);    // negative particles momenta (electron is the first)
    
    } else {
        
        InTree -> SetBranchAddress("Px"             , &PpX);
        InTree -> SetBranchAddress("Py"             , &PpY);
        InTree -> SetBranchAddress("Pz"             , &PpZ);
        InTree -> SetBranchAddress("Px_e"           , &Px_e);
        InTree -> SetBranchAddress("Py_e"           , &Py_e);
        InTree -> SetBranchAddress("Pz_e"           , &Pz_e);
        
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TSchemeEG2::CreateOutTree(TString file_name, TString suffix){
    //Create a new file + a clone of old tree in new file
    OutFile = new TFile(Form("/Users/erezcohen/Desktop/DataMining/3NSRC/SchemedData/%s_%s.root",file_name.Data(),suffix.Data()),"recreate");
//    OutFile = new TFile(Form("/home/erez/DataMining/3NSRC/SchemedData/%s_%s.root",file_name.Data(),suffix.Data()),"recreate");
    OutTree = InTree -> CloneTree(0);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TSchemeEG2::WriteOutFile(){
    Printf("\nOutTree has %d entries",(int)OutTree->GetEntries());
    OutTree -> AutoSave();
    OutFile -> Write();
    OutFile -> Close();
    delete InFile;
    delete OutFile;
}

