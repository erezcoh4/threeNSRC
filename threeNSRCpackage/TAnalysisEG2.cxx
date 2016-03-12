//  TAnalysisEG2.cpp
//  Created by Erez Cohen on 1/31/15.

#include "TAnalysisEG2.h"
#define rad2deg TMath::RadToDeg()

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2::TAnalysisEG2( TTree * fInTree, TString fName )
:TPlots( fInTree , fName ){
    Tree        = fInTree;
    Nentries    = Tree -> GetEntries();
    printf("Initiating process for %s with %d Nentries\n",fName.Data(),Nentries);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TAnalysisEG2::TAnalysisEG2( TString FileName , bool mixed )
:TPlots(Form("/Users/erezcohen/Desktop/DataMining/AnaFiles/%s",FileName.Data()),(mixed)?"MixedTree":"anaTree",FileName,false){
    FileName = FileName;
    file = new TFile(Form("/Users/erezcohen/Desktop/DataMining/AnaFiles/%s",FileName.Data()));
    Mixed   = mixed;
    Tree    = (Mixed) ? (TTree*)file->Get("MixedTree") : (TTree*)file->Get("anaTree");
    Nentries = Tree -> GetEntries();
    printf("Initiating process for analysis with %d Nentries\n",Nentries);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::Plot_3p_Momemta( TCut cut ){
    SetBranches();
    int N3p = Tree->GetEntries(cut) , n3p = 0;
    Printf("From (%s),\nWill plot %d trios in this cut:",cut.GetTitle(),N3p);
    int Cctr = 0 , EventPlotCtr = 0;
    TCanvas * cTrio[50];
    for ( int entry = 0 ; entry < Tree->GetEntries() ; entry++ ){
        if( Tree -> Draw("Xb",cut,"goff",1,entry) ){
            Tree -> GetEntry( entry );
            SetTrioInitialState();
            n3p ++ ;
            if ( EventPlotCtr%6 == 0 ) {
                Cctr++;
                cTrio[Cctr-1] = CreateCanvas(Form("c%d",Cctr-1),"Divide",3,2);
            }
            EventPlotCtr++ ;
            cTrio[Cctr-1] -> cd ( EventPlotCtr - 6*(Cctr-1) );
            Plot3pMomenta( entry );
            cTrio[Cctr-1] -> cd ( EventPlotCtr - 6*(Cctr-1) + 1);
            PrintMomentaSet( entry , n3p );
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::Plot3pMomenta(int entry){
    // 3D
    Frame3D( Form("Entry %d, p(c.m.) = %.2f GeV/c",entry,(protons.at(0)+protons.at(1)+protons.at(2)).Mag()) );
    Plot3DAxesSystem();
    Line3D( protons.at(0)   , 2     , 4);
    Line3D( protons.at(1)   , 30    , 4);
    Line3D( protons.at(2)   , 4     , 4);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::Print_3p_Momemta( TCut cut ){
    ShemeTreeToCut(cut);
    SetBranches();
    for ( int entry = 0 ; entry < Tree->GetEntries() ; entry++ ){
        Tree -> GetEntry( entry );
        SetTrioInitialState();
        PrintMomentaSet( entry , Tree->GetEntries() );
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::PrintMomentaSet(int entry, int n){
    Printf("Event %d (number %d)",entry,n);
    for ( size_t p  = 0; p < protons.size() ; p++ ){
        Printf("initial p(%lu) = %f , %f , %f [α=%.1f]"
               , p+1 , protons.at(p).x(), protons.at(p).y(), protons.at(p).z(),alpha[p]);
    }
    PrintOut3Vector( Pcm , "p(c.m.)" );
    if (protons.size() == 3){
        SHOW(ipAngle);
        SHOW(oopAngle);
    }
    PrintLine();
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::SetTrioInitialState( bool DoMix ){
    if (!protons.empty())        protons.erase(protons.begin(),protons.end());
    protons.push_back( p->at(0).Vect() - q->Vect() );
    Pmiss       = protons.at(0);
    Pcm         = protons.at(0);
    Prec        = new TVector3();
    if (p->size() > 1) {
       for (size_t i = 1 ; i < p->size()  ; i++){
            *Prec       += p->at(i).Vect();
            Pcm         += p->at(i).Vect();
            protons     . push_back( p->at(i).Vect() );
        }
       if (p->size() == 3) {
            // In and Out of plane angle that p(miss) forms with p(rec), where p(rec)=p(2)+p(3) is in the plane formed by p(2) & p(3)
            // n  = protons.at(1).Cross(protons.at(2)).Unit();    // perpendicular to the plane formed by p(2) & p(3)
            oopAngle    = 90 - TMath::RadToDeg()*( protons.at(0).Angle(protons.at(1).Cross(protons.at(2)).Unit()) );
            // Pmiss Projection On p2-p3 plane = Pmiss - (Pmiss·n) n;
            ipAngle     = TMath::RadToDeg()*( protons.at(0) - (protons.at(0).Dot(protons.at(1).Cross(protons.at(2)).Unit()))*protons.at(1).Cross(protons.at(2)).Unit()).Angle(*Prec);
        }
    }
    if (!DoMix) {
        for (size_t i = 0 ; i < p->size()  ; i++)
            pCTOFCut[i] = pPID[i] * pCut[i] ; // pPID - polynomial proton id based on CTOF. pCut: momentum < 2.8 GeV/c
    }
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::ShemeTreeToCut( TCut cut ){
    TTree * tmpTree = Tree -> CloneTree(0);
    int N2File = Tree->GetEntries(cut);
    for ( int i = 0 ; i < Tree -> GetEntries() ; i++ ){
        if ( Tree -> Draw("Xb" , cut , "goff" , 1 , i) ){
            Tree -> GetEntry( i );
            tmpTree -> Fill();
            if (tmpTree->GetEntries()%(N2File/20)==0)
                Printf("[%0.f%%]",100.*(float)(tmpTree->GetEntries())/N2File);
        }
    }
    Tree = tmpTree;
    Printf("%lld entries passed the cut",Tree -> GetEntries());
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::Write2File( TTree * OutTree , TCut cut ){
    SetBranches();
    SetOutputTreeBranches(OutTree);
    for ( int entry = 0 ; entry < Tree -> GetEntries() ; entry++ ){
        Tree -> GetEntry( entry );
        SetTrioInitialState();
        OutTree -> Fill();
        if (OutTree->GetEntries()%10==0) {
            Printf("[%0.f%%]",100.*(float)(OutTree->GetEntries())/Tree->GetEntries());
            PrintMomentaSet( entry , OutTree->GetEntries() );
        }
    }
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::Mix2File(TTree * OutTree, TCut cut){ // mix all trios, all done in q-Pmiss system
    SetBranches( true );
    SetOutputTreeBranches(OutTree);
    TLorentzVector  mixed_q         , mixed_p1      , mixed_p2      , mixed_p3;
    Float_t         mixed_CTOFCut[3], mixed_pEdep[3],mixed_alpha[10] , mixed_alpha_q , mixed_SumAlpha;
    Printf("will mix %lld entries (%lldx%lldx%lld)"
           ,Tree->GetEntries()*(Tree->GetEntries()-1)*(Tree->GetEntries()-2),Tree -> GetEntries(),Tree -> GetEntries()-1,Tree -> GetEntries()-2);
    
    for ( int i = 0 ; i < Tree -> GetEntries() ; i++ ){
        Tree -> GetEntry( i );
        mixed_p1        = p->at(0);
        mixed_CTOFCut[0]= pCTOFCut[0];
        mixed_pEdep[0]  = pEdep[0];
        mixed_alpha[0]  = alpha[0];
        mixed_alpha_q   = alpha_q;
        mixed_SumAlpha  = mixed_alpha[0] - mixed_alpha_q;
        
        
        for ( int j = 0 ; j < Tree -> GetEntries() ; j++ ){
            if (j!=i) {
                Tree -> GetEntry( j );
                mixed_p2        = p->at(1);
                mixed_CTOFCut[1]= pCTOFCut[1];
                mixed_pEdep[1]  = pEdep[1];
                mixed_alpha[1]  = alpha[1];
                mixed_SumAlpha += mixed_alpha[1];
                
                for ( int k = 0 ; k < Tree -> GetEntries() ; k++ ){
                    if (k!=i && k!=j) {
                        Tree -> GetEntry( k );
                        mixed_p3        = p->at(1);
                        mixed_CTOFCut[2]= pCTOFCut[2];
                        mixed_pEdep[2]  = pEdep[2];
                        mixed_alpha[2]  = alpha[2];
                        mixed_SumAlpha += mixed_alpha[2];
                       
                        
                        Tree -> GetEntry( i ); // retain all electron variables from original event....
                        if (!p->empty())     p->erase(p->begin(),p->end());
                        p -> push_back(mixed_p1);
                        p -> push_back(mixed_p2);
                        p -> push_back(mixed_p3);
                        for (int l = 0 ;  l < Tree -> GetEntries() ; l++ ) {
                            alpha[l]    = mixed_alpha[l];
                            pCTOFCut[l] = mixed_CTOFCut[l];
                            pEdep[l]    = mixed_pEdep[l];
                        }
                        SumAlpha    = mixed_SumAlpha;

                        SetTrioInitialState(true);
                        if (OutTree->GetEntries()%50==0)
                            PrintMomentaSet( OutTree->GetEntries() , OutTree->GetEntries() );
                        
                        OutTree -> Fill();
                    }
                }
            }
        }
    }
    Printf("Filled %lld Mixed entries",OutTree->GetEntries());
}








//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::SetBranches(bool DoMix){ // Dec-21,2015
    p = 0;    q = 0;    Prec = 0;   // Pcm = 0;    Pmiss = 0;      
    Tree -> SetBranchAddress("Np"       , &Np);
    Tree -> SetBranchAddress("P"        , &p);
    Tree -> SetBranchAddress("q"        , &q);
    Tree -> SetBranchAddress("alpha"    , &alpha);
    Tree -> SetBranchAddress("alpha_q"  , &alpha_q);
    Tree -> SetBranchAddress("Xb"       , &Xb);
    Tree -> SetBranchAddress("ThetaPQ"  , &ThetaPQ);
    Tree -> SetBranchAddress("PoverQ"   , &PoverQ);
    Tree -> SetBranchAddress("SumAlpha" , &SumAlpha);
    Tree -> SetBranchAddress("pEdep"    , &pEdep);
    Tree -> SetBranchAddress("pCTOF"    , &pCTOF);
    if (DoMix) {
        Tree -> SetBranchAddress("pCTOFCut"     , &pCTOFCut);
    } else {
        Tree -> SetBranchAddress("pPID"     , &pPID);
        Tree -> SetBranchAddress("pCut"     , &pCut);
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TAnalysisEG2::SetOutputTreeBranches( TTree * OutTree ){
    OutTree -> Branch("A"       ,&A                 ,"A/I");                    // mass number
    OutTree -> Branch("Xb"      ,&Xb                ,"Xb/F");
    OutTree -> Branch("SumAlpha",&SumAlpha          ,"SumAlpha/F");
    OutTree -> Branch("ThetaPQ" ,&ThetaPQ           ,"ThetaPQ/F");
    OutTree -> Branch("PoverQ"  ,&PoverQ            ,"PoverQ/F");
    OutTree -> Branch("Np"      ,&Np                ,"Np/I");
    OutTree -> Branch("P"                           ,&p);                   // final 4-momenta - std::vector<TLorentzVector>
    OutTree -> Branch("protons"                     ,&protons);             // initial 3-momenta - std::vector<TVector3>
    OutTree -> Branch("q"       ,"TLorentzVector"   ,&q);
    OutTree -> Branch("Pcm"     ,"TVector3"         ,&Pcm);
    OutTree -> Branch("Pmiss"   ,"TVector3"         ,&Pmiss);
    OutTree -> Branch("Prec"    ,"TVector3"         ,&Prec);
    OutTree -> Branch("oopAngle",&oopAngle          ,"oopAngle/F");         // out-of-plane Angle of p(miss) outside the plane formed by p2 & p3
    OutTree -> Branch("ipAngle" ,&ipAngle           ,"ipAngle/F");         // in-plane Angle of p(miss) with p(rec)
    OutTree -> Branch("alpha"   ,&alpha             ,"alpha[Np]/F");
    OutTree -> Branch("alpha_q" ,&alpha_q           ,"alpha_q/F");
    OutTree -> Branch("pEdep"   ,&pEdep             ,"pEdep[Np]/F");
    OutTree -> Branch("pCTOF"   ,&pCTOF             ,"pCTOF[Np]/F");
    OutTree -> Branch("pCTOFCut",&pCTOFCut          ,"pCTOFCut[Np]/I");
}































//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TH2F * TAnalysisEG2::GetCTOFSignalToNoise(TCut cut, float XMin, float XMax){ // Dec-30, 2015
    // (a) project the 2D histogram into each of the axes (p2 and p3),
    // (b) calculate signal to noise ratio of each histogram under the peak
    const Int_t Nbins   = 100;
    int Niterations     = 5;
    TH2F * h2 = H2WithProjections( "pCTOF[1]","pCTOF[2]", cut , Nbins,-20,20,Nbins,-20,20
                                  ,Form("positive p 2 vs. positive p 3")
                                  ,"CTOF p 2 [ns]","CTOF p 3 [ns]");
    TH1F * hP[2] , * hBackground[2] , * hSignal[2];
    TSpectrum * spectrum = new TSpectrum();
    TCanvas * c = CreateCanvas("projections","DivideSquare");
    for (int i = 0 ; i < 2 ; i++ ) {
        c -> cd(i+1);
        hP[i] = (i==0) ? (TH1F*) h2 -> ProjectionX ("hPx") : (TH1F*) h2 -> ProjectionY ("hPy");
        SetFrame(hP[i],Form("positive particle %d",i+2),"CTOF [ns]","");
        hP[i] -> Draw();
        hBackground[i] = (TH1F*) spectrum -> Background (hP[i]);
        SetFrame(hBackground[i],"","","",4,46);
        hBackground[i] -> Draw("same");
        c -> cd(i+3);
        hSignal[i] = new TH1F(Form("hSignal%d",i),"",Nbins,-20,20);
        hSignal[i] -> Add(hP[i],1);
        hSignal[i] -> Add(hBackground[i],-1);
        SetFrame(hSignal[i],Form("background subtracted p %d",i+2),"CTOF [ns]","",38,38,3001);
        hSignal[i] -> SetMinimum(0);
        hSignal[i] -> Draw();
        Line(XMin , 0 , XMin , hSignal[i]->GetMaximum() , 2 , 2 );
        Line(XMax , 0 , XMax , hSignal[i]->GetMaximum() , 2 , 2 );
        Text(-10,0.5*hSignal[i]->GetMaximum()
             ,Form("removed %.1f%% background" ,100 * (1-IntegralH1D(hSignal[i],XMin,XMax) / IntegralH1D(hP[i],XMin,XMax))));
    }
    
    // 2d spectra analysis
    TCanvas * c2 = CreateCanvas("background subtraction 2D","Divide",2,1);
    c2 -> cd(1);
    h2 -> Draw("colz");
    TH2F* hBack2 = new TH2F("hBack","",Nbins,-20,20,Nbins,-20,20);
    TSpectrum2 * s2 = new TSpectrum2();
    float ** source = new float *[Nbins];
    for (int i = 0 ; i < Nbins ; i++) {
        source[i] = new float[Nbins];
        for (int j = 0 ; j < Nbins ; j++) {
            source[i][j] = h2 -> GetBinContent(i+1,j+1);
        }
    }
    s2 -> Background(source,Nbins,Nbins,Niterations,Niterations,TSpectrum2::kBackDecreasingWindow,TSpectrum2::kBackSuccessiveFiltering);
    for (int i = 0; i < Nbins; i++)
        for (int j = 0; j < Nbins; j++)
            hBack2 -> SetBinContent(i + 1,j + 1, source[i][j]);
    TH2F* hSig2D = new TH2F("hSignal2","",Nbins,-20,20,Nbins,-20,20);
    hSig2D -> Add(h2,1);
    hSig2D -> Add(hBack2,-1);
    c2 -> cd(2);
    hSig2D -> Draw("colz");
    Box(XMin , XMin , XMax , XMax , 2 , 0.1 );
    Text(-10,10 ,Form("removed %.1f%% background" ,100 * (1-IntegralH2D(hSig2D,XMin,XMin,XMax,XMax) / IntegralH2D(h2,XMin,XMin,XMax,XMax))));
    return h2;
}


















