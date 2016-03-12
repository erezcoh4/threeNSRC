//  TCalcPhysVars.C
//  Created by Erez Cohen 01/11/15.
#include "TCalcPhysVars.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TCalcPhysVars::TCalcPhysVars( TTree * fInTree , TTree * fOutTree , int fA, bool fRawData  ){
    InTree  = fInTree;
    OutTree = fOutTree;
    RawData = fRawData;
    A       = fA;
    DifferentTargets();
    calc    = new TCalculations();
    plot    = new TPlots();
    random  = new TRandom2();
    InitializeInputTree();
    InitializeOutputTree();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::DifferentTargets(){   // [GeV/c2]
    calc -> TargetMassAndDeltaE( A      , &mA   , &CoulombDeltaE);
    TargetAtRest.SetVectM( TVector3() , mA  );   // Target initially at rest relative to beam
    A_4Vector   .SetVectM( TVector3() , mA  );
    A2mA        = (float)A/mA;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::AcquireAndSortNucleons(int entry,bool DoPrint){
    InTree -> GetEntry(entry);
    if(RawData) {
        Px_e    = N_Px[0];  Py_e = N_Py[0]; Pz_e    = N_Pz[0];
        X_e     = Y_e  = Z_e     = -100;
        Nu      = 5.009 - sqrt(  Px_e*Px_e + Py_e*Py_e + Pz_e*Pz_e );
    }
    
    // electron
    q3Vector.SetXYZ( - Px_e , - Py_e , 5.009 - Pz_e );
    q_phi   = q3Vector.Phi();
    q_theta = q3Vector.Theta();
    eVertex . SetXYZ( X_e , Y_e , Z_e );


    
    if(RawData)
        Q2  = q3Vector.Mag2() - Nu*Nu;
    
    // get protons - energy loss correction and Coulomb corrections
    if (!Pp.empty())    Pp.clear();
    if (!P.empty())     P.clear();
    if (!pBack.empty()) pBack.clear();
    NpBack  = 0;
    
    for (int p = 0 ; p < Np ; p++ ){
        UnsortedPp[p]   . SetXYZ( PpX[p],PpY[p],PpZ[p] );
        UnsortedPp[p]   = calc -> EnergyLossCorrrection(UnsortedPp[p]);
        UnsortedPp[p]   = calc -> CoulombCorrection( UnsortedPp[p] , CoulombDeltaE );
        pMag[p]         = UnsortedPp[p].Mag();
    }
    
    // sort the protons....
    int InDeX[20];
    TMath::Sort(Np , pMag , InDeX);      // sort the protons according to their magnitude
    
    // If we have 3 protons, randomly choose which is p2 and which is p3
    if( (Np==3) && (random -> Uniform() > 0.5) ){ // switch between p2 and p3 with a probablity of 50%
        int TmpInDeX  = InDeX[1];
        InDeX[1]      = InDeX[2];
        InDeX[2]      = TmpInDeX;
    }
    for (int p = 0 ; p < Np ; p++ ){
        Pp          .push_back(UnsortedPp   [InDeX[p]]    );
        pVertex     .push_back(TVector3(unsortedX_proton[InDeX[p]],unsortedY_proton[InDeX[p]],unsortedZ_proton[InDeX[p]]));
        pZ[p]       = unsortedZ_proton      [InDeX[p]];
        CTOF[p]     = unsortedCTOF          [InDeX[p]];
        pEdep[p]    = unsortedpEdep         [InDeX[p]];
        pCut[p]     = unsortedpCut          [InDeX[p]];
        pPID[p]     = unsortedpPID          [InDeX[p]];
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::ComputePhysVariables(){
    
    pLead       = Pp.at(0);
    Pmiss       = pLead - q3Vector;
    //     p/q , 𝜃(p,q)
    PoverQ      = pLead.Mag() / q3Vector.Mag();
    ThetaPQ     = RAD2DEG*( pLead.Angle(q3Vector) );
    
    //     imidiately transform to q-Pmiss-frame, to calculate everything in this frame
    //     q is the z axis, Pmiss is in x-z plane: Pmiss=(Pmiss[x],0,Pmiss[q])
    Pmiss       .RotateZ(-q_phi);
    Pmiss       .RotateY(-q_theta);
    Pmiss_Phi   = Pmiss.Phi();
    Pmiss       .RotateZ(-Pmiss_Phi);
    q3Vector    .RotateZ(-q_phi);
    q3Vector    .RotateY(-q_theta);
    q           .SetVect( q3Vector );
    q           .SetE( Nu );
    
    
    // A(e,e'p)X missing energy
    double Mre  = mA - Np*Mp;            // taking out the 1 proton off the initial target
    double Tp[20];
    double Trec = sqrt(Pmiss.Mag2() + Mre*Mre) - Mre;
    Emiss       =  Nu - Trec;
    
    // c.m. momentum
    // Light cone fractions
    alpha_q     = A2mA * ( Nu - q3Vector.Z() );
    SumAlpha    = SumpBackAlpha = -alpha_q;
    Pcm         = -q3Vector ;
    
    
    // Wtilde = (mA - m'(A - Np - 1)  + q - P0 - P1 - P2 - ... ) = invariant mass of the outgoing hadronic system
    Wtilde      = A_4Vector + q;
    
    
    
    // Loop over the protons....
    for (int p = 0 ; p < Np ; p++){
        RotateVector( & Pp.at(p) );
        Pcm        += Pp.at(p);
        p4momentum  .SetVectM( Pp.at(p) , Mp );
        P           .push_back( p4momentum );
        


        // A(e,e'p)X missing energy
        Tp[p]       = P.at(p).E() - Mp;      // kinetic energy of the proton (p)
        Emiss      -= Tp[p];
        
        
        
        // Light cone fractions
        alpha[p]    = A2mA * ( P.at(p).E() - P.at(p).Pz()  )   ;
        SumAlpha   += alpha[p];
        
        
        
        // backward going protons
        theta_pq[p] = RAD2DEG*( Pp.at(p).Angle(q3Vector) );
        if ( theta_pq[p] > 110 ){
            pBack.push_back( p4momentum );
            pBackAlpha[NpBack]   = A2mA * ( pBack.at(NpBack).E() - pBack.at(NpBack).Pz()  )   ;
            SumpBackAlpha       += pBackAlpha[NpBack];
            Wtilde              -= pBack.at(NpBack);
            NpBack ++;
        }
    }
    
    A_Np_1_4Vector .SetVectM( TVector3() , (A - NpBack - 1)*Mp );  // approximation of the A-NpBack-1 ... system mass
    Wtilde      -= A_Np_1_4Vector;

    
    W2tilde     = fabs(Wtilde.Mag2());
    Xbtilde     = Q2 / (W2tilde + Q2 - Mp2) ;
    if (Xbtilde < 0) {
        Printf(" W2tilde = %f , Q2 = %f,Mp2 = %f , W2tilde + Q2 - Mp2 = %f, Q2 / (W2tilde + Q2 - Mp2) = %f"
               ,W2tilde , Q2 , Mp2,W2tilde + Q2 - Mp2,Q2 / (W2tilde + Q2 - Mp2));
    }
    
    counter++;
    OutTree -> Fill();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::RotateVector( TVector3 * V ){
    // move to q-Pmiss system: q is the z axis, Pmiss is in x-z plane: Pmiss=(Pmiss[x],0,Pmiss[q])
    V -> RotateZ(-q_phi);
    V -> RotateY(-q_theta);
    V -> RotateZ(-Pmiss_Phi);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::PrintData(TString Stage){
    timer.Stop();
    std::cout<< "after " << timer.CpuTime()<< " ns,"<< Stage << "(Entry " << Entry << Form(")-%d protons, %d neutrons",Np,N_n) << std::endl;
    //    plot -> Print4Momentum(q,"q");
    for (int p = 0 ; p < Np ; p++){
        //        plot -> Print4Momentum(P.at(p),Form("p(%d)",p));
        SHOW(pEdep[p]);
        //        plot -> PrintOut3Vector(Form("p(%d) vertex",p),pVertex[p]);
    }
    //    plot -> PrintOut3Vector("Pmiss",Pmiss);
    //    plot -> PrintOut3Vector("Pcm",Pcm);
    SHOW(Emiss);
    SHOW(ThetaPQ);
    SHOW(SumAlpha);
    SHOW(alpha[0]);
    SHOW(W2tilde);
    SHOW(Xbtilde);
    std::cout << "--------------------------" << std::endl ;
    timer.Start();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::InitializeInputTree(){
    InTree -> SetBranchAddress("Xb"                 , &Xb);
    InTree -> SetBranchAddress("P_nmb"              , &Np);
    InTree -> SetBranchAddress("N_nmb"              , &N_n);
    if(RawData) { // For the raw data we have the following branches
        // protons vertex
        InTree -> SetBranchAddress("P_X"            , &unsortedX_proton);
        InTree -> SetBranchAddress("P_Y"            , &unsortedY_proton);
        InTree -> SetBranchAddress("P_Z"            , &unsortedZ_proton);
        InTree -> SetBranchAddress("P_Px"           , &PpX);    // protons momenta
        InTree -> SetBranchAddress("P_Py"           , &PpY);
        InTree -> SetBranchAddress("P_Pz"           , &PpZ);
        InTree -> SetBranchAddress("P_CTOF"         , &unsortedCTOF);    // protons CTOF
        InTree -> SetBranchAddress("P_PID"          , &unsortedpPID);    // positive particles momenta
        InTree -> SetBranchAddress("P_cut"          , &unsortedpCut);    // positive particles momenta
        InTree -> SetBranchAddress("P_Edep"         , &unsortedpEdep);
        InTree -> SetBranchAddress("N_Px"           , &N_Px);    // negative particles momenta (electron is the first)
        InTree -> SetBranchAddress("N_Py"           , &N_Py);    // negative particles momenta (electron is the first)
        InTree -> SetBranchAddress("N_Pz"           , &N_Pz);    // negative particles momenta (electron is the first)
    }
    else {
        InTree -> SetBranchAddress("Q2"             , &Q2);
        InTree -> SetBranchAddress("Nu"             , &Nu);
        InTree -> SetBranchAddress("T_nmb"          , &N_t);
        InTree -> SetBranchAddress("Px_e"           , &Px_e);
        InTree -> SetBranchAddress("Py_e"           , &Py_e);
        InTree -> SetBranchAddress("Pz_e"           , &Pz_e);
        InTree -> SetBranchAddress("W"              , &W);
        InTree -> SetBranchAddress("X_e"            , &X_e);
        InTree -> SetBranchAddress("Y_e"            , &Y_e);
        InTree -> SetBranchAddress("Z_e"            , &Z_e);
        InTree -> SetBranchAddress("CTOF"           , &unsortedCTOF);
        InTree -> SetBranchAddress("Px"             , &PpX);
        InTree -> SetBranchAddress("Py"             , &PpY);
        InTree -> SetBranchAddress("Pz"             , &PpZ);
        InTree -> SetBranchAddress("X"              , &unsortedX_proton);
        InTree -> SetBranchAddress("Y"              , &unsortedY_proton);
        InTree -> SetBranchAddress("Z"              , &unsortedZ_proton);
    }
    Nentries    = InTree -> GetEntries();
    counter     = 0;
    std::cout << "Initialized Input InTree TCalcPhysVars, Nentries = " <<  Nentries << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCalcPhysVars::InitializeOutputTree(){
    // Integer branches
    OutTree -> Branch("Np"                  ,&Np                    ,"Np/I");                   // number of protons in event
    OutTree -> Branch("NpBack"              ,&NpBack                ,"NpBack/I");               // number of protons going back
    OutTree -> Branch("N_n"                 ,&N_n                   ,"N_n/I");                  // number of neutrons in event
    
    // TLorentzVector branches
    OutTree -> Branch("P"                   ,&P);                   // all momentum 4-vectors together - in array
    OutTree -> Branch("pBack"               ,&pBack);               // 4-momenta of backward going protons
    OutTree -> Branch("q"                   ,"TLorentzVector"       ,&q);
    OutTree -> Branch("Wtilde"              ,"TLorentzVector"       ,&Wtilde);

    // TVector3 branches
    OutTree -> Branch("Pp"                  ,&Pp);                  // all momenutm 3-vectors together - in array
    OutTree -> Branch("pVertex"             ,&pVertex);             // all vertex 3-vectors together - in array
    OutTree -> Branch("Pmiss"               ,"TVector3"             ,&Pmiss);
    OutTree -> Branch("Pcm"                 ,"TVector3"             ,&Pcm);
    OutTree -> Branch("eVertex"             ,"TVector3"             ,&eVertex);

    // Float_t branches
    OutTree -> Branch("Xb"                  ,&Xb                    , "Xb/F");
    OutTree -> Branch("Q2"                  ,&Q2                    , "Q2/F");
    OutTree -> Branch("pPID"                ,pPID                   , "pPID[Np]/I");
    OutTree -> Branch("pEdep"               ,pEdep                  , "pEdep[Np]/F");
    OutTree -> Branch("pCut"                ,pCut                   , "pCut[Np]/I");
    OutTree -> Branch("pZ"                  ,&pZ                    , "pZ[Np]/F");
    OutTree -> Branch("W"                   ,&W                     , "W/F");
    OutTree -> Branch("pCTOF"               ,&CTOF                  , "pCTOF[Np]/F");
    OutTree -> Branch("Emiss"               ,&Emiss                 , "Emiss/F");
    OutTree -> Branch("PoverQ"              ,&PoverQ                , "PoverQ/F");
    OutTree -> Branch("alpha"               ,&alpha                 , "alpha[Np]/F");            //light cone fraction
    OutTree -> Branch("alpha_q"             ,&alpha_q               , "alpha_q/F");
    OutTree -> Branch("SumAlpha"            ,&SumAlpha              , "SumAlpha/F");            // sum of light cone fractions
    OutTree -> Branch("W2tilde"             ,&W2tilde               , "W2tilde/F");            // sum of light cone fractions
    OutTree -> Branch("Xbtilde"             ,&Xbtilde               , "Xbtilde/F");            // sum of light cone fractions
    OutTree -> Branch("ThetaPQ"             ,&ThetaPQ               , "ThetaPQ/F");
    OutTree -> Branch("theta_pq"            ,&theta_pq              , "theta_pq[Np]/F");
    OutTree -> Branch("pBackAlpha"          ,&pBackAlpha            , "pBackAlpha[NpBack]/F");      //light cone fraction of backward protons
    OutTree -> Branch("SumpBackAlpha"       ,&SumpBackAlpha         , "SumpBackAlpha/F");            // sum of light cone fractions
    
    counter = 0;
    timer.Start();
    std::cout << "Initialized Output Tree TCalcPhysVars" << std::endl;
}







