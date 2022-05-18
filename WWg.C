#define WWg_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TLorentzVector.h"
#include "get_btag_scale.C"
#include "math.h"
using namespace std;
void WWg::Loop(TString year,TString name)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TString y;
   if(year=="16") y="16post";
   else y=year;
   TFile*f1 = new TFile("/home/pku/guanz/HLT/HLT_mva_less/results/HLTSF_20"+y+".root");
   TH1D*h1 = (TH1D*)f1->Get("HLTSF_20"+y);
   TH1D*h2 = (TH1D*)f1->Get("HLTSF_UP_ERR_20"+y);
   TH1D*h3 = (TH1D*)f1->Get("HLTSF_DOWN_ERR_20"+y);
   

   TFile*f2 = new TFile("~/cms/btag_eff/"+y+"/btag_eff_20"+y+".root");
   TH1D*hbeff = (TH1D*)f2->Get("h_eff_b");
   TH1D*hceff = (TH1D*)f2->Get("h_eff_c");
   TH1D*hleff = (TH1D*)f2->Get("h_eff_light");
   Long64_t nbytes = 0, nb = 0;
   Bool_t BSL=0,LEP=0,PHOTON=0,missET=0;
   int tot=0;double pfmet,pfmetphi,puppimet,puppimetphi,phiVlep;
   TLorentzVector Zp4, photonp4, jet1p4, jet2p4,lep1p4,lep2p4;
   HLT_SF=1;HLT_SF_Up=1;HLT_SF_Down=1;
   //Btag_Weight = 1; 
   float test=1,test1=1;
      UInt_t newb=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout << name.Contains("fakeL")<<endl;
      bool isdata=0,isMC=0;
      UInt_t newb=0;
      float pmc1=1;
      float pmc2=1;
      float pdata1=1;
      float pdata2=1;
      float pdata1_up=1;
      njet_30_jesup=0;njet_30_jesdown=0;nbjet_20_jesup=0;nbjet_20_jesdown=0;njet_30_jerup=0;njet_30_jerdown=0;nbjet_20_jerup=0;nbjet_20_jerdown=0;njet_30_test=0;
      float pdata2_up=1;
      float pdata1_down=1;
      float pdata2_down=1;
      float pmc=1;
      float pdata=1;
      float pdata_up =1;
      float pdata_down=1;
      float pi=3.1415926;
      float deta,dphi;
      Btag_Weight = 1; 
      HLT_SF=1;HLT_SF_Up=1;HLT_SF_Down=1;
      if( (name.Contains("Muon")==0 && name.Contains("Ele")==0 && name.Contains("MET")==0 && name.Contains("EGamma")==0 &&  name.Contains("plj")==0 && name.Contains("fakeL")==0) ) isMC=1;
      else isdata=1;
           double value,value1;
           if(year.Contains("18")){ value=0.2783;value1=0.1208;}
           else if(year.Contains("17")){ value=0.3040;value1=0.1355;}
           else if(year.Contains("16")){
		   if(name.Contains("pre")){ value=0.2598;value1=0.2027;}
		   else{value=0.2489;value1=0.1918;}
	   }
       //if( isMC ){
       if( isMC ){
           cout<<jentry<<" "<<nentries <<endl;
		   btag_weight_medium_b     = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"nominal_b");
		   btag_weight_medium_c      = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"nominal_c");
		   btag_weight_medium_light      = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"nominal_light");
		   btag_weight_medium_b_up   = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"bup");
		   btag_weight_medium_c_up   = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"cup");
		   btag_weight_medium_light_up   = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"lightup");
	   btag_weight_medium_b_down = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"bdown");
	   btag_weight_medium_c_down = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"cdown");
	   btag_weight_medium_light_down = get_btag_scale(hbeff,hceff,hleff,value,nJet,Jet_btagDeepFlavB,Jet_hadronFlavour,Jet_btagSF_deepjet_M,Jet_btagSF_deepjet_M_up,Jet_btagSF_deepjet_M_down,Jet_pt,Jet_eta,"lightdown");
    }
  
  
     /* 
        cout << "Btag weight B now: "<<btag_weight_medium_b <<" " << btag_weight_medium_light <<endl;
        cout << "Btag weight C now: "<<btag_weight_medium_c <<" " << btag_weight_medium_light <<endl;
        cout << "Btag weight light up  now: "<<btag_weight_medium_light_up <<endl;
        cout << "Btag weight light down now: "<<btag_weight_medium_light_down <<endl;
        cout << "Btag weight b up  now: "<<btag_weight_medium_b_up <<endl;
        cout << "Btag weight c up  now: "<<btag_weight_medium_c_up <<endl;
        cout << "Btag weight b down now: "<<btag_weight_medium_b_down <<endl;
        cout << "Btag weight c down now: "<<btag_weight_medium_c_down <<endl;
        cout << "-------------------" <<endl;  
      */
      //cout << "in this event there are " << nJet << "Jets!!!" <<endl;
      //cout << "in this event 1th Jet btagSF is " << Jet_btagSF_deepjet_shape[1] << "!!!" <<endl;
      //
      //
      if (isMC){
     while ( newb <= nJet ) 
     { 
        if (abs(Jet_eta[newb]) > 4.7) 
        {newb+=1;continue; }
        deta=Jet_eta[newb]-photoneta ;
        dphi=Jet_phi[newb]-photonphi;
        //cout << deta << " " <<dphi << " newb: "<<newb<<endl;
        if (dphi > pi/2)  dphi -= 2 * pi;
        else if (dphi < -pi/2) dphi += 2 * pi;
       // cout << sqrt(deta*deta+dphi*dphi) <<endl;
        if( sqrt(deta*deta+dphi*dphi) < 0.4 )
        {newb+=1;continue; }
        
        deta=Jet_eta[newb]-lep1eta ;
        dphi=Jet_phi[newb]-lep1phi;
        if (dphi > pi/2)  dphi -= 2 * pi;
        else if (dphi < -pi/2) dphi += 2 * pi;
        if( sqrt(deta*deta+dphi*dphi) < 0.4 ) 
        {newb+=1;continue; }

        deta=Jet_eta[newb]-lep2eta ;
        dphi=Jet_phi[newb]-lep2phi;
        if (dphi > pi/2)  dphi -= 2 * pi;
        else if (dphi < -pi/2) dphi += 2 * pi;
        if( sqrt(deta*deta+dphi*dphi) < 0.4 )
        {newb+=1;continue; }
        
        //cout << Jet_jetId[newb] <<" update: "  << (Jet_jetId[newb] >> 1 & 1) <<endl;
         if (Jet_jetId[newb] >> 1 & 1){
          if (y=="16post") {
         if (Jet_pt_jesTotalUp[newb] > 30) {
             njet_30_jesup = njet_30_jesup+1;
             }
         if (Jet_pt_jesTotalDown[newb] > 30) {
             njet_30_jesdown = njet_30_jesdown+1;
             }
         if (y=="16post" && Jet_btagDeepFlavB[newb] >0.2489 && Jet_pt_jesTotalUp[newb] > 20){
               nbjet_20_jesup = nbjet_20_jesup+1;   
         }
         if (y=="16post" && Jet_btagDeepFlavB[newb] >0.2489 && Jet_pt_jesTotalDown[newb] > 20){
               nbjet_20_jesdown = nbjet_20_jesdown+1;   
         }
         
         if (Jet_pt_jerUp[newb] > 30) {
             njet_30_jerup = njet_30_jerup+1;
             }
         if (Jet_pt_jerDown[newb] > 30) {
             njet_30_jerdown = njet_30_jerdown+1;
             }
         if (y=="16post" && Jet_btagDeepFlavB[newb] >0.2489 && Jet_pt_jerUp[newb] > 20){
               nbjet_20_jerup = nbjet_20_jerup+1;   
         }
         if (y=="16post" && Jet_btagDeepFlavB[newb] >0.2489 && Jet_pt_jesTotalDown[newb] > 20){
               nbjet_20_jerdown = nbjet_20_jerdown+1;   
         }


         }
          if (y=="16pre") {
         if (Jet_pt_jesTotalUp[newb] > 30) {
             njet_30_jesup = njet_30_jesup+1;
             }
         if (Jet_pt_jesTotalDown[newb] > 30) {
             njet_30_jesdown = njet_30_jesdown+1;
             }
         if (y=="16pre" && Jet_btagDeepFlavB[newb] >0.2598 && Jet_pt_jesTotalUp[newb] > 20){
               nbjet_20_jesup = nbjet_20_jesup+1;   
         }
         if (y=="16pre" && Jet_btagDeepFlavB[newb] >0.2598 && Jet_pt_jesTotalDown[newb] > 20){
               nbjet_20_jesdown = nbjet_20_jesdown+1;   
         }
         
         if (Jet_pt_jerUp[newb] > 30) {
             njet_30_jerup = njet_30_jerup+1;
             }
         if (Jet_pt_jerDown[newb] > 30) {
             njet_30_jerdown = njet_30_jerdown+1;
             }
         if (y=="16pre" && Jet_btagDeepFlavB[newb] >0.2598 && Jet_pt_jerUp[newb] > 20){
               nbjet_20_jerup = nbjet_20_jerup+1;   
         }
         if (y=="16pre" && Jet_btagDeepFlavB[newb] >0.2598 && Jet_pt_jerDown[newb] > 20){
               nbjet_20_jerdown = nbjet_20_jerdown+1;   
         }
}
          if (y=="17") {
         if (Jet_pt_jesTotalUp[newb] > 30) {
             njet_30_jesup = njet_30_jesup+1;
             }
         if (Jet_pt_jesTotalDown[newb] > 30) {
             njet_30_jesdown = njet_30_jesdown+1;
             }
         if (y=="17" && Jet_btagDeepFlavB[newb] >0.3040 && Jet_pt_jesTotalUp[newb] > 20){
               nbjet_20_jesup = nbjet_20_jesup+1;   
         }
         if (y=="17" && Jet_btagDeepFlavB[newb] >0.3040 && Jet_pt_jesTotalDown[newb] > 20){
               nbjet_20_jesdown = nbjet_20_jesdown+1;   
         }
         
         if (Jet_pt_jerUp[newb] > 30) {
             njet_30_jerup = njet_30_jerup+1;
             }
         if (Jet_pt_jerDown[newb] > 30) {
             njet_30_jerdown = njet_30_jerdown+1;
             }
         if (y=="17" && Jet_btagDeepFlavB[newb] >0.3040 && Jet_pt_jerUp[newb] > 20){
               nbjet_20_jerup = nbjet_20_jerup+1;   
         }
         if (y=="17" && Jet_btagDeepFlavB[newb] >0.3040 && Jet_pt_jerDown[newb] > 20){
               nbjet_20_jerdown = nbjet_20_jerdown+1;   
         }

}
          if (y=="18") {
         if (Jet_pt[newb] > 30) {
             njet_30_test = njet_30_test+1;
             }
         if (Jet_pt_jesTotalUp[newb] > 30) {
             njet_30_jesup = njet_30_jesup+1;
             }
         if (Jet_pt_jesTotalDown[newb] > 30) {
             njet_30_jesdown = njet_30_jesdown+1;
             }
         if (y=="18" && Jet_btagDeepFlavB[newb] >0.2783 && Jet_pt_jesTotalUp[newb] > 20){
               nbjet_20_jesup = nbjet_20_jesup+1;   
         }
         if (y=="18" && Jet_btagDeepFlavB[newb] >0.2783 && Jet_pt_jesTotalDown[newb] > 20){
               nbjet_20_jesdown = nbjet_20_jesdown+1;   
         }
         
         if (Jet_pt_jerUp[newb] > 30) {
             njet_30_jerup = njet_30_jerup+1;
             }
         if (Jet_pt_jerDown[newb] > 30) {
             njet_30_jerdown = njet_30_jerdown+1;
             }
         if (y=="18" && Jet_btagDeepFlavB[newb] >0.2783 && Jet_pt_jerUp[newb] > 20){
               nbjet_20_jerup = nbjet_20_jerup+1;   
         }
         if (y=="18" && Jet_btagDeepFlavB[newb] >0.2783 && Jet_pt_jerDown[newb] > 20){
               nbjet_20_jerdown = nbjet_20_jerdown+1;   
         }
      }
     }
      newb += 1;
     }
      }
      cout << "njets_diff:" << (njet_30_test- njets30_pc) << " njet_30_jesup " << njet_30_jesup<<endl; 
      LEP = channel==1 && fabs(lep1_pid)==13 && fabs(lep2_pid)==11 && lep1_charge*lep2_charge<0 && drll>0.5 && lep1pt>15 && lep2pt>15 && fabs(lep1eta) < 2.4 && fabs(lep2eta) < 2.5 && n_loose_ele==1 && n_loose_mu==1 && ptll>30 && mll>20;
      PHOTON = n_photon>0 && photonet > 20. && ( (fabs(photoneta) < 1.4442) || ( fabs(photoneta)<2.5 && fabs(photoneta)>1.566 ) ) && drl1a>0.5 && drl2a>0.5;
      if(name.Contains("TTJets")) PHOTON =  PHOTON && photon_isprompt==1 && photon_gen_matching==4 && lep1_isprompt && lep2_isprompt;
      if( !( LEP && PHOTON  ) )
	      continue;

      lep1p4.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.105666);
      lep2p4.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.000511);
      phiVlep=(lep1p4+lep2p4).Phi();
      photonp4.SetPtEtaPhiM(photonet, photoneta, photonphi,0);

      if( name.Contains("Muon")||name.Contains("Ele")||name.Contains("fakeL") ||name.Contains("plj") ){
	      pfmet=MET_T1_pt;pfmetphi=MET_T1_phi;
      }
      else{
	      pfmet=MET_T1Smear_pt;pfmetphi=MET_T1Smear_phi;
      }
      puppimet=PuppiMET_pt;puppimetphi=PuppiMET_phi;

      mT_pf=sqrt(2*(ptll*pfmet*(1-cos(phiVlep-pfmetphi) ) ) );
      mT_puppi=sqrt(2*(ptll*puppimet*(1-cos(phiVlep-puppimetphi) ) ) );

      if (isMC){



      if(lep1pt>lep2pt){
	      mT2_pf=sqrt(2*(lep2pt*pfmet*(1-cos(lep2phi-pfmetphi) ) ) );
	      mT2_puppi=sqrt(2*(lep2pt*puppimet*(1-cos(lep2phi-puppimetphi) ) ) );
	      ml1g=(lep1p4+photonp4).M();
	      ml2g=(lep2p4+photonp4).M();
              HLT_SF=h1->GetBinContent(h1->FindBin(lep1pt,lep2pt));
              HLT_SF_Up=h2->GetBinContent(h2->FindBin(lep1pt,lep2pt));
              HLT_SF_Down=h3->GetBinContent(h3->FindBin(lep1pt,lep2pt));
      }
      else{
	      mT2_pf=sqrt(2*(lep1pt*pfmet*(1-cos(lep1phi-pfmetphi) ) ) );
	      mT2_puppi=sqrt(2*(lep1pt*puppimet*(1-cos(lep1phi-puppimetphi) ) ) );
	      ml1g=(lep2p4+photonp4).M();
	      ml2g=(lep1p4+photonp4).M();
              HLT_SF=h1->GetBinContent(h1->FindBin(lep2pt,lep1pt));
              HLT_SF_Up=h2->GetBinContent(h2->FindBin(lep2pt,lep1pt));
              HLT_SF_Down=h3->GetBinContent(h3->FindBin(lep2pt,lep1pt));
      }} 

      if (isdata){
      if(lep1pt>lep2pt){
	      mT2_pf=sqrt(2*(lep2pt*pfmet*(1-cos(lep2phi-pfmetphi) ) ) );
	      mT2_puppi=sqrt(2*(lep2pt*puppimet*(1-cos(lep2phi-puppimetphi) ) ) );
	      ml1g=(lep1p4+photonp4).M();
	      ml2g=(lep2p4+photonp4).M();
              HLT_SF=1;
              HLT_SF_Up=h2->GetBinContent(h2->FindBin(lep1pt,lep2pt));
              HLT_SF_Down=h3->GetBinContent(h3->FindBin(lep1pt,lep2pt));
      }
      else{
	      mT2_pf=sqrt(2*(lep1pt*pfmet*(1-cos(lep1phi-pfmetphi) ) ) );
	      mT2_puppi=sqrt(2*(lep1pt*puppimet*(1-cos(lep1phi-puppimetphi) ) ) );
	      ml1g=(lep2p4+photonp4).M();
	      ml2g=(lep1p4+photonp4).M();
              HLT_SF=1;
              HLT_SF_Up=h2->GetBinContent(h2->FindBin(lep2pt,lep1pt));
              HLT_SF_Down=h3->GetBinContent(h3->FindBin(lep2pt,lep1pt));
      }} 
      //if(jentry%5000==0) cout<< name<<"isdata?:"<< isdata <<" - "<<jentry<<" "<<nentries<<" "<<scalef<<" "<<HLT_SF<<endl;

      ExTree->Fill();
      tot++;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<nentries<<" "<<tot<<endl;
   f1->Close();
}
