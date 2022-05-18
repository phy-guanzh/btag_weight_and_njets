#include "WWg.h"
float WWg::get_btag_scale(TH1D*hbeff,TH1D*hceff,TH1D*hleff,double value,int nJet, float Jet_btagDeepFlavB[30],int Jet_hadronFlavour[30], float Jet_btagSF_deepjet[30],float Jet_btagSF_deepjet_up[30],float Jet_btagSF_deepjet_down[30],float Jet_pt[30], float Jet_eta[30],TString type){
	float MCweight_b=1;
	float MCweight_c=1;
	float MCweight_light=1;
	float Dataweight_b=1;
	float Dataweight_c=1;
	float Dataweight_light=1;
	float Dataweight_b_up=1,Dataweight_light_up=1;
	float Dataweight_c_up=1,Dataweight_b_down=1;
	float Dataweight_c_down=1,Dataweight_light_down=1;
	float btag_weight_nominal_b=1,btag_weight_nominal_c=1,btag_weight_nominal_light=1,btag_weight_light_up=1,btag_weight_b_up=1,btag_weight_c_up=1,btag_weight_light_down=1,btag_weight_b_down=1,btag_weight_c_down=1;
	for(int i=0;i<nJet;i++){
         //if  (Jet_pt[i]<20 || abs(Jet_eta[i]> 2.4)) continue;
		if(fabs(Jet_hadronFlavour[i])==5){
			if(Jet_btagDeepFlavB[i]>value){
				MCweight_b=MCweight_b*hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]));
				Dataweight_b=Dataweight_b*hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet[i];
				Dataweight_b_up=Dataweight_b_up*hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_up[i];
				Dataweight_b_down=Dataweight_b_down*hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_down[i];
                //if (Jet_pt[i]< 10 )cout << "beff: pt<10   "<<hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i])) <<endl;
                //if (Jet_pt[i]< 20 )cout << "beff: pt<20   "<<hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i])) <<endl;
                if (abs(Jet_eta[i]> 2.4) )cout << "beff:abs_eta>2.4   "<<hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i])) <<endl;
                if (abs(Jet_eta[i]> 2.5) )cout << "beff: abs_eta>2.5   "<<hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i])) <<endl;
			}
			else{
				MCweight_b=MCweight_b*(1-hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i])));
				Dataweight_b=Dataweight_b*(1-hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet[i]);
				Dataweight_b_up=Dataweight_b_up*(1-hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_up[i]);
				Dataweight_b_down=Dataweight_b_down*(1-hbeff->GetBinContent(hbeff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_down[i]);
			}
		}
		else if(fabs(Jet_hadronFlavour[i])==4){
			if(Jet_btagDeepFlavB[i]>value){
				MCweight_c=MCweight_c*hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]));
				Dataweight_c=Dataweight_c*hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet[i];
				Dataweight_c_up=Dataweight_c_up*hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_up[i];
				Dataweight_c_down=Dataweight_c_down*hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_down[i];
			}
			else{
				MCweight_c=MCweight_c*(1-hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i])));
				Dataweight_c=Dataweight_c*(1-hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet[i]);
				Dataweight_c_up=Dataweight_c_up*(1-hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_up[i]);
				Dataweight_c_down=Dataweight_c_down*(1-hceff->GetBinContent(hceff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_down[i]);
			}
		}
		else{
			if(Jet_btagDeepFlavB[i]>value){
				MCweight_light=MCweight_light*hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]));
				Dataweight_light=Dataweight_light*hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet[i];
				Dataweight_light_up=Dataweight_light_up*hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_up[i];
				Dataweight_light_down=Dataweight_light_down*hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_down[i];
			}
			else{
				MCweight_light=MCweight_light*(1-hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i])));
		    	Dataweight_light=Dataweight_light*(1-hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet[i]);
				Dataweight_light_up=Dataweight_light_up*(1-hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_up[i]);
				Dataweight_light_down=Dataweight_light_down*(1-hleff->GetBinContent(hleff->FindBin(Jet_pt[i],Jet_eta[i]))*Jet_btagSF_deepjet_down[i]);
			}
		}
	}
	if(1>0){
        //cout << "MCweight_b " <<MCweight_b << " Dataweight_b " << Dataweight_b <<" Dataweight_b_up "<<Dataweight_b_up << " Dataweight_b_down "<< Dataweight_b_down <<endl;
		btag_weight_nominal_b=Dataweight_b/MCweight_b;
		btag_weight_nominal_c=Dataweight_c/MCweight_c;
		btag_weight_nominal_light=Dataweight_light/MCweight_light;
		btag_weight_light_up=Dataweight_light_up/MCweight_light;
		btag_weight_light_down=Dataweight_light_down/MCweight_light;
		btag_weight_b_up=Dataweight_b_up/MCweight_b;
		btag_weight_b_down=Dataweight_b_down/MCweight_b;
		btag_weight_c_up=Dataweight_c_up/MCweight_c;
		btag_weight_c_down=Dataweight_c_down/MCweight_c;
        if (btag_weight_nominal_b !=btag_weight_nominal_b ) btag_weight_nominal_b=1;
        if (btag_weight_nominal_light !=btag_weight_nominal_light ) btag_weight_nominal_light=1;
        if (btag_weight_nominal_c !=btag_weight_nominal_c ) btag_weight_nominal_c=1;
        if (btag_weight_light_up !=btag_weight_light_up ) btag_weight_light_up =1;
        if (btag_weight_light_down !=btag_weight_light_down ) btag_weight_light_down =1;
        if (btag_weight_c_up !=btag_weight_c_up ) btag_weight_c_up =1;
        if (btag_weight_c_down !=btag_weight_c_down ) btag_weight_c_down =1;
        if (btag_weight_b_up !=btag_weight_b_up ) btag_weight_b_up =1;
        if (btag_weight_b_down !=btag_weight_b_down ) btag_weight_b_down =1;
}
	if(type.Contains("bup"))
		return btag_weight_b_up;
	else if(type.Contains("bdown"))
		return btag_weight_b_down;
	else if(type.Contains("cup"))
		return btag_weight_c_up;
	else if(type.Contains("cdown"))
		return btag_weight_c_down;
	else if(type.Contains("lightup"))
		return btag_weight_light_up;
	else if(type.Contains("lightdown"))
		return btag_weight_light_down;
	else if (type.Contains("nominal_b"))
		return btag_weight_nominal_b;
	else if (type.Contains("nominal_c"))
		return btag_weight_nominal_c;
	else if (type.Contains("nominal_light"))
		return btag_weight_nominal_light;
    else 
        return 0;
}

