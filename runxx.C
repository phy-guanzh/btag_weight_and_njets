#include "WWg.C"
#include <iostream>
#include <fstream>
using namespace std;

void runxx(TString year, TString dir) {
ifstream infile("file"+year);
string buffer; 
TString infilename;

int k=1;

while (k>0){
getline (infile, buffer) ;
infilename = buffer;
if(infilename.Contains("root")==0) {k=-2; continue;}
if(infilename.Contains("end")==1) {k=-2; break;}
//TString filename ="cutlep-"+infilename;
TString filename ="cutla-"+infilename;
TString outname="cutla-"+infilename;

cout<<dir<<infilename<<" -> "<<outname<<endl;

TFile *file1 =new TFile(dir+filename);
TTree *tree1 = (TTree*) file1->Get("Events");
WWg m1(tree1,outname);
cout<<outname<<endl;
m1.Loop(year,outname);
m1.endJob();
 
}
}

int main(){
	//TString dir16 ="/home/pku/anying/cms/rootfiles/WWg/2016/";
	//TString dir17 ="/home/pku/anying/cms/rootfiles/WWg/2017/";
	//TString dir18 ="/home/pku/guanz/cms/updatejetmet/";
	TString dir18 ="/home/pku/guanz/cms/updatejetmet/";
	TString dir17 ="/home/pku/guanz/cms/updatejetmet/";
	TString dir16 ="/home/pku/guanz/cms/updatejetmet/";
//	TString dir16 ="~/gitsig/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/WWG/MVA/2016post/WWG_selector/";
	//TString dir18 ="~/gitsig/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/WWG/MVA/2016post/WWG_selector/";
	//TString dir17 ="~/gitsig/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/WWG/MVA/2016post/WWG_selector/17/";
	//TString dir17 ="/home/pku/guanz/cms/mva_rootfiles/";
	//TString dir16 ="/home/pku/guanz/cms/all_mva_rootfiles/";
	//runxx("18",dir18);
	//runxx("17",dir17);
	//runxx("16",dir16);
	runxx("16pre",dir16);
//	runxx("16",dir16);
	return 1;
}

