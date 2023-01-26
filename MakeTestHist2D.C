#include <iostream>
#include <sstream>
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include <TStyle.h>
#include"TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include <ROOT/RDataFrame.hxx>
#include "ana_const.h"
using namespace std;
void MakeTestHist2D(){
	gStyle->SetOptFit(1);
	TH2::SetDefaultSumw2();
	ROOT::RDataFrame dmyDf_mc1_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test1.root"); 
	ROOT::RDataFrame dmyDf_mc2_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test2.root"); 
	ROOT::RDataFrame dmyDf_mc3_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test3.root"); 
	ROOT::RDataFrame dmyDf_mc4_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test4.root"); 

	auto ReWeight = [](float cosThetaCS, float PhiCS, float wt) { return wt*(1+ 0.6*cosThetaCS*cosThetaCS+(2*0.2)*(sqrt(1-cosThetaCS*cosThetaCS)*cosThetaCS*cos(PhiCS))+0.2*(1-cosThetaCS*cosThetaCS)*cos(2*PhiCS)*0.5)/(1+ cosThetaCS*cosThetaCS);};

	auto dmyDf_mc1 = dmyDf_mc1_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});
	auto dmyDf_mc2 = dmyDf_mc2_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});
	auto dmyDf_mc3 = dmyDf_mc3_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});
	auto dmyDf_mc4 = dmyDf_mc4_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});

	cout << "count total events file 1: "<< *dmyDf_mc1.Count()<<endl;
	cout << "count total events file 2: "<< *dmyDf_mc2.Count()<<endl;
	cout << "count total events file 3: "<< *dmyDf_mc3.Count()<<endl;
	cout << "count total events file 4: "<< *dmyDf_mc4.Count()<<endl;
	
	int events= 45000;
	int steps1=*dmyDf_mc1.Count()/events;
	int steps2=*dmyDf_mc2.Count()/events;
	int steps3=*dmyDf_mc3.Count()/events;
	int steps4=*dmyDf_mc4.Count()/events;

	int range1 = *dmyDf_mc1.Count()/steps1;
	int range2 = *dmyDf_mc2.Count()/steps2;
	int range3 = *dmyDf_mc3.Count()/steps3;
	int range4 = *dmyDf_mc4.Count()/steps4;

	cout << "counts in each step1"<< *dmyDf_mc1.Count()/steps1 <<endl;
	cout << "counts in each step2"<< *dmyDf_mc2.Count()/steps2 <<endl;
	cout << "counts in each step3"<< *dmyDf_mc3.Count()/steps3 <<endl;
	cout << "counts in each step3"<< *dmyDf_mc4.Count()/steps4 <<endl;
	TFile *myfile = new TFile("testData2D.root","recreate");

	////For the first file
	int split=0;
	for (int i_split = 1; i_split <=steps1 ; i_split++){	
		int begin_= 1+(i_split-1)*range1;
		int end_= range1* i_split;
		cout << "begin event: "<< begin_ << "end event: "<< end_<<endl;

		auto mc_reco = dmyDf_mc1.Range(begin_,end_).Histo2D({"phi_theta_mc_4pi", "phi_theta_mc_4pi", nBinPhi,phiMin, phiMax, nBinCosth, costhMin, costhMax},"phi", "costh","ParReweight");
		char name[100];
		cout << "file number in root: "<< i_split << endl;
		sprintf(name, "%s%i","costh_test_",i_split);
		TH2D* h = (TH2D*)mc_reco->Clone();    
		h->SetName(name); h->Write(name,TObject::kWriteDelete);
	}
	//////for the second file
	for (int i_split = 1; i_split <=steps2 ; i_split++){    
		int begin_= 1+(i_split-1)*range2;
		int end_= range2* i_split;
		cout << "begin event: "<< begin_ << "end event: "<< end_<<endl;
		auto mc_reco = dmyDf_mc2.Range(begin_,end_).Histo2D({"phi_theta_mc_4pi", "phi_theta_mc_4pi", nBinPhi,phiMin, phiMax, nBinCosth, costhMin, costhMax},"phi", "costh","ParReweight");



		char name[100];
		cout << "file number in root: "<< i_split+steps1 << endl;
		sprintf(name, "%s%i","costh_test_",i_split+steps1);
		TH2D* h = (TH2D*)mc_reco->Clone();        
		h->SetName(name); h->Write(name,TObject::kWriteDelete);
	}

	//////for the third file
	for (int i_split = 1; i_split <=steps3 ; i_split++){    
		int begin_= 1+(i_split-1)*range3;
		int end_= range3* i_split;
		cout << "begin event: "<< begin_ << "end event: "<< end_<<endl;
		auto mc_reco = dmyDf_mc3.Range(begin_,end_).Histo2D({"phi_theta_mc_4pi", "phi_theta_mc_4pi", nBinPhi,phiMin, phiMax, nBinCosth, costhMin, costhMax},"phi", "costh","ParReweight");
		char name[100];
		cout << "file number in root: "<< i_split+steps1+steps2 << endl;
		sprintf(name, "%s%i","costh_test_",i_split+steps1+steps2);
		TH2D* h = (TH2D*)mc_reco->Clone();    
		h->SetName(name); h->Write(name,TObject::kWriteDelete);
	}

	 //////for the 4th file
        for (int i_split = 1; i_split <=steps4 ; i_split++){    
                int begin_= 1+(i_split-1)*range4;
                int end_= range4* i_split;
                cout << "begin event: "<< begin_ << "end event: "<< end_<<endl;
                auto mc_reco = dmyDf_mc4.Range(begin_,end_).Histo2D({"phi_theta_mc_4pi", "phi_theta_mc_4pi", nBinPhi,phiMin, phiMax, nBinCosth, costhMin, costhMax},"phi", "costh","ParReweight");
                char name[100];
                cout << "file number in root: "<< i_split+steps1+steps2+ steps3 << endl;
                sprintf(name, "%s%i","costh_test_",i_split+steps1+steps2+steps3);
                TH2D* h = (TH2D*)mc_reco->Clone();    
                h->SetName(name); h->Write(name,TObject::kWriteDelete);
        } 
	//Closing the file
	myfile->Close();
	myfile->ls();
}
