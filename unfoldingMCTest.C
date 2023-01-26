#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include <TRandom3.h>
#include "ana_const.h"
//if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
//endif
//R__LOAD_LIBRARY(../../RooUnfold/libRooUnfold.so)
R__LOAD_LIBRARY(/Users/forhadhossain/test/RooUnfold/build/libRooUnfold.dylib)
using namespace std;
using std::cout;
using std::endl;
void unfoldingMCTest(){
	//gStyle->SetStatFont(63);
	//gStyle->SetStatFontSize(18);
	//gROOT->ForceStyle();
	//gStyle->SetStatW(.6);
	//gStyle->SetStatH(.5);
	//gStyle->SetPalette(kBird);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	gStyle->SetMarkerStyle(20);
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	ROOT::EnableImplicitMT();
	//loading all root files:
	ROOT::RDataFrame Df_mc1("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test1.root");
	ROOT::RDataFrame Df_train("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Train.root");
	ROOT::RDataFrame Df_4pi("result_mc", "data/DY_4pi_GMC_4piAccfactor.root");

	//Reconstructed cut for the MC
        string passed = cutRecoMC;
        cout << "cut: "<< passed << endl;
        auto Df_train2 = Df_train.Define("passed", passed);
	//Getting Acceptance Factor:
	auto hist2d_mc_4pi = Df_4pi.Histo2D({"mc_4pi", "mc_4pi", nBinPhi,phiMin, phiMax, nBinCosth, costhMin, costhMax},"true_phi", "true_costh","weight");
	auto hist2d_mc_reco = Df_mc1.Histo2D({"mc_reco", "mc_reco", nBinPhi,phiMin, phiMax, nBinCosth, costhMin, costhMax},"phi", "costh","weight");
	hist2d_mc_4pi->Scale(1./hist2d_mc_4pi->Integral()); 
	TH2*acc_factor = (TH2D*)hist2d_mc_reco->Clone();
	acc_factor->Divide(hist2d_mc_4pi.GetPtr());
	acc_factor->Scale(1./acc_factor->Integral()); 
	//Building the response matrix:
	TH2D* phi_theta_gen = new TH2D("phi_theta_gen", "MC Generated #phi_{CS} Vs. cos#theta_{CS};#phi_{CS};cos#theta_{CS}",nBinPhi,phiMin,phiMax,nBinCosth,costhMin,costhMax);
	TH2D* phi_theta_reco = new TH2D("phi_theta_reco", "MC Reconstructed  #phi_{CS} Vs. cos#theta_{CS};#phi_{CS};cos#theta_{CS}", nBinPhi,phiMin,phiMax,nBinCosth,costhMin,costhMax);
	RooUnfoldResponse response2D(phi_theta_reco, phi_theta_gen);
	Df_train2.Foreach([&response2D](float reco_phi, float reco_costh, float tr_phi, float tr_costh, float wt, bool pass, int trig){ if(pass==true){response2D.Fill(reco_phi,reco_costh,tr_phi,tr_costh,wt); } else{ response2D.Miss(tr_phi,tr_costh,wt);} }, {"phi","costh","true_phi","true_costh","weight" , "passed", "fpga1"} );

	auto* R = response2D.HresponseNoOverflow();
	auto* c5 = new TCanvas();
	R->SetStats(0);
	R->Draw("colztext");
	c5->SaveAs("plot/testresponse.png");
	bool fitOK = false;
	float cos_low=-0.5; float cos_high=0.5;
	TF2* fit2D = new TF2("fit2D", "[0] * ( 1 + [1]*y*y + 2*[2]*sqrt(1-y*y)*y*cos(x) + [3]*(1-y*y)*cos(2*x)/2.) ", -M_PI, M_PI,cos_low, cos_high);
	fit2D->SetParNames("A", "#lambda","#mu","#nu");
	fit2D->SetParameters(1,1,0,0);


	TF2* fit2D_bbb = new TF2("fit2D_bbb", "[0] * ( 1 + [1]*y*y + 2*[2]*sqrt(1-y*y)*y*cos(x) + [3]*(1-y*y)*cos(2*x)/2.) ", -M_PI, M_PI,cos_low, cos_high);
	fit2D_bbb->SetParNames("A", "#lambda","#mu","#nu");
	fit2D_bbb->SetParameters(1,1,0,0);

	float fit_limit=.5;
	TF1* fcos = new TF1("fcos", "[0] * ( 1 + [1] * x*x ) ",-fit_limit,fit_limit);
	int NTest=54;
	TH2D* reco_test[NTest];
	TH2D* unfolded_bbb_[NTest];
	TFile *acc_file = TFile::Open("testData2D.root","read");
	//TFile *acc_file = TFile::Open("testData67.root","read"); // real data
	int nitr=10;

	float nu_high=0.5;
	TH1D* test_lambda_bayes = new TH1D("test_lambda_bayes", "test_lambda_bayes", 20, -0.5, 1.8);
	TH1D* test_lambda_bbb = new TH1D("test_lambda_bbb", "test_lambda_bbb", 20, -0.5, 1.8);

	TH1D* test_mu_bayes = new TH1D("test_mu_bayes", "test_mu_bayes", 20, -0.3, nu_high);
	TH1D* test_mu_bbb = new TH1D("test_mu_bbb", "test_mu_bbb", 20, -0.3, nu_high);

	TH1D* test_nu_bayes = new TH1D("test_nu_bayes", "test_nu_bayes", 20, -0.3, nu_high);
	TH1D* test_nu_bbb = new TH1D("test_nu_bbb", "test_nu_bbb", 20, -0.3, nu_high);


	TCanvas * c1= new TCanvas("c1","c1", 1800,1200);
	c1->Divide(3,2);

	for(int k=1; k<=NTest; k++){

		TH1D* lambda_bayes = new TH1D("lambda_bayes", "Unfolded #lambda results in different iterations; Iterations; #lambda",nitr,0.5,nitr+0.5);
		TH1D* mu_bayes = new TH1D("mu_bayes", "Unfolded #mu results in different iterations; Iterations; #mu",nitr,0.5,nitr+0.5);
		TH1D* nu_bayes = new TH1D("nu_bayes", "Unfolded #nu results in different iterations; Iterations; #nu",nitr,0.5,nitr+0.5);
		TH1D* chi2NDF_bayes = new TH1D("chi2NDF_bayes","Unfolded #chi^2/NDF results in different iterations; Iterations; #Unfolded #chi^2/NDF",nitr,0.5,nitr+0.5);
		TH1D* lambda_bbb = new TH1D("lambda_bbb", "Unfolded #lambda results in different iterations; Iterations; #lambda",nitr,0.5,nitr+0.5);
		TH1D* mu_bbb= new TH1D("mu_bbb", "Unfolded #mu results in different iterations; Iterations; #mu",nitr,0.5,nitr+0.5);
		TH1D* nu_bbb = new TH1D("nu_bbb", "Unfolded #nuresults in different iterations; Iterations; #nu",nitr,0.5,nitr+0.5);
		TH1D* chi2NDF_bbb = new TH1D("chi2NDF_bbb","Unfolded #chi^2/NDF results in different iterations; Iterations; #Unfolded #chi^2/NDF",nitr,0.5,nitr+0.5);

		reco_test[k-1] = (TH2D*)acc_file->Get(Form("costh_test_%i",k));
		for(int ii =2; ii<=nitr; ii++){
			c5->cd();
			reco_test[k-1]->Draw("colz");
			c5->SaveAs("plot/recoT.png");
			RooUnfoldBayes unfoldBayes(&response2D, reco_test[k-1],ii);
			//RooUnfoldBinByBin unfoldBayes(&response2D, reco_test[k-1]);
			TH2D* unfolded_bayes = (TH2D*) unfoldBayes.Hunfold();
			float  tempChisq=0.0;
			int iter=0;
			fitOK=false;
			while(!fitOK)
			{   
				iter = iter+1;
				unfolded_bayes->Fit("fit2D","R");
				tempChisq = fit2D->GetChisquare() / fit2D->GetNDF();
				fitOK = (tempChisq < 2.0 || iter >10 ) ? true : false;
			}

			if(tempChisq < 2.0 && ii==4 ) { 
				test_lambda_bayes->Fill(fit2D->GetParameter(1));
				test_mu_bayes->Fill(fit2D->GetParameter(2));
				test_nu_bayes->Fill(fit2D->GetParameter(3));
			}  

			lambda_bayes->SetBinContent(ii,fit2D->GetParameter(1));
			mu_bayes->SetBinContent(ii,fit2D->GetParameter(2));
			nu_bayes->SetBinContent(ii,fit2D->GetParameter(3));

			lambda_bayes->SetBinError(ii,fit2D->GetParError(1));
			mu_bayes->SetBinError(ii,fit2D->GetParError(2));
			nu_bayes->SetBinError(ii,fit2D->GetParError(3));
			chi2NDF_bayes->SetBinContent(ii, fit2D->GetChisquare() / fit2D->GetNDF());
			delete  unfolded_bayes;
			if (gROOT && gROOT->FindObjectAny("unfolded_bayes"))delete gROOT->FindObjectAny("unfolded_bayes");
		}


		//Acceptance correction by the acceptance factor:

		cout << "applying the acc factor: ==================================="<<endl;
		TH2D* unfolded_bbb = (TH2D*)reco_test[k-1]->Clone();
		unfolded_bbb->Divide(acc_factor);

		int iter=15;
		fitOK=false;
		float tempChisq = 0.0;

		while(!fitOK)
		{    
			iter = iter+1;
			unfolded_bbb->Fit("fit2D_bbb","R0");
			tempChisq = fit2D_bbb->GetChisquare() / fit2D_bbb->GetNDF();
			fitOK = (tempChisq < 2.0 || iter >20 ) ? true : false;
		}   

		if(tempChisq < 2.0 ) {
			test_lambda_bbb->Fill(fit2D_bbb->GetParameter(1));
			test_mu_bbb->Fill(fit2D_bbb->GetParameter(2));
			test_nu_bbb->Fill(fit2D_bbb->GetParameter(3));
		}

		lambda_bbb->SetBinContent(4,fit2D_bbb->GetParameter(1));
		mu_bbb->SetBinContent(4,fit2D_bbb->GetParameter(2));
		nu_bbb->SetBinContent(4,fit2D_bbb->GetParameter(3));

		lambda_bbb->SetBinError(4,fit2D_bbb->GetParError(1));
		mu_bbb->SetBinError(4,fit2D_bbb->GetParError(2));
		nu_bbb->SetBinError(4,fit2D_bbb->GetParError(3));
		chi2NDF_bbb->SetBinContent(4, fit2D_bbb->GetChisquare() / fit2D_bbb->GetNDF());

		c1->cd(1);
		chi2NDF_bayes->Draw("HIST");
		chi2NDF_bbb->Draw("HIST same");
		chi2NDF_bbb->SetLineColor(kRed);

		c1->cd(2);
		lambda_bayes->SetLineColor(kBlue);
		lambda_bayes->SetMaximum(2.0);
		lambda_bayes->Draw("E");
		lambda_bayes->SetLineWidth(3);
		lambda_bayes->SetFillStyle(0);

		lambda_bbb->Draw("E2 same");
		lambda_bbb->SetLineColor(kRed);
		lambda_bbb->SetFillStyle(0);
		c1->cd(3);
		mu_bayes->SetLineColor(kBlue);
		mu_bayes->SetLineWidth(3);
		mu_bayes->SetFillStyle(0);
		mu_bayes->Draw("E");

		mu_bbb->Draw("E2 same");
		mu_bbb->SetLineColor(kRed);
		mu_bbb->SetFillStyle(0);

		c1->cd(4);
		nu_bayes->SetMarkerColor(kBlue);
		nu_bayes->SetLineColor(kBlue);
		nu_bayes->SetLineWidth(3);
		nu_bayes->Draw("E");
		nu_bayes->SetFillStyle(0);

		nu_bbb->Draw("E2 same");
		nu_bbb->SetLineColor(kRed);
		nu_bbb->SetFillStyle(0);
		c1->cd(5);
		unfolded_bbb->Draw("colz0");
		unfolded_bbb->SetTitle("Acceptance Factor Bin-By-Bin;#phi_{CS}; cos#theta_{CS}");
		c1->cd(6);
		RooUnfoldBayes unfoldBayes2(&response2D, reco_test[k-1],4);	
		TH2D* unfolded = (TH2D*) unfoldBayes2.Hunfold();
		unfolded->Fit("fit2D","R0");	

		tempChisq = fit2D->GetChisquare() / fit2D->GetNDF();
		/*
		   if(tempChisq < 2.0) { 
		   test_lambda_bayes->Fill(fit2D->GetParameter(1));
		   test_mu_bayes->Fill(fit2D->GetParameter(2));
		   test_nu_bayes->Fill(fit2D->GetParameter(3));
		   }   

		 */
		unfolded->Draw("colz");
		unfolded->SetTitle("Bayesian Iterative (with iteration =4);#phi_{CS}; cos#theta_{CS}");

		//c1->SaveAs(Form("plot/lmn_0.8_0.3_0.4/test_%i.png",k));
		//c1->SaveAs(Form("plot/lmn_0.6_0.2_0.2/test_%i.png",k));
		c1->SaveAs(Form("plot/test_%i.png",k));
		//c1->SaveAs(Form("plot/lmn_1_0_0/test_%i.png",k));
		//gStyle->SetOptStat(1111);

		delete lambda_bayes;
		delete mu_bayes;
		delete nu_bayes;
		delete lambda_bbb;
		delete mu_bbb;
		delete nu_bbb;
		delete unfolded;
		delete chi2NDF_bayes;
		delete chi2NDF_bbb;
		delete unfolded_bbb;
	}	

	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1111);

	TCanvas* c4 = new TCanvas("c4","c4",1800,600);
	c4->Divide(3,1);
	c4->cd(1);
	test_lambda_bayes->Draw("HIST");
	test_lambda_bayes->SetMaximum(35);
	test_lambda_bayes->SetLineColor(kBlue);
	test_lambda_bayes->SetLineWidth(2);
	test_lambda_bayes->SetTitle(";#lambda; Yield");

	test_lambda_bbb->Draw("HIST sames");
	test_lambda_bbb->SetLineColor(kRed);
	test_lambda_bbb->SetLineWidth(2);
/*
	TPaveStats *stats1 = (TPaveStats*)test_lambda_bayes->FindObject("stats");

	if (stats1 != NULL){
		stats1->SetY1NDC(.695+.200); stats1->SetY2NDC(.495+.200); stats1->SetTextSize(0.025); stats1->SetX1NDC(0.75); stats1->SetX2NDC(.995);
		stats1->SetTextColor(kBlue);
	}           


	TPaveStats *stats2 = (TPaveStats*)test_lambda_bbb->FindObject("stats");
	if (stats2 != NULL){
		stats2->SetY1NDC(.695); stats2->SetY2NDC(.495); stats2->SetTextSize(0.025); stats2->SetX1NDC(0.75); stats2->SetX2NDC(.995);
		stats2->SetTextColor(kRed); 
	}       
*/
	c4->cd(2);
	test_mu_bayes->Draw("HIST");
	test_mu_bayes->SetMaximum(35);
	test_mu_bayes->SetLineColor(kBlue);
	test_mu_bayes->SetLineWidth(2);
	test_mu_bayes->SetTitle(";#mu; Yield");

	test_mu_bbb->Draw("HIST sames");
	test_mu_bbb->SetLineColor(kRed);
	test_mu_bbb->SetLineWidth(2);
/*
	TPaveStats *stats3 = (TPaveStats*)test_mu_bayes->FindObject("stats");

	if (stats3 != NULL){
		stats3->SetY1NDC(.695+.200); stats3->SetY2NDC(.495+.200); stats3->SetTextSize(0.025); stats3->SetX1NDC(0.75); stats3->SetX2NDC(.995);
		stats3->SetTextColor(kBlue);
	}

	TPaveStats *stats4 = (TPaveStats*)test_mu_bbb->FindObject("stats");
	if (stats4 != NULL){
		stats4->SetY1NDC(.695); stats4->SetY2NDC(.495); stats4->SetTextSize(0.025); stats4->SetX1NDC(0.75); stats4->SetX2NDC(.995);
		stats4->SetTextColor(kRed);
	}

*/
	c4->cd(3);
	test_nu_bayes->Draw("HIST");
	test_nu_bayes->SetMaximum(35);
	test_nu_bayes->SetLineColor(kBlue);
	test_nu_bayes->SetLineWidth(2);
	test_nu_bayes->SetTitle(";#nu; Yield");

	test_nu_bbb->Draw("HIST sames");
	test_nu_bbb->SetLineColor(kRed);
	test_nu_bbb->SetLineWidth(2);
/*
	TPaveStats *stats5 = (TPaveStats*)test_nu_bayes->FindObject("stats");
	if (stats5 != NULL){
		stats5->SetY1NDC(.695+.200); stats5->SetY2NDC(.495+.200); stats5->SetTextSize(0.025); stats5->SetX1NDC(0.75); stats5->SetX2NDC(.995);
		stats5->SetTextColor(kBlue);
	}   
	TPaveStats *stats6 = (TPaveStats*)test_nu_bbb->FindObject("stats");
	if (stats6 != NULL){
		stats6->SetY1NDC(.695); stats6->SetY2NDC(.495); stats6->SetTextSize(0.025); stats6->SetX1NDC(0.75); stats6->SetX2NDC(.995);
		stats6->SetTextColor(kRed);
	}   
*/
	c4->SaveAs("plot/test.png");
	//c4->SaveAs("plot/lmn_0.8_0.3_0.4/test.png");
	//c4->SaveAs("plot/lmn_0.6_0.2_0.2/test.png");
	//c4->SaveAs("plot/lmn_1_0_0/test.png");
	//c4->SaveAs("plot/lmn_0.8_0.3_0.4/test.C");
	//c4->SaveAs("plot/lmn_0.6_0.2_0.2/test.C");
	//c4->SaveAs("plot/lmn_1_0_0/test.C");



}	
