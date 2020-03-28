#define plot_cxx
#include "plot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMath.h>

//void Draw_draw_dphi(TF1* tf, TH2F* dphi_dist[6][2], int nth_SW, int Region);   
	
void plot::Loop()
{
	gStyle->SetOptStat(0);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH2F *dphi_me0_D12_1 = new TH2F("","",100,0,100,1000,-0.3,0.3);
   TH2F *deta_me0_D12_1 = new TH2F("","",100,0,100,1000,-0.3,0.3);

   TH2F *dphi_dist = new TH2F("","",1000,0,100,1000,-0.2,0.2);
   TH2F *deta_dist = new TH2F("","",1000,0,100,1000,-0.2,0.2);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//   for (Long64_t jentry=0; jentry<10;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  if( me0SegNum == 0 || fRecHitN == 0 ) continue;

	  float closest_dr = 9999;
	  int closest_me0 = -1;
	  int Nme0 = me0SegPosX->size();

	  for(int i = 0; i < Nme0; i++){
		  TVector3 me0SegPos;
		  me0SegPos.SetXYZ(me0SegPosX->at(i),me0SegPosY->at(i),me0SegPosZ->at(i));
		  float phi = me0SegPos.Phi();
		  float eta = me0SegPos.Eta();
		  float dphi = deltaPhi(genPartPhi->at(0),phi);

		  float current_dr = sqrt(pow(dphi,2)+pow(genPartEta->at(0)-eta,2));
		  //if( genPartPt->at(0) < 20 ) continue;
		  if( current_dr < closest_dr ){
			  closest_dr = current_dr;
			  closest_me0 = i;
		  }
	  }
	  if(closest_dr == 9999. || closest_me0 == -1) continue;

	  int region = 0;

	  TVector3 me0SegPos;
	  me0SegPos.SetXYZ(me0SegPosX->at(closest_me0),me0SegPosY->at(closest_me0),me0SegPosZ->at(closest_me0));
	  float eta = me0SegPos.Eta();
	  float phi = me0SegPos.Phi();
	  if(fabs(eta) > 2.0 && fabs(eta) <= 2.4 ) region = 1;
	  if(fabs(eta) > 2.4 && fabs(eta) < 2.8 )  region = 2;
	  if(fabs(eta) < 2.0 || fabs(eta) > 2.8 ) continue;

	  if( region == 0 ) continue;

	  std::vector<TVector3> D1;
	  std::vector<TVector3> D2;
	  std::vector<TVector3> D3;
	  std::vector<TVector3> D4;
	  std::vector<TVector3> D5;

	  for(int i = 0; i < fRecHitN; i++){
		  if( fRecHitDisk->at(i)== 1 ) D1.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 2 ) D2.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 3 ) D3.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 4 ) D4.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 5 ) D5.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  
	  }
	  for(int i = 0; i < D2.size(); i++){
		  for(int j = 0; j < D3.size(); j++){
			  float phi12 = (D3.at(j) - D2.at(i)).Phi();
			  float eta12 = (D3.at(j) - D2.at(i)).Eta();
			  float dphi = phi12 - phi;
			  if( dphi > M_PI ) dphi = dphi - 2*M_PI;
			  if( dphi <= -M_PI ) dphi = dphi + 2*M_PI;
			  if(region == 2 && genPartPt->at(0) > 10 )dphi_dist->Fill(genPartPt->at(0),dphi);
			  if(region == 2 && genPartPt->at(0) > 10 )deta_dist->Fill(genPartPt->at(0),eta - eta12);
		  }
	  } 
   }// end of event

   TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);
/////// Setting the parameter limits //////////
//   median_fitFunc->SetParLimits(0,-0.0013,0);
//   median_fitFunc->SetParLimits(1,0.3,0.32);
//   median_fitFunc->SetParLimits(2,-1,0.0);
//   median_fitFunc->SetParLimits(3,-0.02,-0.01);
//   median_fitFunc->SetParLimits(4,1.0,1.05);

    median_fitFunc->SetParameter(0,0.00106159);
	median_fitFunc->SetParameter(1,0.00288571);
	median_fitFunc->SetParameter(2,0.);
	median_fitFunc->SetParameter(3,0.);
	median_fitFunc->SetParameter(4,0.);
//   median_fitFunc->SetParLimits(0,0.0025,0.0026);
//   median_fitFunc->SetParLimits(1,0,0.5);
//   median_fitFunc->SetParLimits(2,-0.5,0.1);
//   median_fitFunc->SetParLimits(3,0.2,0.4);
//   median_fitFunc->SetParLimits(4,0.4,0.6);

   //TFile *file = new TFile("../ME0_median.root","read");
   TFile *file = new TFile("../roi_median.root","read");
   //file->GetXaxis()->SetRangeUser(0,100);
   
   int nth_SW = 1;
   int Region = 2;
   int eta_or_phi = 1; //eta = 1, phi = 2

	TString nth_SW_;
	nth_SW_.Form("%d",nth_SW);
	
	TString Region_;
	Region_.Form("%d",Region);

   TGraphErrors* median_;

   TCanvas *c1 = new TCanvas();
   if( eta_or_phi ==1 ) median_ = (TGraphErrors*)gDirectory->Get("ME0_Pix_deta_SW_"+nth_SW_+"_eta_region_"+Region_+"_median");
   if( eta_or_phi ==2 ) median_ = (TGraphErrors*)gDirectory->Get("ME0_Pix_dphi_SW_"+nth_SW_+"_eta_region_"+Region_+"_median");
   median_->GetYaxis()->SetRangeUser(-0.2, 0.2);
   median_->GetXaxis()->SetRangeUser(0, 100);
   median_->SetTitle("");
   median_->GetXaxis()->SetTitle("p_{T(gen)}");
   if(eta_or_phi == 1) median_->GetYaxis()->SetTitle("#Delta#eta");
   if(eta_or_phi == 2) median_->GetYaxis()->SetTitle("#Delta#phi");
   median_->Draw("apX");


   median_fitFunc->SetLineColor(kRed);
   median_fitFunc->SetLineStyle(2);
   median_fitFunc->SetLineWidth(3);
   median_->Fit(median_fitFunc,"0W");

   if(eta_or_phi == 1){
   //deta_dist->GetYaxis()->SetTitle("#Delta#phi");
   deta_dist->Draw("same scat=1");
   }

   if(eta_or_phi == 2){
   //dphi_dist->GetYaxis()->SetTitle("#Delta#phi");
   dphi_dist->Draw("same scat=1");
   }

   double ab[5] = {0};
   median_fitFunc->GetParameters(ab);

   int binSize = 91;

   double x1[binSize], y1[binSize], x2[binSize], y2[binSize];
   for(int j = 0; j<binSize; j++){
	   x1[j]=j+10;
	   x2[j]=10.;
	   y1[j] = ab[0]*pow(j+10,0)+ab[1]*pow(j+10,ab[2])*exp(-pow(j+10,ab[3])+ab[4]);
	   if(eta_or_phi == 1 ) y2[j] = 0.05;
	   if(eta_or_phi == 2 ) y2[j] = 0.02;
   }
  // if(eta_or_phi == 1) y1[0] = 0.;
  // if(eta_or_phi == 2) y1[0] = 0.02;
   for(int i = 0; i < 20; i++){
	   cout<<"y1["<<i<<"] = "<<y1[i]<<endl;
   }
   c1->SetGrid();
   TGraphErrors *gr2 = new TGraphErrors(binSize,x1,y1,x2,y2);
   gr2->GetXaxis()->SetRangeUser(10,100);
   gr2->GetXaxis()->SetLimits(10,100);
//   auto axis = gr2->GetXaxis();
//   axis->SetLimits(10.,200);
   gr2->SetLineColor(kRed);
   gr2->SetLineWidth(3);
   gr2->SetFillColorAlpha(kYellow,0.20);
   gr2->SetFillStyle(3001);

   gr2->Draw("sameE3");

   median_fitFunc->Draw("lsame");


   if(eta_or_phi == 1) c1->SaveAs("final_plot/deta_me0_pix_sw_"+nth_SW_+"_region_"+Region_+".png");
   if(eta_or_phi == 2) c1->SaveAs("final_plot/dphi_me0_pix_sw_"+nth_SW_+"_region_"+Region_+".png");

//   c1->Close();

}
/*
void Draw_draw_dphi(TF1* tf, TH2F* dphi_dist[6][2], int nth_SW, int Region){   
	TString nth_SW_;
	nth_SW_.Form("%d",nth_SW);
	
	TString Region_;
	Region_.Form("%d",Region);

	TH2F *hist;
	hist = (TH2F*)dphi_dist[nth_SW][Region]->Clone(); 

   tf->SetLineColor(kRed);
   tf->SetLineStyle(2);
   median_->Fit(tf,"0W");

   double ab[5] = {0};
   tf->GetParameters(ab);

   int binSize = 200;

   double x1[200], y1[200], x2[200], y2[200];
   for(int j = 2; j<200; j++){
	   x1[j]=j;
	   x2[j]=0.;
	   y1[j] = ab[0]*pow(j,0)+ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3])+ab[4]);
	   y2[j] = 0.03;
   }
   TGraphErrors *gr2 = new TGraphErrors(binSize,x1,y1,x2,y2);
   gr2->SetLineColor(kRed);
   gr2->SetFillColorAlpha(kBlue,0.20);
   gr2->SetFillStyle(3001);

   gr2->Draw("sameE3");

   tf->Draw("lsame");
   
   return 0;
}
*/
