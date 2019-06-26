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

   TH2F *dphi_dist = new TH2F("","",1000,0,100,1000,-0.1,0.1);
   TH2F *deta_dist = new TH2F("","",1000,0,100,1000,-0.1,0.1);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//   for (Long64_t jentry=0; jentry<1000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  if( me0SegNum == 0 || fRecHitN == 0 ) continue;
      if( genPartPt->at(1) < 10 ) continue;
	
	  float closest_dr = 9999;
	  int closest_me0 = -1;
	  int Nme0 = me0SegPosX->size();

	  for(int i = 0; i < Nme0; i++){
		  TVector3 me0SegPos;
		  me0SegPos.SetXYZ(me0SegPosX->at(i),me0SegPosY->at(i),me0SegPosZ->at(i));
		  float phi = me0SegPos.Phi();
		  float eta = me0SegPos.Eta();
		  float dphi = deltaPhi(genPartPhi->at(1),phi);

		  dphi_me0_D12_1->Fill(genPartPt->at(1),dphi);

		  float current_dr = sqrt(pow(dphi,2)+pow(genPartEta->at(1)-eta,2));
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

	  std::vector<TVector3> BS;
	  std::vector<TVector3> D1;
	  std::vector<TVector3> D2;
	  std::vector<TVector3> D3;
	  std::vector<TVector3> D4;
	  std::vector<TVector3> D5;

	  BS.push_back( TVector3(0,0,0));
	  for(int i = 0; i < fRecHitN; i++){
		  if( fRecHitDisk->at(i)== 1 ) D1.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 2 ) D2.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 3 ) D3.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 4 ) D4.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  if( fRecHitDisk->at(i)== 5 ) D5.push_back( TVector3(fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)) );
		  
	  }
	  for(int i = 0; i < D4.size(); i++){
			  float phi12 = ( D4.at(i)/**/- BS.at(0)/**/).Phi();
			  float Dphi = phi12 - phi;
			  if( Dphi > M_PI ) Dphi = Dphi - 2*M_PI;
			  if( Dphi <= -M_PI ) Dphi = Dphi + 2*M_PI;
			  if(region == 2 )dphi_dist->Fill(genPartPt->at(1),Dphi);
	  } 
   }// end of event

   TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);
   median_fitFunc->SetParLimits(2, -0.82, -0.80);
   median_fitFunc->SetParLimits(3, 0.17, 0.18);
   median_fitFunc->SetParLimits(4, 2.4, 2.5);

   TFile *file = new TFile("../roi_median.root","read");
//   TFile *file = new TFile("../ROI_offi.root","read");
   
   int nth_SW = 4;
   int Region = 2;

	TString nth_SW_;
	nth_SW_.Form("%d",nth_SW);
	
	TString Region_;
	Region_.Form("%d",Region);

   TGraphErrors* median_;

   TCanvas *c1 = new TCanvas();
   median_ = (TGraphErrors*)gDirectory->Get("Pixel_"+nth_SW_+"eta_region"+Region_+"_median");
   median_->GetYaxis()->SetRangeUser(-0.1, 0.1);
   median_->GetYaxis()->SetTitle("#Delta#phi");
   median_->GetXaxis()->SetRangeUser(0,100);
   median_->GetXaxis()->SetTitle("p_{T}(gen)");
//   median_->Draw("apX");


//   median_fitFunc->SetParLimits(2, -3, 3);
//   median_fitFunc->SetParLimits(4, -2, 2);
   median_fitFunc->SetLineColor(kRed);
   median_fitFunc->SetLineStyle(2);
   median_->Fit(median_fitFunc,"0W");

   dphi_dist->GetXaxis()->SetTitle("gen pt [GeV]");
   dphi_dist->GetYaxis()->SetTitle("#Delta#phi");
   dphi_dist->Draw("same scat=1");

   double ab[5] = {0};
   median_fitFunc->GetParameters(ab);

   int binSize = 100;

   double x1[90], y1[90], x2[90], y2[90];
   for(int j = 10; j<100; j++){
	   x1[j]=j;
	   x2[j]=0.;
	   y1[j] = ab[0]*pow(j,0)+ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3])+ab[4]);
	   y2[j] = 0.01;
   }
   c1->SetGrid();
   TGraphErrors *gr2 = new TGraphErrors(binSize,x1,y1,x2,y2);
   gr2->SetLineColor(kRed);
   gr2->SetFillColorAlpha(kYellow,0.20);
   gr2->SetFillStyle(3001);
   auto axis=gr2->GetXaxis();
   axis->SetLimits(10,100);

   gr2->Draw("sameE3");

   median_fitFunc->Draw("lsame");
//   Draw_draw_dphi(median_fitFunc, dphi_dist[i][j], i, j);


   c1->SaveAs("offi_plots/ROI_"+nth_SW_+"_region_"+Region_+".png");

   c1->Close();

   TCanvas *c2 = new TCanvas();
   dphi_me0_D12_1->GetXaxis()->SetTitle("#Delta#phi (gen,me)");
   dphi_me0_D12_1->Draw();
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
