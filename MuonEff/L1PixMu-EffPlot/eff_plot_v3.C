#define eff_plot_cxx
#include "eff_plot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

void eff_plot::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   int nbins = 40; float x1 = 2.0; float x2 = 2.8;
   
   TH1F* hEG_denom = new TH1F("hEG_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hEG_nom = new TH1F("hEG_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hEG_nom->Sumw2();
   hEG_denom->Sumw2();

   TH1F* hPix_denom = new TH1F("hPix_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hPix_nom = new TH1F("hPix_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hPix_nom->Sumw2();
   hPix_denom->Sumw2();

   TH1F* hPix4_denom = new TH1F("hPix4_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hPix4_nom = new TH1F("hPix4_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hPix4_nom->Sumw2();
   hPix4_denom->Sumw2();

   TH1F* hTrk_denom = new TH1F("hTrk_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hTrk_nom = new TH1F("hTrk_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hTrk_nom->Sumw2();
   hTrk_denom->Sumw2();

   TH1F *me0_phi = new TH1F("","",100,-4,4);
   TH1F *me0_eta = new TH1F("","",100,1.9,3);

   TH1F *me0_dphi = new TH1F("","",500,-0.2,0.2);
   TH1F *me0_deta = new TH1F("","",500,-0.2,0.2);
   TH1F *me0_dr   = new TH1F("","",500,-0.2,0.2);

   TH1F *cl_me0_dphi = new TH1F("","",500,-0.2,0.2);
   TH1F *cl_me0_deta = new TH1F("","",500,-0.2,0.2);
   TH1F *cl_me0_dr   = new TH1F("","",500,-0.2,0.2);

   int nbins_pt = 45; float x1_pt = 10.; float x2_pt = 100.;   
//   int nbins_pt = 20; float x1_pt = 10.; float x2_pt = 50.;   
   TH1F* hEG_pt = new TH1F("hEG_pt",";P_{T} of generated electron; Efficiency",nbins_pt,x1_pt,x2_pt);
   TH1F* hPix_pt = new TH1F("hPix_pt",";P_{T} of generated electron [GeV]; Efficiency",nbins_pt,x1_pt,x2_pt);

   int width_bit[27] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000, 0x10000, 0x20000, 0x40000, 0x80000,0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000};

   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  cout<<"==== "<<jentry<<" ===="<<endl;
	  cout<<"gen Eta = "<<nt_genEta<<endl;

      if(nt_genPt < 10.) continue;
	  if( fabs(nt_genEta) < 2.0 || fabs(nt_genEta) > 2.8 ) continue;
//      if( fabs(nt_genEta) > 3.0 ) continue;

	  cout<<"after cut"<<endl;
	  cout<<"gen Eta = "<<nt_genEta<<endl;

      int eg_size = ntME0Et->size();
      int gen_matched_eg = -1;
      float dr = 999.;

	  for(int i = 0; i < eg_size; i++){
		  float deta = nt_genEta - ntME0Eta->at(i);
		  float dphi = nt_genPhi - ntME0Phi->at(i);
		  if(dphi > M_PI) dphi = dphi - 2*M_PI;
		  if(dphi <= -M_PI) dphi = dphi + 2*M_PI;

		  float dR = sqrt(pow(dphi,2)+pow(deta,2));

		  me0_dphi->Fill(dphi);
		  me0_deta->Fill(deta);
		  me0_dr  ->Fill(dR);
	  }
      // find the closest L1 egamma to gen electron
      for(int i = 0; i < eg_size; i++){
		  float dphi = nt_genPhi - ntME0Phi->at(i);
		  if(dphi > M_PI) dphi = dphi - 2*M_PI;
		  if(dphi <= -M_PI) dphi = dphi + 2*M_PI;
         float temp_dr = sqrt(pow(dphi,2)+pow(nt_genEta-ntME0Eta->at(i),2));
         if( temp_dr < dr){
           dr = temp_dr;
           gen_matched_eg = i;
         }
      }// eg loop
     
      float pt_err = 0.; 
      if(gen_matched_eg != -1) {pt_err = fabs(nt_genPt-ntME0Et->at(gen_matched_eg))/nt_genPt; }
      if( gen_matched_eg != -1 /* && pt_err < 0.5 */){

		  if((trigger_bit_width->at(gen_matched_eg)&width_bit[0])==width_bit[0]){
			  hPix_nom->Fill(nt_genEta, 1.);
			  float deta = nt_genEta - ntME0Eta->at(gen_matched_eg);
			  float dphi = nt_genPhi - ntME0Phi->at(gen_matched_eg);
			  if(dphi > M_PI) dphi = dphi - 2*M_PI;
			  if(dphi <= -M_PI) dphi = dphi + 2*M_PI;

			  float dR = sqrt(pow(dphi,2)+pow(deta,2));
			  cl_me0_dphi->Fill(dphi);
			  cl_me0_deta->Fill(deta);
			  cl_me0_dr->Fill(deta);

			  me0_phi->Fill(ntME0Phi->at(gen_matched_eg));
			  me0_eta->Fill(ntME0Eta->at(gen_matched_eg));
		  }

        hEG_nom->Fill(nt_genEta, 1.);
      }

      hEG_denom->Fill(nt_genEta, 1.);
      hPix_denom->Fill(nt_genEta, 1.);

   }// events loop

   // divide the two histograms to plot efficiency 
   TGraphAsymmErrors* hEG = new TGraphAsymmErrors(hEG_nom, hEG_denom,"B");
   TGraphAsymmErrors* hPix = new TGraphAsymmErrors(hPix_nom, hPix_denom,"B");

/****************************************************************/
   for(int bin = 0; bin < nbins; bin++){
	   float nom_event = hPix_nom->GetBinContent(bin);
	   float den_event = hPix_denom->GetBinContent(bin);
	   cout<<bin<<"th bin (Eta "<<2.0+0.02*bin<<" ): nom / denom = "<<nom_event<<" / "<<den_event<<endl;
   }
/*
   double sum=0;
   double x, xx;
   double y, yy;
   for(int i = 0; i < 20; i++){
	   hPix_nom->GetPoint(i,x,y);
	   hPix_denom->GetPoint(i,xx,yy);
	   cout<<"bin : "<<i<<", eta = "<<x<<endl;
	   cout<<"nominator: "<<y<<"denominator: "<<yy<<endl;
	  // sum += y;
   }
   cout<<"ave = "<<sum/20<<endl;
*/
 //  TCanvas *cc = new TCanvas();
 //  hPix_denom->SetFillColor(kRed);
 //  hPix_denom->Draw();
 //  hPix_nom->SetFillColor(kBlue);
 //  hPix_nom->Draw("same");


   TCanvas *c1 = new TCanvas("c1","c1",1200,900);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(1); // axis width, default is 1
   c1->SetTopMargin(0.1);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.15);
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);
   c1->cd();

   //hPix->SetTitle("CMSSW_10_1_0_pre3, Phase 2, <PU>=0");
   hPix->GetXaxis()->SetTitle("#eta_{gen} ");
   hPix->GetXaxis()->SetTitleOffset(1.);
   hPix->GetXaxis()->SetTitleSize(0.055);
   //hPix->GetXaxis()->SetNdivisions(510);
   //hPix->GetYaxis()->SetNdivisions(512);
   hPix->GetXaxis()->SetLabelSize(0.05);
   hPix->GetYaxis()->SetLabelSize(0.05);
   hPix->GetXaxis()->SetRangeUser(1.9, 3.0);
   hPix->GetYaxis()->SetRangeUser(0.0, 1.05);
   hPix->GetYaxis()->SetTitle("Efficiency");
   hPix->GetYaxis()->SetTitleOffset(1.2);
   hPix->GetYaxis()->SetTitleSize(0.055);

   hPix->SetMarkerColor(kBlue);
   hPix->SetLineColor(kBlue);
   hPix->SetLineWidth(1);
   hPix->SetMarkerStyle(29);
   hPix->SetMarkerSize(1.8);
   hPix->Draw("ape");


   hEG->SetMarkerColor(kRed);
   hEG->SetLineColor(kRed);
   hEG->SetLineWidth(1);
   hEG->SetMarkerStyle(20);
   hEG->SetMarkerSize(1.);
   hEG->Draw("pe same");

   TLegend *Lgd = new TLegend(0.25, 0.9, 0.90, 0.95);

   Lgd-> SetNColumns(4);
   Lgd->SetFillColor(0);
   Lgd->SetTextFont(42);
   Lgd->SetTextSize(0.035);
   Lgd->SetBorderSize(0);
   Lgd->SetFillStyle(0);
   //Lgd->AddEntry(hEG,"Phase-2 L1 EG(barrel ECAL/ HGCAL)","lp");
   Lgd->AddEntry(hEG,"Phase-2 L1 ME0","lp");
   Lgd->AddEntry(hPix,"Pixel matching","lp");
   //Lgd->AddEntry(hTrk,"L1 track matching","lp");
   //Lgd->AddEntry(hPix4,"Pixel matching(4 hits)","lp");
   Lgd->Draw();

   //TLatex t(0.1,1.28,"CMSSW_10_1_0_pre3, Phase 2, <PU>=200");
   TLatex t(0.,1.11,"CMSSW_10_1_5, Phase 2, <PU>=200");
   t.SetTextSize(0.035);
   t.Draw();

  // TLatex pt_cut(2.,0.2,"p_{T}^{gen} > 10 GeV");
  // pt_cut.SetTextSize(0.035);
  // pt_cut.Draw();

   TCanvas *c3 = new TCanvas();
   me0_phi->GetXaxis()->SetTitle("#phi");
   me0_phi->GetYaxis()->SetTitle("# of closest me0 segment");
   me0_phi->Draw();

   TCanvas *c4 = new TCanvas();
   me0_eta->GetXaxis()->SetTitle("#eta");
   me0_eta->GetYaxis()->SetTitle("# of closest me0 segment");
   me0_eta->Draw();

   TCanvas *c5 = new TCanvas();
   me0_deta->GetXaxis()->SetTitle("#Delta #eta");
   me0_deta->GetYaxis()->SetTitle("# of me0 segment");
   me0_deta->Draw();

   TCanvas *c6 = new TCanvas();
   me0_dphi->GetXaxis()->SetTitle("#Delta #phi");
   me0_dphi->GetYaxis()->SetTitle("# of me0 segment");
   me0_dphi->Draw();

   TCanvas *c7 = new TCanvas();
   me0_dr->GetXaxis()->SetTitle("#Delta R");
   me0_dr->GetYaxis()->SetTitle("# of me0 segment");
   me0_dr->Draw();

   TCanvas *c8 = new TCanvas();
   cl_me0_deta->GetXaxis()->SetTitle("#Delta #eta");
   cl_me0_deta->GetYaxis()->SetTitle("# of closest me0 segment");
   cl_me0_deta->Draw();

   TCanvas *c9 = new TCanvas();
   cl_me0_dphi->GetXaxis()->SetTitle("#Delta #phi");
   cl_me0_dphi->GetYaxis()->SetTitle("# of closest me0 segment");
   cl_me0_dphi->Draw();

   TCanvas *c10 = new TCanvas();
   cl_me0_dr->GetXaxis()->SetTitle("#Delta R");
   cl_me0_dr->GetYaxis()->SetTitle("# of closest me0 segment");
   cl_me0_dr->Draw();

   c1->Print("test_Eff_200PU.png");
 //  cc->SaveAs("nom_denom.png");
}
