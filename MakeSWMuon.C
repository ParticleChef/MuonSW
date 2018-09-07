#define MakeSWMuon_cxx
#include "MakeSWMuon.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>

void MakeSWMuon::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MakeSWMuon.C
//      root> MakeSWMuon t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

	TFile *file = new TFile("roi_median.root","recreate");

    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float low_et = 0;
   float high_et = 0;

//   int binSize = 90;

   vector<float> x[5][2], median[5][2], x_err[5][2], median_err[5][2];

   for( int nth = 0; nth < 189; nth++){
	   low_et = 10. + nth;
	   high_et = low_et + 1.;
	   cout << "et: " << low_et << endl;

	   vector<float> dPhi_l1me0_pix[5][2];
//	   float abs_pterr;
	   float closest_me0Et;
	   float closest_me0Eta;
	   float closest_me0Phi;
	   float closest_me0_dr;
   
   Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      //find closest me0 object to the gen
      int me0N_ = me0SegPosX->size(); // !!!!!! Use me0SegPosX instead me0 Pt !!!!!!
      if(me0N_ == 0 ) continue; //skip if there is no l1 me objects at all


      // closest me0
      float closest_dr = 9999;
      int closest_me0 = 0;
      int me0Count = 0;

      for( int i = 0; i < me0N_; i++){
	      TVector3 me0SegPos;
	      me0SegPos.SetXYZ(me0SegPosX->at(i), me0SegPosY->at(i), me0SegPosZ->at(i));
	      float me0Phi = me0SegPos.Phi();
	      float me0Eta = me0SegPos.Eta();

	      float dPhi = deltaPhi(propgenElPartPhi->at(0), me0Phi);
	      float current_dr = sqrt( pow(dPhi,2) + pow(propgenElPartEta->at(0) - me0Eta,2) );
	      if( genPartPt->at(0) > high_et || genPartPt->at(0) < low_et ) continue; // !!!!! Use genPartPt instead me0 Pt !!!!!
	      me0Count++;
	      if(current_dr < closest_dr){
		      closest_dr = current_dr;
		      closest_me0 = i;
	      }
      }// end of loop to find the closest me0 the gen


      if(closest_dr == 9999) continue;
      if(me0Count == 0) continue;

      cout<<"me0N_ : "<< me0N_ <<endl;
      cout<<"genPartPt = "<<genPartPt->at(0)<<endl;
      cout<<"closest_me0 = "<<closest_me0<<endl;
      cout<<"closest_dr = "<<closest_dr<<endl;
      TVector3 emvector;

//      abs_pterr = fabs(genPartPt->at(0)-propgenElPartPt->at(0))/genPartPt->at(0); // !!!!! Use propgenElPartPt instead closest_me0 Pt !!!!!
      closest_me0Et = propgenElPartPt->at(0); // !!!!! Use propgenElPartPt instead closest_me0 Pt !!!!!
      TVector3 me0SegPos;
      me0SegPos.SetXYZ(me0SegPosX->at(closest_me0), me0SegPosY->at(closest_me0), me0SegPosZ->at(closest_me0));
      closest_me0Eta = me0SegPos.Eta();
      closest_me0Phi = me0SegPos.Phi();
      closest_me0_dr = closest_dr;
      emvector.SetXYZ(me0SegPosX->at(closest_me0),me0SegPosY->at(closest_me0), me0SegPosZ->at(closest_me0));

      // store pixel hits into vectors
      std::vector<TVector3> first_disk_hits;
      std::vector<TVector3> second_disk_hits;
      std::vector<TVector3> third_disk_hits;
      std::vector<TVector3> fourth_disk_hits;
      std::vector<TVector3> fifth_disk_hits;
      std::vector<int> hitted_layers;

      int layers[6] = {}; // initialize as 0 for each event, save number of hits for each pixel layer
      layers[0] = 1; // beam spot

      int fpix_size = fRecHitGx->size();

      cout<<"fpix_size = "<<fpix_size<<endl;

      for( int a = 0; a < fpix_size; a++){

	      cout<<"fRecHitDisk = "<<fRecHitDisk->at(a)<<endl;

	      if( fRecHitDisk->at(a) == 1 ){
		      layers[1]++;
		      first_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 2 ){
		      layers[2]++;
		      second_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 3 ){
		      layers[3]++;
		      third_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 4 ){
		      layers[4]++;
		      fourth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
	      if( fRecHitDisk->at(a) == 5 ){
		      layers[5]++;
		      fifth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) ));
	      }
      }

      int n_pixels = 0; // initialize as 0 for each event
      hitted_layers.push_back(0); // 0 mean beam spot i.e., (0,0,0)
      for( int i = 1; i < 5; i++){
      	if( layers[i] != 0 ){
      		hitted_layers.push_back(i); // check if hits on each barrel or disk exists
      		n_pixels++;

      	}
      }

      cout<<"n_pixels = "<<n_pixels<<endl;

      if( n_pixels >= 1 ){
      	for( std::vector<int>::iterator first_hit = hitted_layers.begin()+1; first_hit != hitted_layers.end(); first_hit++){
		for ( int k = 0; k < layers[*first_hit]; k++){
			double dPhi = 0;
			double R = 0;

			TVector3 pixel_vector;

			if( *first_hit == 1 ) pixel_vector = first_disk_hits[k];
			if( *first_hit == 2 ) pixel_vector = second_disk_hits[k];
			if( *first_hit == 3 ) pixel_vector = third_disk_hits[k];
			if( *first_hit == 4 ) pixel_vector = fourth_disk_hits[k];
			if( *first_hit == 5 ) pixel_vector = fifth_disk_hits[k];

			dPhi = deltaPhi(pixel_vector.Phi(), closest_me0Phi);

			if( *first_hit == 1){
				if( fabs(closest_me0Eta) >= 2.0 && fabs(closest_me0Eta) < 2.4  ) dPhi_l1me0_pix[0][0].push_back(dPhi);
				if( fabs(closest_me0Eta) >= 2.4 && fabs(closest_me0Eta) < 2.8  ) dPhi_l1me0_pix[0][1].push_back(dPhi);
			}
			if( *first_hit == 2){
				if( fabs(closest_me0Eta) >= 2.0 && fabs(closest_me0Eta) < 2.4  ) dPhi_l1me0_pix[1][0].push_back(dPhi);
				if( fabs(closest_me0Eta) >= 2.4 && fabs(closest_me0Eta) < 2.8  ) dPhi_l1me0_pix[1][1].push_back(dPhi);
			}
			if( *first_hit == 3){
				if( fabs(closest_me0Eta) >= 2.0 && fabs(closest_me0Eta) < 2.4  ) dPhi_l1me0_pix[2][0].push_back(dPhi);
				if( fabs(closest_me0Eta) >= 2.4 && fabs(closest_me0Eta) < 2.8  ) dPhi_l1me0_pix[2][1].push_back(dPhi);
			}
			if( *first_hit == 4){
				if( fabs(closest_me0Eta) >= 2.0 && fabs(closest_me0Eta) < 2.4  ) dPhi_l1me0_pix[3][0].push_back(dPhi);
				if( fabs(closest_me0Eta) >= 2.4 && fabs(closest_me0Eta) < 2.8  ) dPhi_l1me0_pix[3][1].push_back(dPhi);
			}
			if( *first_hit == 5){
				if( fabs(closest_me0Eta) >= 2.0 && fabs(closest_me0Eta) < 2.4  ) dPhi_l1me0_pix[4][0].push_back(dPhi);
				if( fabs(closest_me0Eta) >= 2.4 && fabs(closest_me0Eta) < 2.8  ) dPhi_l1me0_pix[4][1].push_back(dPhi);
			}
		}
	}
      }
   }// event loop

   for( int pixel = 0; pixel < 5; pixel++){
   	for( int i = 0; i < 2; i++){
		std::sort (dPhi_l1me0_pix[pixel][i].begin(), dPhi_l1me0_pix[pixel][i].end());

		if( dPhi_l1me0_pix[pixel][i].size() != 0 ){
			x[pixel][i].push_back(low_et); median[pixel][i].push_back(getMedian(dPhi_l1me0_pix[pixel][i]));
			x_err[pixel][i].push_back(0.); median_err[pixel][i].push_back(getMedianErr(dPhi_l1me0_pix[pixel][i]));
		}
	}
   }
   }// et loop

   TGraphErrors* median_gr[5][2];
   for( int j = 0; j < 5; j++){
	   for( int i = 0; i < 2; i++){
		   int point_size = x[j][i].size();
		   median_gr[j][i] = new TGraphErrors(point_size, &x[j][i][0], &x_err[j][i][0], &median_err[j][i][0]);

		   median_gr[j][i]->SetMarkerStyle(24);
		   median_gr[j][i]->SetMarkerSize(0.3);
		   TString pixel_;
		   pixel_.Form("%d", j + 1);
		   TString eta_region;
		   eta_region.Form("%d", i + 1);
		   median_gr[j][i]->SetTitle(pixel_+"th_Pixel_"+eta_region+"_median");
		   median_gr[j][i]->SetName("Pixel_"+pixel_+"eta_region"+eta_region+"_median");

		   median_gr[j][i]->Write();
	   }
   }
   file->Write();
}
