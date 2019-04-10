void draw_median_fit(const std::string& median){

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  c1->SetRightMargin(0.1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetGridy();
  c1->SetGridx();
  c1->cd();

  cout << "string " << median << endl; 
  TGraphErrors* median_;
  median_ = (TGraphErrors*)gDirectory->Get((median).c_str() );

  median_->GetYaxis()->SetRangeUser(-0.1, 0.05);
  median_->SetTitle("");
  median_->Draw("apX");

}

void doFits(const std::string& median){

 TGraph* median_;
 median_ = (TGraph*)gDirectory->Get((median).c_str());

 TString nth_sw = (median).c_str();

 TF1 *median_fitFunc = new TF1("funcUp","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 200.);
 median_fitFunc->SetLineColor(kRed);
 median_fitFunc->SetLineStyle(2);

 median_->Fit(median_fitFunc,"0W");

 double ab[5] = {0};
 median_fitFunc->GetParameters(ab);
// cout<<ab[0]*pow(2,0)+ab[1]*pow(2,ab[2])*exp(-pow(2,ab[3])+ab[4]);

 int binSize = 100;

 double x1[binSize], y1[binSize], x2[binSize], y2[binSize];
 for(int j = 2; j<binSize; j++){
	 x1[j]=j;
	 x2[j]=0.;
	 y1[j] = ab[0]*pow(j,0)+ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3])+ab[4]);
	 y2[j] = 0.01;
 }
 TGraphErrors *gr2 = new TGraphErrors(binSize,x1,y1,x2,y2);
 gr2->SetLineColor(kBlue);
 gr2->SetFillColorAlpha(kBlue,0.20);
 gr2->SetFillStyle(3001);

 gr2->Draw("sameE3");
 median_fitFunc->Draw("lsame");

 TCanvas* Current = gPad->GetCanvas();
// Current->SaveAs(nth_sw+"_roi.pdf");
// Current->Print("ROI_offi/"+nth_sw+"_roi.png"); //save in directory
 Current->Print(nth_sw+"_roi.png");

 delete Current;

 // save the fit parameters
 ofstream fit_result;
 char fit_parameter[50];
 sprintf(fit_parameter, "./roi_offi_.txt");

 fit_result.open(fit_parameter);

 fit_result << endl;
 fit_result << "if( dp_de_dr == 0 && i == 6  && up_down == 0){" <<endl;
 for( int i=0; i < 5; i++){
     fit_result << "p[" << i << "] = " << median_fitFunc->GetParameter(i) << ";" << endl;
 }
 //fit_result << "return p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);" << endl;
 fit_result << "}" << endl;
 fit_result << endl;

 fit_result.close();

}
