{

  TFile f(Form("%s",gApplication->Argv(4) ));

  TTree* t = (TTree*)f.Get("FEB_1");

  std::vector<double>* FEB_1_hitsChannel;
  std::vector<double>* FEB_1_hitAmpl;

  t->SetBranchAddress("FEB_1_hitsChannel",&FEB_1_hitsChannel);
  t->SetBranchAddress("FEB_1_hitAmpl",&FEB_1_hitAmpl);

  TH1F* h[100];
  for(int i=0;i<100;i++){
    h[i] = new TH1F("","Gain; ADC; Counts",500,0,1000);
  
  }

  for(int i=0;i<t->GetEntries()-1; i++){

    t->GetEntry(i);

    for(int j=0;j< FEB_1_hitsChannel->size()/2;j++){
      int chn = (*FEB_1_hitsChannel)[j];
      double ampl = (*FEB_1_hitAmpl)[j];
      if(ampl > 0)
        h[chn] ->Fill(ampl);
      //cout<<chn<<" "<<ampl<<endl;
    }
  }

  TCanvas* c1 = new TCanvas();
  c1->Divide(6,6);
  for(int i=0;i<32;i++){
    c1->cd(i+1);
    h[i] -> Draw();
  }
  TCanvas* c2 = new TCanvas();
  c2->Divide(6,6);
  for(int i=0;i<32;i++){
    c2->cd(i+1);
    h[i+32] -> Draw();
  }
  TCanvas* c3 = new TCanvas();
  c3->Divide(6,6);
  for(int i=0;i<32;i++){
    c3->cd(i+1);
    h[i+64] -> Draw();
  }

  new TCanvas();
  h[64]->Draw();
  new TCanvas();
  h[65]->Draw();

}
