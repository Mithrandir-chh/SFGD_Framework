// ROOT libraries
#include "TROOT.h"
#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObject.h"
#include <TStyle.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TClass.h>
#include <TChain.h>
#include <TLegend.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Fit/Fitter.h>
#include <TRandom3.h>
#include <TEfficiency.h>

//C++ libraries
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <array>
#include <time.h>
#include <map>

#include "TSystem.h"
#include "TMacro.h"
#include <sstream>
#include <iterator>
#include <algorithm>
#include "TImage.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "Rtypes.h"
#include <TPaveLabel.h>

void auto_check() {

  int nBoard = 18;
  TFile f(Form("%s", gApplication->Argv(1)));
  TTree* tt[nBoard];
  Int_t FEBs[nBoard] = {0,1,2,3,4,8,9,10,11,12,16,17,18,19,20,24,25,26,27};
  tt[0] = (TTree*)f.Get("FEB_0");
  tt[1] = (TTree*)f.Get("FEB_1");
  tt[2] = (TTree*)f.Get("FEB_2");
  tt[3] = (TTree*)f.Get("FEB_3");
  tt[4] = (TTree*)f.Get("FEB_4");
  tt[5] = (TTree*)f.Get("FEB_8");
  tt[6] = (TTree*)f.Get("FEB_9");
  tt[7] = (TTree*)f.Get("FEB_10"); 
  tt[8] = (TTree*)f.Get("FEB_11");
  tt[9] = (TTree*)f.Get("FEB_12");
  tt[10] = (TTree*)f.Get("FEB_16");
  tt[11] = (TTree*)f.Get("FEB_17");
  tt[12] = (TTree*)f.Get("FEB_18");
  tt[13] = (TTree*)f.Get("FEB_19");
  tt[14] = (TTree*)f.Get("FEB_20");
  tt[15] = (TTree*)f.Get("FEB_24");
  tt[16] = (TTree*)f.Get("FEB_25");
  tt[17] = (TTree*)f.Get("FEB_26");
  tt[18] = (TTree*)f.Get("FEB_27");
  TH2F* h2[30];
  TH1F* h1[30];
  for (Int_t i=0;i<30;i++){
    h2[i] = new TH2F("","Channel vs HG Counts; Channel number; HG Counts",400,0,2000,96,0,96);
    h1[i] = new TH1F("","Hit time from Spill; Hit time from spill; Counts",100,0,4.E8);
  }

  TCanvas *cc = new TCanvas("cc", "cc", 800, 800);
  cc->Divide(5,6);
  {
    std::cout<<"plotting .. "<<std::endl;
    cc->cd(1);
    tt[0]->Draw("FEB_0_hitsChannel:FEB_0_hitAmpl>> h2[0]","FEB_0_hitsChannel>0 && FEB_0_hitAmpl>0","colz");
    cc->cd(2);
    tt[1]->Draw("FEB_1_hitsChannel:FEB_1_hitAmpl>> h2[1]","FEB_1_hitsChannel>0 && FEB_1_hitAmpl>0","colz");
    cc->cd(3);
    tt[2]->Draw("FEB_2_hitsChannel:FEB_2_hitAmpl>> h2[2]","FEB_2_hitsChannel>0 && FEB_2_hitAmpl>0","colz");
    cc->cd(4);
    tt[3]->Draw("FEB_3_hitsChannel:FEB_3_hitAmpl>> h2[3]","FEB_3_hitsChannel>0 && FEB_3_hitAmpl>0","colz");
    cc->cd(5);
    tt[4]->Draw("FEB_4_hitsChannel:FEB_4_hitAmpl>> h2[4]","FEB_4_hitsChannel>0 && FEB_4_hitAmpl>0","colz");
    cc->cd(6);
    tt[5]->Draw("FEB_8_hitsChannel:FEB_8_hitAmpl>> h2[5]","FEB_8_hitsChannel>0 && FEB_8_hitAmpl>0","colz");
    cc->cd(7);
    tt[6]->Draw("FEB_9_hitsChannel:FEB_9_hitAmpl>> h2[6]","FEB_9_hitsChannel>0 && FEB_9_hitAmpl>0","colz");
    cc->cd(8);
    tt[7]->Draw("FEB_10_hitsChannel:FEB_10_hitAmpl>> h2[7]","FEB_10_hitsChannel>0 && FEB_10_hitAmpl>0","colz");
    cc->cd(9);
    tt[8]->Draw("FEB_11_hitsChannel:FEB_11_hitAmpl>> h2[8]","FEB_11_hitsChannel>0 && FEB_11_hitAmpl>0","colz");
    cc->cd(10);
    tt[9]->Draw("FEB_12_hitsChannel:FEB_12_hitAmpl>> h2[9]","FEB_12_hitsChannel>0 && FEB_12_hitAmpl>0","colz");
    cc->cd(11);
    tt[10]->Draw("FEB_16_hitsChannel:FEB_16_hitAmpl>> h2[10]","FEB_16_hitsChannel>0 && FEB_16_hitAmpl>0","colz");
    cc->cd(12);
    tt[11]->Draw("FEB_17_hitsChannel:FEB_17_hitAmpl>> h2[11]","FEB_17_hitsChannel>0 && FEB_17_hitAmpl>0","colz");
    cc->cd(13);
    tt[12]->Draw("FEB_18_hitsChannel:FEB_18_hitAmpl>> h2[12]","FEB_18_hitsChannel>0 && FEB_18_hitAmpl>0","colz");
    cc->cd(14);
    tt[13]->Draw("FEB_19_hitsChannel:FEB_19_hitAmpl>> h2[13]","FEB_19_hitsChannel>0 && FEB_19_hitAmpl>0","colz");
    cc->cd(15);
    tt[14]->Draw("FEB_20_hitsChannel:FEB_20_hitAmpl>> h2[14]","FEB_20_hitsChannel>0 && FEB_20_hitAmpl>0","colz");
    cc->cd(16);
    tt[15]->Draw("FEB_24_hitsChannel:FEB_24_hitAmpl>> h2[15]","FEB_24_hitsChannel>0 && FEB_24_hitAmpl>0","colz");    
    cc->cd(17);
    tt[16]->Draw("FEB_25_hitsChannel:FEB_25_hitAmpl>> h2[16]","FEB_25_hitsChannel>0 && FEB_25_hitAmpl>0","colz");
    cc->cd(18);
    tt[17]->Draw("FEB_26_hitsChannel:FEB_26_hitAmpl>> h2[17]","FEB_26_hitsChannel>0 && FEB_26_hitAmpl>0","colz");
    cc->cd(19);
    tt[18]->Draw("FEB_27_hitsChannel:FEB_27_hitAmpl>> h2[18]","FEB_27_hitsChannel>0 && FEB_27_hitAmpl>0","colz");

  }
  cc->Update();
  cc->WaitPrimitive();
  for(int i=0;i<nBoard;i++){
    tt[i]->Reset();
  }  
  f.Close();  
  exit(0);
}

int main(int argc, char* argv[]) {
  TApplication theApp("woo",&argc, argv);
  auto_check();
  theApp.Run();
  exit(0);
  return 0;
}
