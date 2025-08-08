#define THIS_NAME nEventDisplay
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS

#include "../src/tools/global_header.hh"
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

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool    batch          = false;
bool    IsMC           = true;
bool    IsReversed     = false; 
bool    IsTwenty       = false;
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "~/Desktop/templateAna_Out.root";
TString fileIn         = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
int IsPrototype        = true;
int SFGD_X             = 204;
int SFGD_Y             = 56;
int SFGD_Z             = 188;
bool useClean          = false;
bool useNN             = false;
int dTimeIni           = -100;
int dTimeFin           = -100;
int evtIni             = 0;
int evtFin             = 0;

TFile*          FileInput;
TFile*          FileOutput;
TTree*          dataOut;
TTree*          data;
ND280SFGDEvent* inputEvent;
ND280SFGDEvent* recoEvent;
TBranch*        recoBranch; 
TBranch*        inputBranch;
int             nEvents;
Event*          unpackEvent;

#include "../src/tools/global_tools.cc" 
#include "../src/tools/reconstruction_methods.hh" 

void nEventDisplay() {

  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent, 
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  std::vector<ND280SFGDHit*> mppc_hits;

  // Specify directory for output Event Display PNG
  //string foutDir = GetDir(fileIn.Data()) + "/EventDisplay_PNG/";
  //string foutDir2 = GetDir(fileIn.Data()) + "/hitTimes_PNG/";
  //string createFolder = "mkdir -p ";
  //string createFolder2 = "mkdir -p ";
  //createFolder += foutDir.c_str();
  //createFolder2 += foutDir2.c_str();
  //system(createFolder.c_str());
  //system(createFolder2.c_str());

  //ostringstream foutPNGnum;
  //string foutPNG;

  cout<<"total events = "<<nEvents<<endl;
  cout<<"displaying "<<evtFin-evtIni<<" events"<<endl;

  // Filling Histograms_______________________________________________________________________________

  TH2F *hXY = new TH2F("XY view", "XY view; X position [cm];Y position [cm]; Hit time from spill [2.5 ns]", 24, 0, 24, 8, 0, 8);
  TH2F *hXZ = new TH2F("XZ view", "XZ view; X position [cm];Z position [cm]; Hit time from spill [2.5 ns]", 24, 0, 24, 48, 0, 48);
  TH2F *hZY = new TH2F("ZY view", "ZY view; Z position [cm];Y position [cm]; Hit time from spill [2.5 ns]", 48, 0, 48, 8, 0, 8);

  Int_t nNeutrons = 0;
  Int_t QmaxX = 0;
  Int_t QmaxY = 0;
  Int_t QmaxZ = 0;
  Int_t QmaxNum = 0;
  Int_t eventsMissed = 0;
  Int_t prevEventTime = 0;
  gStyle->SetOptStat(0);
  //loop over events
  for (int iev=evtIni; iev<evtFin; iev++){
    if (iev == maxEvents-1) break;

    //cout<<"_Plotting events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;

    mppc_hits = getEventMPPChits(iev);

    //cout << unpackEvent->GetFEB12LeadTime() << endl;
    //skip FEB12 triggers that are within 1.8us of the previous trigger
    if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 700)
      continue;

    Int_t Qmax0 = 0;
    Int_t Qmax1 = 0;
    Int_t Qmax2 = 0;
    Int_t QmaxEl0[2] = {0, 0};
    Int_t QmaxEl1[2] = {0, 0};
    Int_t QmaxEl2[2] = {0, 0};
    Int_t highChargeHits = 0;

    //Skip event if there are less than 5 hits in the event
    if (mppc_hits.size()<3) 
    {
      eventsMissed += 1;
      continue;
    }
    //cout << "mppc hits = " << mppc_hits.size() << endl; 

    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      //time cut 
      if (hit->GetDt() < -326 || hit->GetDt() > 360) continue;
      //cut on charge reduces crosstalk hits
      if (hit->GetPE() > 15){
        highChargeHits++;
        if (hit->GetView() == 0){
          hXY->Fill(hit->GetX(),hit->GetY(),hit->GetTfromSpill());
          if (hit->GetPE() > Qmax0){
             Qmax0 = hit->GetPE();
             QmaxEl0[0] = hit->GetX();
             QmaxEl0[1] = hit->GetY();
          }
        }
        if (hit->GetView() == 1){
          hXZ->Fill(hit->GetX(), hit->GetZ(), hit->GetTfromSpill());
          if (hit->GetPE() > Qmax1){
            Qmax1 = hit->GetPE();
            QmaxEl1[0] = hit->GetX();
            QmaxEl1[1] = hit->GetZ();
          }
        }
        else{
          hZY->Fill(hit->GetZ(), hit->GetY(), hit->GetTfromSpill());
          if (hit->GetPE() > Qmax2)
          {
            Qmax2 = hit->GetPE();
            QmaxEl2[0] = hit->GetZ();
            QmaxEl2[1] = hit->GetY();
          }
        }
        //cout << hit->GetSpillTime()*2.5 << endl;
      }
      prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
    }

    if (Qmax0 > 0) {
      QmaxX += QmaxEl0[0];
      QmaxY += QmaxEl0[1];
      QmaxNum += 1;
    }

    // Plotting //Events__________________________________________________//______________________________

    if (highChargeHits < 3) {
      eventsMissed += 1;
      continue;
    }
    cout << "Micropulses missed: " << eventsMissed << endl;
    eventsMissed = 0;

    TCanvas *cc = new TCanvas("cc", "cc", 800, 800);

    cc->Divide(2,2);

    cc->cd(1);
    hXY->Draw("COLZ");
    hXY->SetMinimum(hXY->GetMaximum()-5);
    hXY->GetZaxis()->SetTitleOffset(2.3);
    gPad->SetRightMargin(0.2);
    cc->cd(2);
    hXZ->Draw("COLZ");
    hXZ->SetMinimum(hXZ->GetMaximum()-5);
    hXZ->GetZaxis()->SetTitleOffset(2.3);
    gPad->SetRightMargin(0.2);
    cc->cd(3);
    hZY->Draw("COLZ");
    hZY->SetMinimum(hZY->GetMaximum()-5);
    hZY->GetZaxis()->SetTitleOffset(2.3);
    gPad->SetRightMargin(0.2);
    cc->Update();
    cc->WaitPrimitive();

    hXY->Reset();
    hXZ->Reset();
    hZY->Reset();

    //cc->SaveAs((foutDir+"/event_"+to_string(iev)+".png").c_str());

    for (size_t i=0; i<mppc_hits.size();i++) delete mppc_hits[i];
  }
  cout<<"_Plotting events...100% done"<<endl;

  cout << "Number of neutrons: " << nNeutrons << endl;
  cout << "Average XY position: (" << QmaxX/QmaxNum << ", " << QmaxY/QmaxNum << ")" << endl;

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
