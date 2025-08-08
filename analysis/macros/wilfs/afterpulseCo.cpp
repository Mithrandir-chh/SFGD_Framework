#define THIS_NAME afterpulseCo
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
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLegend.h"
#include "TImage.h"
#include "TLeaf.h"
#include "TColor.h"
#include "TPaletteAxis.h"
#include <vector>

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool batch = false;
bool IsMC = true;
bool IsReversed = false;
bool IsTwenty = false;
int maxEvents = std::numeric_limits<int>::max();
;
int maxSelEvents = std::numeric_limits<int>::max();
;
int selEvents = 0;
float start_time = clock();
bool RM_CROSSTALK = false;
bool SHOW_TRUE = true;
bool SHOW_RECO = true;
TString fileOut = "~/Desktop/templateAna_Out.root";
TString fileIn = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
bool IsPrototype = true;
int SFGD_X = 204;
int SFGD_Y = 56;
int SFGD_Z = 188;
bool useClean = false;
bool useNN = false;
int dTimeIni = -100;
int dTimeFin = -100;
int evtIni = 0;
int evtFin = 0;

TFile *FileInput;
TFile *FileOutput;
TTree *dataOut;
TTree *data;
ND280SFGDEvent *inputEvent;
ND280SFGDEvent *recoEvent;
TBranch *recoBranch;
TBranch *inputBranch;
int nEvents;
Event *unpackEvent;

#include "../src/tools/global_tools.cc"
#include "../src/tools/reconstruction_methods.hh"

void afterpulseCo()
{

  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent,
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  std::vector<ND280SFGDHit *> mppc_hits;

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
  //string foutPNG// set t2k style for plots
  
  TString localStyleName = "T2K";
  Int_t localWhichStyle = 3;
  TStyle *t2kstyle = SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->ForceStyle(t2kstyle);;

  cout << "total events = " << nEvents << endl;
  cout << "displaying " << evtFin - evtIni << " events" << endl;

  // Initialise Histograms_______________________________________________________________________________

  Int_t prevEventTime = -999;
  Int_t FEBs[18] = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};
  // Int_t FEBmap[28] = {0,1,2,3,4,0,0,0,5,6,7,8,0,0,0,0,9,10,11,12,13,0,0,0,14,15,16,17};
  // Int_t FEBmap1[28] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
  // Int_t FEBmap2[28] = {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,3,0,0,0};
  // Int_t FEBmap3[28] = {0,0,0,0,1,0,0,0,2,3,4,5,0,0,0,0,0,0,6,7,8,0,0,0,0,9,10,11};
  // Double_t binSize = 0.01;

  TH1D* aftCh = new TH1D("AftCh",";Afterpulse Charge [p.e.];Number of Entries", 500,0,500);
  TH1D* bef = new TH1D("bef",";Hit before afterpulse;Number of Entries", 5,0,5);
  TH1D* ch8 = new TH1D("ch8","FEB8;Channel Number;Number of Entries", 96,0,96);
  TH1D* ch11 = new TH1D("ch11","FEB11;Channel Number;Number of Entries", 96,0,96);
  TH1D* zeroCh8 = new TH1D("ZeroCh8","FEB8;Channel Number;Number of Entries", 96,0,96);
  TH1D* zeroCh11 = new TH1D("ZeroCh11","FEB11;Channel Number;Number of Entries", 96,0,96);
  
  // gStyle->SetOptStat(0);
  //loop over events
  for (int iev = evtIni; iev < evtFin; iev++)
  {
    if (iev == maxEvents - 1)
      break;
    cout << "_Processing events..." << (int)(iev * 100 / (Double_t)(evtFin - evtIni)) << "% done         \r" << flush;

    mppc_hits = getEventMPPChits(iev);

    if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 250){
      for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
      continue;
    }

    //loop over hits
    for (size_t ihit = 0; ihit < mppc_hits.size(); ihit++)
    {
      ND280SFGDHit *hit = mppc_hits[ihit];
      int before = 0;
      //ignore hits with zero or negative charge
      if (hit->GetPE() <= 0.) continue; 

      if (hit->GetFEB() != 8 && hit->GetFEB() != 11) continue;

      if (hit->GetDt() < -87 || hit->GetDt() > -70) continue;

      //check for previous hits in same channel
      for (size_t ihit2 = 0; ihit2 < mppc_hits.size(); ihit2++){
        ND280SFGDHit *hit2 = mppc_hits[ihit2];
        // if (hit2->GetPE() <= 0.) continue;
        if (hit->GetFEB() == hit2->GetFEB() && hit->GetCh() == hit2->GetCh())
          if (hit2->GetDt() > -120 && hit2->GetDt() < -87) before++;
      }
      // cout << hit->GetFEB() << " " << hit->GetCh() << endl;

      if (hit->GetFEB() == 8){
        ch8->Fill(hit->GetCh());
        if (before == 0) zeroCh8->Fill(hit->GetCh());
      }
      else{
         ch11->Fill(hit->GetCh());
          if (before == 0) zeroCh11->Fill(hit->GetCh());
      }

      aftCh->Fill(hit->GetPE());
      bef->Fill(before);
    }
    prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
    //release memory used for hit information
    for (size_t i=0;i<mppc_hits.size();i++) delete mppc_hits[i];
  }
  cout << "_Processing events...100% done" << endl;

  // Plotting Events__________________________________________________//______________________________

  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(0);
  TCanvas *cc1 = new TCanvas("cc1","cc1",0,0,800,630);
  cc1->cd();
  // TF1 *f1 = new TF1("f1","[0]+[1]*x",500,2500);
  // HGvLG->Fit(f1,"R");
  aftCh->Draw();
  // f1->Draw("same");
  cc1->SaveAs("aftCh.C");
  cc1->Clear();
  bef->Draw();
  cc1->SaveAs("bef.C");
  cc1->Clear();
  ch8->Draw();
  cc1->SaveAs("ch8.C");
  cc1->Clear();
  ch11->Draw();
  cc1->SaveAs("ch11.C");
  cc1->Clear();
  zeroCh8->Draw();
  cc1->SaveAs("zeroCh8.C");
  cc1->Clear();
  zeroCh11->Draw();
  cc1->SaveAs("zeroCh11.C");
  cc1->Clear();

  //cc->SaveAs((foutDir+"/event_"+to_string(iev)+".png").c_str());

  // FileOutput->Close();
  // FileInput->Close();
  writeOutput();
  handleEndOfExecution();
}
