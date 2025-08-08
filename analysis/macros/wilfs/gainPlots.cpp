#define THIS_NAME gainPlots
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
#include <fstream>

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

void gainPlots()
{

  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent,
  // recoBranch, inputBranch, nEvents, unpackEvent

  std::vector<ND280SFGDHit *> mppc_hits;

  //ostringstream foutPNGnum;
  //string foutPNG// set t2k style for plots
  
  TString localStyleName = "T2K";
  Int_t localWhichStyle = 3;
  TStyle *t2kstyle = SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->ForceStyle(t2kstyle);;

  // Initialise Histograms_______________________________________________________________________________

  Int_t FEBs[18] = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};

  TH1D* gain1 = new TH1D("gain1",";Gain [ADC/p.e.];Number of Entries", 50, 0, 50);
  TH1D* gain2 = new TH1D("gain2",";Gain [ADC/p.e.];Number of Entries", 50, 0, 50);
  TH1D* gain3 = new TH1D("gain3",";Gain [ADC/p.e.];Number of Entries", 50, 0, 50);

  for (int mcr=0; mcr<4; mcr++){
    for (int slot=0; slot<5; slot++){
      if (mcr != 0 && mcr != 2 && slot == 4) continue;
      ifstream inFile(("../../../data_preprocessing/bin/calib_files/HG50/HG_PE/MCR"+to_string(mcr)+"_Slot_"+to_string(slot)+"_gain.txt").c_str());
      if (!inFile) cout << "Could not find file" << endl;
      int i = 0;
      double x;
      while (inFile >> x){
        if (i%4 == 1){
          if (x != 0){
            if (mcr == 1){
              if (slot == 0 || slot == 3)
                gain3->Fill(x);
              else
                gain2->Fill(x);
            }
            else if(mcr == 0){
              if (slot == 3 || slot == 4)
                gain2->Fill(x);
              else
                gain1->Fill(x);
            }
            else gain1->Fill(x);
          }
        }
        i++;
      }
      inFile.close();
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *cc = new TCanvas("cc","cc",800,630);
  gain1->SetLineColor(kGreen);
  gain1->SetLineStyle(1);
  gain1->SetLineWidth(3);
  gain2->SetLineColor(kBlue);
  gain2->SetLineStyle(7);
  gain2->SetLineWidth(3);
  gain3->SetLineColor(kRed);
  gain3->SetLineStyle(2);
  gain3->SetLineWidth(3);
  // gain1->Scale(1/gain1->Integral("width"));
  gain1->Draw("hist");
  cc->Update();
  // cc->SaveAs("gain1.C");
  // cc->Clear();
  // gain2->Scale(1/gain2->Integral("width"));
  gain2->Draw("same hist");
  cc->Update();
  // cc->SaveAs("gain2.C");
  // cc->Clear();
  // gain3->Scale(1/gain3->Integral("width"));
  gain3->Draw("same hist");
  TLegend *leg = new TLegend(5,30,20,200);
  leg->AddEntry(gain1,"Type 1","l");
  leg->AddEntry(gain2,"Type 2","l");
  leg->AddEntry(gain3,"Type 3","l");
  leg->Draw();
  cc->Update();
  TPaveText *pave1 = new TPaveText(5,30,20,200);
  Double_t entries1 = gain1->GetEntries();
  Double_t mean1 = gain1->GetMean();
  Double_t std1 = gain1->GetStdDev();
  pave1->AddText("");
  pave1->AddText(Form("Entries %g",entries1));
  pave1->AddText(Form("Mean %.2f ns",mean1));
  pave1->AddText(Form("Std Dev %.2f ns",std1));
  pave1->AddLine(0,0.33,1,0.33);
  Double_t entries2 = gain2->GetEntries();
  Double_t mean2 = gain2->GetMean();
  Double_t std2 = gain2->GetStdDev();
  pave1->AddText("");
  pave1->AddText(Form("Entries %g",entries2));
  pave1->AddText(Form("Mean %.2f ns",mean2));
  pave1->AddText(Form("Std Dev %.2f ns",std2));
  pave1->AddLine(0,0.66,1,0.66);
  Double_t entries3 = gain3->GetEntries();
  Double_t mean3 = gain3->GetMean();
  Double_t std3 = gain3->GetStdDev();
  pave1->AddText("");
  pave1->AddText(Form("Entries %g",entries3));
  pave1->AddText(Form("Mean %.2f ns",mean3));
  pave1->AddText(Form("Std Dev %.2f ns",std3));
  pave1->Paint("NDC");
  pave1->Draw("same");
  cc->Modified(); cc->Update();
  cc->SaveAs("gain.C");


  //cc->SaveAs((foutDir+"/event_"+to_string(iev)+".png").c_str());

  // FileOutput->Close();
  // FileInput->Close();
  handleEndOfExecution();
}
