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
#include "TFitResult.h" 

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool batch = false;
bool IsMC = true;
bool IsReversed = false;
bool IsTwenty = false;
size_t maxEvents = std::numeric_limits<size_t>::max();
;
size_t maxSelEvents = std::numeric_limits<size_t>::max();
;
int selEvents = 0;
float start_time = clock();
bool RM_CROSSTALK = false;
bool SHOW_TRUE = true;
bool SHOW_RECO = true;
TString fileOut = "~/Desktop/templateAna_Out.root";
TString fileIn = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
int IsPrototype = true;
int SFGD_X = 204;
int SFGD_Y = 56;
int SFGD_Z = 188;
bool useClean = false;
bool useNN = false;
int dTimeIni = -100;
int dTimeFin = -100;
size_t evtIni = 0;
size_t evtFin = 0;

TFile *FileInput;
TFile *FileOutput;
TTree *dataOut;
TTree *data;
ND280SFGDEvent *inputEvent;
ND280SFGDEvent *recoEvent;
TBranch *recoBranch;
TBranch *inputBranch;
size_t nEvents;
Event *unpackEvent;

#include "../src/tools/global_tools.cc"
#include "../src/tools/reconstruction_methods.hh"

void nEventDisplay()
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
  //string foutPNG;

  cout << "total events = " << nEvents << endl;
  cout << "displaying " << evtFin - evtIni << " events" << endl;

  // Filling Histograms_______________________________________________________________________________

  Int_t prevEventTime = -999;

  TH1F *extinctionProfile = new TH1F("Neutron Extinction Profile", "Z layer of the first hit of each neutron event;Z [cm];Number of Entries", 47,1,48);
  
  //gStyle->SetOptStat(0);
  //loop over events
  for (size_t iev = evtIni; iev < evtFin; iev++)
  {
    if (iev == maxEvents - 1)
      break;
    cout << "_Processing events..." << (int)(iev * 100 / (Double_t)(evtFin - evtIni)) << "% done         \r" << flush;

    //cout<<"_Plotting events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;

    Int_t highChargeHits = 0;
    Double_t zLayer = 49;

    mppc_hits = getEventMPPChits(iev);

    //cout << unpackEvent->GetFEB12LeadTime() << endl;
    //skip FEB12 triggers that are within 1.8us of the previous trigger
    if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 700)
      continue;

    //Skip event if there are less than 3 hits in the event
    if (mppc_hits.size() < 3)
    {
      continue;
    }
    //cout << "mppc hits = " << mppc_hits.size() << endl;

    //loop over hits
    for (size_t ihit = 0; ihit < mppc_hits.size(); ihit++)
    {
      ND280SFGDHit *hit = mppc_hits[ihit];
      //time window for neutrons
      if (hit->GetDt() < -335 || hit->GetDt() > 356)
        continue;
      //fiducial cuts and fibre travel time calculation
      if (hit->GetView() == 0) 
      {
        if (hit->GetX() < 6 || hit->GetX() > 10 || hit->GetY() < 3 || hit->GetY() > 6) continue;
      }
      if (hit->GetView() == 1)
      {
        if (hit->GetX() < 6 || hit->GetX() > 10) continue;
      }
      if (hit->GetView() == 2)
      {
        if (hit->GetY() < 3 || hit->GetY() > 6) continue;
      }
      //cut on charge reduces crosstalk hits and type 3 mppc dark counts
      if (hit->GetPE() > 40)
      {
        highChargeHits += 1; 
        if (hit->GetView() == 1)
        {
          if (zLayer > hit->GetZ())
          {
            //earliestHitTime = hit->GetTfromSpill();
            zLayer = hit->GetZ();
            //extinctionProfile->Fill(hit->GetZ(),1);
          }
        }
      }
    }
    if (highChargeHits>2){
      extinctionProfile->Fill(zLayer,1);
    }
    prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
    for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
  }
  cout << "_Processing events...100% done" << endl;

  // Plotting //Events__________________________________________________//______________________________

  TCanvas *cc = new TCanvas("cc","cc",0,0,800,800);
  //cc->cd();
  cc->SetRightMargin(0.09);
  cc->SetLeftMargin(0.15);
  //extinctionProfile->Draw();
  TF1 *f1 = new TF1("f1","[0]*exp(-[1]*x)");
  //f1->SetParLimits(1, 0.03, 0.04);
  TFitResultPtr r = extinctionProfile->Fit("f1","SLM");
  Double_t extCoeff = r->Parameter(1);
  cout << "Extinction Coefficient = " << extCoeff << endl;

  cc->Update();
  cc->WaitPrimitive();

  //cc->Clear();
  //timeSpec->Draw();

  //cc->Update();
  //cc->WaitPrimitive();

  //cc->SaveAs((foutDir+"/event_"+to_string(iev)+".png").c_str());

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
