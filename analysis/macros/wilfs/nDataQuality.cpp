#define THIS_NAME nDataQuality
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
bool    batch          = true;
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

void nDataQuality() {

  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent, 
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  std::vector<ND280SFGDHit*> mppc_hits;

  cout << "Processing " << evtFin-evtIni << " of " << nEvents << " events" << endl;

  // Create Histograms________________________________________________________________________________

  Int_t FEB[18] = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};
  TH1F *spillsPerFEB = new TH1F("SpillsPerFEB","Spills per FEB;FEB Number;Number of Spills", 28,0,28);
  TH1F *hitsPerFEB = new TH1F("hitsPerFEB","Hits per FEB;FEB Number;Number of Hits",28,0,28);
  TH2F *trigTimeVHitTime = new TH2F("trigTimeVHitTime","FEB12 hit time since spill versus all other FEB hit times since spill; FEB12 hit time since spill [2.5 ns]; Hit time since spill [2.5 ns]; Number of entries", 3000, 0, 300000, 3000, 0, 300000);
  TH1F *FEBsync[28];
  TH1F *hitsPerFEBch[28];
  TH1F *eventsPerSpillAtFEB[28];
  TH1F *hitsPerChPerSpill[28];
  TH1F *hitsPerFEBchDist[28];
  TH1F *hitTimeSinceSpill[28];
  for (size_t i=0; i<18; i++){
    FEBsync[FEB[i]] = new TH1F(("FEB"+to_string(FEB[i])+" synchronisation").c_str(),("FEB"+to_string(FEB[i])+" synchronisation;GTrigTime since spill [2.5 ns];Number of Entries").c_str(),6000,0,6000);
    hitsPerFEBch[FEB[i]] = new TH1F(("FEB"+to_string(FEB[i])+" hits per channel").c_str(),("FEB"+to_string(FEB[i])+" hits per channel;Channel Number;Number of Hits").c_str(),96,0,96);
    eventsPerSpillAtFEB[FEB[i]] = new TH1F(("FEB"+to_string(FEB[i])+" events per spill").c_str(),("FEB"+to_string(FEB[i])+" events per spill;Number of events per spill;Number of entries").c_str(),1000,0,20000);
    hitsPerChPerSpill[FEB[i]] = new TH1F(("FEB"+to_string(FEB[i])+" distribution of hits per channel per spill").c_str(),("FEB"+to_string(FEB[i])+" hits per channel per spill;Number of hits at channel;Number of entries").c_str(),10000,0,10000);
    hitsPerFEBchDist[FEB[i]] = new TH1F(("FEB"+to_string(FEB[i])+" distribution of hits per channel").c_str(),("FEB"+to_string(FEB[i])+" distribution of hits per channel;Number of hits at channel;Number of entries").c_str(),100,0,100000);
    hitTimeSinceSpill[FEB[i]] = new TH1F(("FEB"+to_string(FEB[i])+" hit time since spill start").c_str(),("FEB"+to_string(FEB[i])+" hit time since spill start;Hit time since spill start [2.5 ns];Number of entries").c_str(),300000,0,300000);
  }

  // Filling Histograms_______________________________________________________________________________

  Int_t prevSpill[28] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  Int_t prevGTrigTag[28] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  Int_t eventsPerSpill[28] = {0};
  Int_t hitsPerChPerSpillarray[28][96] = {0};

  //loop over events
  for (size_t iev=evtIni; iev<evtFin; iev++){
    if (iev == maxEvents-1) break;
    cout<<"_Processing events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;
   
    mppc_hits = getEventMPPChits(iev);

    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];

      if (hit->GetSpillTag() != prevSpill[(int)hit->GetFEB()])
      { 
        spillsPerFEB->Fill(hit->GetFEB(),1);
        eventsPerSpill[(int)hit->GetFEB()] = 0;
        for (size_t ich=0; ich<96; ich++) 
        {
          hitsPerChPerSpill[(int)hit->GetFEB()]->Fill(hitsPerChPerSpillarray[(int)hit->GetFEB()][ich], 1);
          hitsPerChPerSpillarray[(int)hit->GetFEB()][ich] = 0;
          eventsPerSpillAtFEB[(int)hit->GetFEB()]->Fill(eventsPerSpill[(int)hit->GetFEB()], 1);
        }
      }
      prevSpill[(int)hit->GetFEB()] = hit->GetSpillTag();
      
      if (hit->GetGTrigTag() != prevGTrigTag[(int)hit->GetFEB()]) FEBsync[(int)hit->GetFEB()]->Fill(abs(hit->GetGTrigTime()-hit->GetSpillTime()),1);
      prevGTrigTag[(int)hit->GetFEB()] = hit->GetGTrigTag();

      hitTimeSinceSpill[(int)hit->GetFEB()]->Fill(hit->GetTfromSpill());

      hitsPerFEB->Fill((int)hit->GetFEB(),1);
      hitsPerFEBch[(int)hit->GetFEB()]->Fill(hit->GetCh(),1);
      hitsPerChPerSpillarray[(int)hit->GetFEB()][(int)hit->GetCh()] += 1;
      trigTimeVHitTime->Fill(hit->GetTfromSpill() - hit->GetDt(), hit->GetTfromSpill(), 1);
    

      eventsPerSpill[(int)hit->GetFEB()] += 1;

    }
    for (size_t i=0;i<mppc_hits.size();i++) delete mppc_hits[i];
  }
  cout<<"_Processing events...100% done"<<endl;

  // Plotting Events________________________________________________________________________________

  gStyle->SetOptStat(0);
  TDirectory *FEBsyncDir = FileOutput->mkdir("FEB Synchronisation");
  TDirectory *hitsPerFEBchDir = FileOutput->mkdir("Hits per FEB channel");
  TDirectory *hitsPerFEBchDistDir = FileOutput->mkdir("Distribution of channel hits per FEB");
  TDirectory *hitsPerChPerSpillDir = FileOutput->mkdir("Distribution of channel hits per FEB per spill");
  TDirectory *eventsPerSpillAtFEBDir = FileOutput->mkdir("Events per spill at FEB");
  TDirectory *hitTimeSinceSpillDir = FileOutput->mkdir("Hit time since spill per FEB");

  //TCanvas *cc = new TCanvas("cc","cc", 0, 0, 1800, 900);

  FEBsyncDir->cd();
  for (size_t i=0; i<18; i++){
    Int_t maxBin = FEBsync[FEB[i]]->FindLastBinAbove(0,1);
    FEBsync[FEB[i]]->GetXaxis()->SetRange(0, maxBin);
    FEBsync[FEB[i]]->Write(("FEB"+to_string(FEB[i])+" synchronisation").c_str());
  }

  hitsPerFEBchDir->cd();
  for (size_t i=0; i<18; i++){
    hitsPerFEBchDir->cd();
    hitsPerFEBch[FEB[i]]->Write(("FEB"+to_string(FEB[i])+" hits per channel").c_str());
    for (Int_t ich=0; ich<96; ich++)
      hitsPerFEBchDist[FEB[i]]->Fill(hitsPerFEBch[FEB[i]]->GetBinContent(ich),1);
    hitsPerFEBchDistDir->cd();
    hitsPerFEBchDist[FEB[i]]->GetXaxis()->SetRange(0,hitsPerFEBchDist[FEB[i]]->FindLastBinAbove()+10);
    hitsPerFEBchDist[FEB[i]]->Write(("FEB"+to_string(FEB[i])+" distribution of hits per channel").c_str());
  }

  hitsPerChPerSpillDir->cd();
  for (size_t i = 0; i < 18; i++)
  {
    Int_t maxBin = hitsPerChPerSpill[FEB[i]]
                       ->FindLastBinAbove(0, 1);
    hitsPerChPerSpill[FEB[i]]->GetXaxis()->SetRange(0, maxBin + 10);
    hitsPerChPerSpill[FEB[i]]->Write(("FEB_" + to_string(FEB[i]) + " distribution of hits per channel per spill").c_str());
  }

  eventsPerSpillAtFEBDir->cd();
  for (size_t i=0; i<18; i++){
    Int_t maxBin = eventsPerSpillAtFEB[FEB[i]]
    ->FindLastBinAbove(0,1);
    eventsPerSpillAtFEB[FEB[i]]->GetXaxis()->SetRange(0, maxBin+10);
    eventsPerSpillAtFEB[FEB[i]]->Write(("FEB_"+to_string(FEB[i])+" events per spill").c_str());
  }

  hitTimeSinceSpillDir->cd();
  for (size_t i=0; i<18; i++){
    Int_t maxBin = hitTimeSinceSpill[FEB[i]]->FindLastBinAbove(0,1);
    hitTimeSinceSpill[FEB[i]]->GetXaxis()->SetRange(0, maxBin+10);
    hitTimeSinceSpill[FEB[i]]->Write(("FEB"+to_string(FEB[i])+" hit time since spill").c_str());
  }

  FileOutput->cd();
  hitsPerFEB->Write("Hits per FEB");
  spillsPerFEB->Write("Spills per FEB");
  trigTimeVHitTime->Write("Trig time vs hit time");

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
