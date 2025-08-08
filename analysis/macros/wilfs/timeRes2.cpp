#define THIS_NAME nEventsMap
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS
#define _CRTDBG_MAP_ALLOC

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
int     maxEvents = std::numeric_limits<int>::max();;
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
size_t evtIni             = 0;
size_t evtFin             = 0;

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

void nEventsMap() {


  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent,
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  cout<<"total events = "<<nEvents<<endl;
  cout<<"displaying "<<evtFin-evtIni<<" events"<<endl;

  // Create Histograms________________________________________________________________________________

  TH2F *EventsMapXY = new TH2F("allEvents_hitMapXY","Events per channel in XY plane;X position [cm];Y position [cm]; Number of hits",24,0,24, 8,0,8);
  TH2F *EventsMapXZ = new TH2F("allEvents_hitMapXZ","Events per channel in XZ plane;X position [cm];Z position [cm]; Number of hits",24,0,24, 48,0,48);
  TH2F *EventsMapZY = new TH2F("allEvents_hitMapZY","Events per channel in ZY plane;Z position [cm];Y position [cm]; Number of hits",48,0,48, 8,0,8);
  TH2F *EventsMapXYQ = new TH2F("allEvents_hitMapXYQ","Average charge per channel in XY plane;X position [cm];Y position [cm]; Average charge [p.e.]",24,0,24, 8,0,8);
  TH2F *EventsMapXZQ = new TH2F("allEvents_hitMapXZQ","Average charge per channel in XZ plane;X position [cm];Z position [cm]; Average charge [p.e.]",24,0,24, 48,0,48);
  TH2F *EventsMapZYQ = new TH2F("allEvents_hitMapZYQ","Average charge per channel in ZY plane;Z position [cm];Y position [cm]; Average charge [p.e.]",48,0,48, 8,0,8);

  TH1F *h_dt = new TH1F("h_dt","dt dist.; dt [2.5ns]; number of events", 800,-400,400);

  // Filling Histograms_______________________________________________________________________________

  //Int_t nNeutrons = 0;
  //Int_t QmaxX = 0;
  //Int_t QmaxY = 0;
  //Int_t QmaxZ = 0;
  //Int_t QmaxNum = 0;
  Int_t prevEventTime = 0;
  //loop over events
  for (size_t iev=evtIni; iev<evtFin; iev++){
    if (iev == maxEvents-1) break;
    cout<<"_Processing events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;
    
    std::vector<ND280SFGDHit*> mppc_hits;
    mppc_hits = getEventMPPChits(iev);

    //skip FEB12 triggers that are within 1.8us of the previous trigger
    if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 700)
      continue;

    if (mppc_hits.size() < 20 || unpackEvent->GetMaxCharge() < 50 || unpackEvent->GetEmptyZ() > 2 || unpackEvent->GetXStdDev() > 0.5 || unpackEvent->GetYStdDev() > 0.5){
      for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
      continue;
    }

    // Int_t Qmax0 = 0;
    // Int_t Qmax1 = 0;
    // Int_t Qmax2 = 0;
    // Int_t QmaxEl0[2] = {0, 0};
    // Int_t QmaxEl1[2] = {0, 0};
    // Int_t QmaxEl2[2] = {0, 0};
    //Int_t highChargeHits = 0;

    //Skip event if there are less than 3 hits in the event
    //if (mppc_hits.size()<10) continue; 
    bool filled1 = false;

    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      if (hit->GetDt() < -360 || hit->GetDt() > 358)
        continue;

      if (hit->GetPE() < 40) continue;

      if (hit->GetView() == 1 && hit->GetZ() == 47 && hit->GetX()==15 && filled1==false){
	std::cout<<hit->GetZ()<<std::endl;
	if(hit->GetDt()<300)
          h_dt->Fill(hit->GetDt());
	else if (hit->GetDt()>300)
	  h_dt->Fill(hit->GetDt()-(698+17));
	filled1 = true;
      }	

      prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
      //cout << hit->GetSpillTime()*2.5 << endl;
    }    //Draw histograms___________________________________________________________________________________

    // if (Qmax0 > 0 && Qmax1 > 0 && Qmax2 > 0) {
      //nNeutrons++;
      //QmaxX += QmaxEl0[0];
      //QmaxY += QmaxEl0[1];
      //QmaxNum += 1;
      //cout << Qmax0 << " " << Qmax1 << " " << Qmax2 << endl;
      //cout << "(" << QmaxEl0[0] << ", " << QmaxEl0[1] << ")" << "(" << QmaxEl1[0] << ", " << QmaxEl1[1] << ")" << "(" << QmaxEl2[0] << ", " << QmaxEl2[1] << ")" << endl;
      // EventsMapXYmaxQ->Fill(QmaxEl0[0],QmaxEl0[1],1);
      // EventsMapXZmaxQ->Fill(QmaxEl1[0],QmaxEl1[1],1);
      // EventsMapZYmaxQ->Fill(QmaxEl2[0],QmaxEl2[1],1);
    // }
    for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
  }
  cout<<"_Processing events...100% done"<<endl;

  //cout << "Number of neutrons: " << nNeutrons << endl;
  //cout << "Average XY position: (" << QmaxX/QmaxNum << ", " << QmaxY/QmaxNum << ")" << endl;

  // Plotting //Events__________________________________________________//______________________________

  TCanvas *cc = new TCanvas();
  h_dt->Draw();

  TF1* g1 = new TF1("g1","gaus",-360,-355);
  TF1* g2 = new TF1("g2","expo",-356,-354);
  TF1* total = new TF1("total","gaus(0)+expo(3)",-358,-354);
  h_dt->Fit(g1,"R");
  h_dt->Fit(g2,"R+");
  g1->SetLineColor(2);
  g2->SetLineColor(1);
  //g1->Draw("same");
  //g2->Draw("same");

  double par[5];
  par[0] = g1->GetParameter(0);
  par[1] = g1->GetParameter(1);
  par[2] = g1->GetParameter(2);
  par[3] = g2->GetParameter(0);
  par[4] = g2->GetParameter(1);

  for(int i=0;i<5;i++) cout<<par[i]<<endl;

  total->SetParameters(par);
  //total->SetParameter(4,3);
  //total->FixParameter(4,3);
  //for(int i=0;i<5;i++)
  //total->FixParameter(i,par[i]);

  h_dt->Fit(total,"R+");
  total->SetLineColor(4);
  total->SetLineWidth(2);
  //total->Draw("same");

/*
  TF1* f1 = new TF1("f1","[0]*exp( -0.5 * ((x-[1])/[2])*((x-[1])/[2])) + [3] * exp(-[4]*(x+360))",-360,-355);
  f1->SetParameter(0,1000);
  f1->SetParameter(1,-357);
  f1->SetParameter(2,1);
  f1->SetParameter(3,1);
  f1->SetParameter(4,1);
  f1->FixParameter(3,1);
  f1->FixParameter(4,1);
  h_dt->Fit("f1");
  f1->SetLineColor(2);
  f1->SetLineWidth(2);
  f1->Draw("same");
*/
  //cc->Update();
  //cc->WaitPrimitive();

  //FileOutput->Close();
  //FileInput->Close();
  //handleEndOfExecution();
}
