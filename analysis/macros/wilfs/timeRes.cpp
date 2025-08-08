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

  TH1F *histTime = new TH1F("histTime","beam center MIP-like X; time difference [2.5 ns]; number of events", 40,-20,20);
  TH1F *histTime2 = new TH1F("histTime","beam center MIP-like sec. half X; time difference [2.5 ns]; number of events", 40,-20,20);
  TH1F *histTime3 = new TH1F("histTime","beam center higher-than MIP X; time difference [2.5 ns]; number of events", 40,-20,20);
  TH1F *histTime4 = new TH1F("histTime","beam center MIP-like Y; time difference [2.5 ns]; number of events", 40,-20,20);

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

    int counter1=0;
    double ctime[1000]={};
    double clayer[1000]={};
    double ccharge[1000]={};
    double tott=0;
    bool filled=false;

    int counter2=0;
    double ctime2[1000]={};
    double clayer2[1000]={};
    double ccharge2[1000]={};
    double tott2=0;
    bool filled2=false;

    int counter3=0;
    double ctime3[1000]={};
    double clayer3[1000]={};
    double ccharge3[1000]={};
    double tott3=0;
    bool filled3=false;

    int counter4=0;
    double ctime4[1000]={};
    double clayer4[1000]={};
    double ccharge4[1000]={};
    double tott4=0;
    bool filled4=false;    

    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      if (hit->GetDt() < -326 || hit->GetDt() > 358)
        continue;

      if (hit->GetPE() < 20) continue;
      if (hit->GetView() == 0 && ! (hit->GetX()==15 && hit->GetY()==4 ) ) break;

      tott += hit->GetPE();
      if(hit->GetView() == 2){
        ctime[counter1] = hit->GetDt(); 
	clayer[counter1] = hit->GetZ(); 
	ccharge[counter1] = hit->GetPE();
	counter1 ++;
      }

      if(hit->GetView() == 2 && hit->GetZ()<24){
        ctime2[counter2] = hit->GetDt();
        clayer2[counter2] = hit->GetZ();
        ccharge2[counter2] = hit->GetPE();
        counter2 ++;
      }

      if(hit->GetView() == 2){
        ctime3[counter3] = hit->GetDt();
        clayer3[counter3] = hit->GetZ();
        ccharge3[counter3] = hit->GetPE();
        counter3 ++;
      }

      if(hit->GetView() == 1){
        ctime4[counter4] = hit->GetDt();
        clayer4[counter4] = hit->GetZ();
        ccharge4[counter4] = hit->GetPE();
        counter4 ++;
      }

      prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
      //cout << hit->GetSpillTime()*2.5 << endl;
    }    //Draw histograms___________________________________________________________________________________

    for(int loop1=0;loop1<100;loop1++){
      for(int loop2=loop1;loop2<100;loop2++){
        if(filled == false && ccharge[loop1]>40 && ccharge[loop1]<60 && abs(clayer[loop1]-clayer[loop2]) == 1 && abs(ccharge[loop1] - ccharge[loop2]) < 5) {
	  histTime->Fill((double) ctime[loop1]-ctime[loop2]);
	  filled = true;
	  //std::cout<<ctime[loop1]<<" "<<ctime[loop2]<<std::endl;
	}
        if(filled2 == false && ccharge2[loop1]>40 && ccharge2[loop1]<60 && abs(clayer2[loop1]-clayer2[loop2]) == 1 && abs(ccharge2[loop1] - ccharge2[loop2]) < 5) {
          histTime2->Fill((double) ctime2[loop1]-ctime2[loop2]);
          filled2 = true;
        }
        if(filled3 == false && ccharge3[loop1]>60 && ccharge3[loop1]<90 && abs(clayer3[loop1]-clayer3[loop2]) == 1 && abs(ccharge3[loop1] - ccharge3[loop2]) < 5) {
          histTime3->Fill((double) ctime3[loop1]-ctime3[loop2]);
          filled3 = true;
        }
        if(filled4 == false && ccharge4[loop1]>40 && ccharge4[loop1]<60 && abs(clayer4[loop1]-clayer4[loop2]) == 1 && abs(ccharge4[loop1] - ccharge4[loop2]) < 5) {
          histTime4->Fill((double) ctime4[loop1]-ctime4[loop2]);
          filled4 = true;
	  //std::cout<<ctime4[loop1]<<" "<<ctime4[loop2]<<std::endl;
        }
	
      }
    }

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
  histTime->Draw();
  //cc->Update();
  //cc->WaitPrimitive();

  TCanvas *cc2 = new TCanvas();
  histTime2->Draw();
  //cc2->Update();
  //cc2->WaitPrimitive();

  TCanvas *cc3 = new TCanvas();
  histTime3->Draw();
  //cc3->Update();
  //cc3->WaitPrimitive();

  TCanvas *cc4 = new TCanvas();
  histTime4->Draw();
  //cc4->Update();
  //cc4->WaitPrimitive();  

  //FileOutput->Close();
  //FileInput->Close();
  //handleEndOfExecution();
}
