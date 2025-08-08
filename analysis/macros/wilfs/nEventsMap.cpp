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
 /*
  TH2F *EventsMapXY = new TH2F("allEvents_hitMapXY","Events per channel in XY plane;X position [cm];Y position [cm]; Number of hits",24,-12,12, 24,-12,12);
  TH2F *EventsMapXZ = new TH2F("allEvents_hitMapXZ","Events per channel in XZ plane;X position [cm];Z position [cm]; Number of hits",24,-12,12, 100,-50,50);
  TH2F *EventsMapZY = new TH2F("allEvents_hitMapZY","Events per channel in ZY plane;Z position [cm];Y position [cm]; Number of hits",100,-50,50, 24,-12,12);
  TH2F *EventsMapXYQ = new TH2F("allEvents_hitMapXYQ","Average charge per channel in XY plane;X position [cm];Y position [cm]; Average charge [p.e.]",24,-12,12, 24,-12,12);
  TH2F *EventsMapXZQ = new TH2F("allEvents_hitMapXZQ","Average charge per channel in XZ plane;X position [cm];Z position [cm]; Average charge [p.e.]",24,-12,12, 100,-50,50);
  TH2F *EventsMapZYQ = new TH2F("allEvents_hitMapZYQ","Average charge per channel in ZY plane;Z position [cm];Y position [cm]; Average charge [p.e.]",100,-50,50, 24,-12,12);
*/  
  TH1F* h_1D1 = new TH1F("Z dist","Z dist; Z position [cm]; Counts",48,0,48);  
  TH1F* h_dt = new TH1F("time", "dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_energy = new TH1F("energy","energy; E [MeV]; Counts",100,0,800);

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

    //cout<<"fRange of event "<<iev<<" is : "<<unpackEvent->GetRange()<<endl;
    //skip FEB12 triggers that are within 1.8us of the previous trigger
    //if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 700)
    //  continue;
/*
    if (mppc_hits.size() < 20 || unpackEvent->GetMaxCharge() < 50 || unpackEvent->GetEmptyZ() > 2 || unpackEvent->GetXStdDev() > 0.5 || unpackEvent->GetYStdDev() > 0.5){
      for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
      continue;
    }
*/
    // Int_t Qmax0 = 0;
    // Int_t Qmax1 = 0;
    // Int_t Qmax2 = 0;
    // Int_t QmaxEl0[2] = {0, 0};
    // Int_t QmaxEl1[2] = {0, 0};
    // Int_t QmaxEl2[2] = {0, 0};
    //Int_t highChargeHits = 0;

    //Skip event if there are less than 3 hits in the event
    //if (mppc_hits.size()<10) continue; 

    //loop over hits
    //
    /*
    int firstZ = 0;
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      if (hit->GetDt() < -326 || hit->GetDt() > 358)
        continue;
      if (hit->GetPE() < 20) continue;
      if (hit->GetView() == 2){
	if(hit-> GetY()<6 && hit->GetY()>3 && firstZ < hit->GetZ()) firstZ = hit->GetZ();
      }
    }
    h_1D1->Fill(48-firstZ);
    */
    //
    Double_t time = 1e10;
    Double_t energy = -1;
    Double_t l = 90;
    Double_t m = 939.5654133;
    Double_t c = 299792458; //speed of light in m/s
    Double_t gammaToF = l/c;
    Double_t offset = 352;

    double timm = 1e9;
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      //time cut 
      if (hit->GetDt() < -326 || hit->GetDt() > 360) continue;
      if (hit->GetPE() > 20)
      {
        cout<<"hit time "<<hit->GetDt()<<endl;
        Double_t ToF = (hit->GetDt()+offset)*2.5*1e-9+gammaToF;
        Double_t v = l/(ToF); //neutron velocity using time of flight and travel distance
        // use the ToF and energy of the earliest high charge hit
        if (time > ToF)
        {
          energy = m/sqrt(1-pow(v/c,2)) - m;
          time = ToF;
	  timm = hit->GetDt();
        }
        //cout<<"ToF, v, energy : "<<ToF<<" "<<v<<" "<<energy<<endl;
      }
    }
    cout<<"time and energy : "<<time<<" "<<energy<<endl;
    h_dt->Fill(timm);
    h_energy ->Fill(energy);

    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      if (hit->GetDt() < -326 || hit->GetDt() > 358)
        continue;
      if (hit->GetPE() < 20) continue;

      if (hit->GetView() == 0){
        //if (hit->GetX()<24 && hit->GetX()>-1 && hit->GetY()<8 && hit->GetY()>-1){
          // if (hit->GetPE() > Qmax0){
            // highChargeHits++;
            // Qmax0 = hit->GetPE();
            // QmaxEl0[0] = hit->GetX();
            // QmaxEl0[1] = hit->GetY();
          // }
          EventsMapXY->Fill(hit->GetX(),hit->GetY(),1);
          EventsMapXYQ->Fill(hit->GetX(),hit->GetY(),hit->GetPE());
        //}
      }
      else if (hit->GetView() == 1){
        //if (hit->GetX()<24 && hit->GetX()>-1 && hit->GetZ()<48 && hit->GetZ()>-1){
          // if (hit->GetPE() > Qmax1){
            // Qmax1 = hit->GetPE();
            // QmaxEl1[0] = hit->GetX();
            // QmaxEl1[1] = hit->GetZ();
          // }
          EventsMapXZ->Fill(hit->GetX(),hit->GetZ(),1);
          EventsMapXZQ->Fill(hit->GetX(),hit->GetZ(),hit->GetPE());
        //}
      }
      else if (hit->GetView() == 2){
        //if (hit->GetZ()<48 && hit->GetZ()>-1 && hit->GetY()<8 && hit->GetY()>-1){
          // if (hit->GetPE() > Qmax2){
            // Qmax2 = hit->GetPE();
            // QmaxEl2[0] = hit->GetZ();
            // QmaxEl2[1] = hit->GetY();
          // }
          EventsMapZY->Fill(hit->GetZ(),hit->GetY(),1);
          EventsMapZYQ->Fill(hit->GetZ(),hit->GetY(),hit->GetPE());
          //if(hit-> GetY()<6 && hit->GetY()>3){
            //h_1D1->Fill(48-hit->GetZ(), hit->GetPE());
	    //break;
	  //}
	//}
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

  gStyle->SetOptStat(0);

  new TCanvas();
  h_1D1->Draw();

  new TCanvas();
  h_dt->Draw();

  new TCanvas();
  h_energy->Draw();

  TCanvas *cc = new TCanvas("cc","cc", 0, 0, 1800, 900);

  cc->Divide(3,2);

  cc->cd(1);
  EventsMapXY->Draw("COLZ");
  EventsMapXY->GetZaxis()->SetTitleOffset(1.5);
  gPad->SetRightMargin(0.15);
  cc->cd(2);
  EventsMapXZ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapXZ->GetZaxis()->SetTitleOffset(1.5);
  cc->cd(3);
  EventsMapZY->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapZY->GetZaxis()->SetTitleOffset(1.5);
  cc->cd(4);
  EventsMapXYQ->Divide(EventsMapXY);
  EventsMapXYQ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapXYQ->GetZaxis()->SetTitleOffset(1.5);
  cc->cd(5);
  EventsMapXZQ->Divide(EventsMapXZ);
  EventsMapXZQ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapXZQ->GetZaxis()->SetTitleOffset(1.5);
  cc->cd(6);
  EventsMapZYQ->Divide(EventsMapZY);
  EventsMapZYQ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapZYQ->GetZaxis()->SetTitleOffset(1.5);

  cc->Update();
  cc->WaitPrimitive();

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
