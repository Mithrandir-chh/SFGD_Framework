#define THIS_NAME nEventsMapUSJ
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
#include "TGraph.h"

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

void nEventsMapUSJ() {


  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent,
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  cout<<"total events = "<<nEvents<<endl;
  cout<<"displaying "<<evtFin-evtIni<<" events"<<endl;

  Double_t l = 90;
  Double_t m = 939.5654133;
  Double_t c = 299792458; //speed of light in m/s
  Double_t gammaToF = l/c;
  //Double_t offset = 352;
  Double_t offset = 337;

  int nbins = 10000;
  Double_t xbins[nbins+1];
  for (size_t ibin=nbins+1; ibin+1>0; ibin--)
  {
    Double_t x = l/(ibin*2.5*1e-9+gammaToF);
    xbins[nbins+1-ibin] = m/sqrt(1-pow(x/c,2)) - m;
  }

  // Create Histograms________________________________________________________________________________

  TH2F *EventsMapXY = new TH2F("allEvents_hitMapXY","allEvents_hitMapXY;X position [cm];Y position [cm]; Number of hits",24,0,24, 8,0,8);
  TH2F *EventsMapXZ = new TH2F("allEvents_hitMapXZ","allEvents_hitMapXZ;X position [cm];Z position [cm]; Number of hits",24,0,24, 48,0,48);
  TH2F *EventsMapZY = new TH2F("allEvents_hitMapZY","allEvents_hitMapZY;Z position [cm];Y position [cm]; Number of hits",48,0,48, 8,0,8);
  TH2F *EventsMapXYUSJ = new TH2F("allEvents_hitMapXYUSJ","allEvents_hitMapXY_US-Japan;X position [cm];Y position [cm]; Number of hits",8,0,8, 8,0,8);
  TH2F *EventsMapXZUSJ = new TH2F("allEvents_hitMapXZUSJ","allEvents_hitMapXZ_US-Japan;X position [cm];Z position [cm]; Number of hits",8,0,8, 32,0,32);
  TH2F *EventsMapZYUSJ = new TH2F("allEvents_hitMapZYUSJ","allEvents_hitMapZY_US-Japan;Z position [cm];Y position [cm]; Number of hits",32,0,32, 8,0,8);

  TH1F* h_dt = new TH1F("time", "USJ dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_dtt = new TH1F("time", "SFGD dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_dt2 = new TH1F("time", "dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_dt3 = new TH1F("time", "dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_dt4 = new TH1F("time", "dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_dt5 = new TH1F("time", "dt; dt [2.5 ns]; Counts", 3000,-1500,1500);
  TH1F* h_dt6 = new TH1F("time", "dt; dt [2.5 ns]; Counts", 3000,-1500,1500);  
  TH1F* h_tspill = new TH1F("timeFromSpill", "Tspill; Tspill [2.5 ns]; Counts", 100e3,0,300e3);
  TH1F* h_energy = new TH1F("energy","energy; E [MeV]; Counts",nbins,xbins);

  TH1F* h_spreadX[10];
  TH1F* h_spreadY[10];
  TH1F* h_spreadXZ[10];
  TH1F* h_spreadYZ[10];
  for (int i=0;i<10;i++){
    h_spreadX[i] = new TH1F("","X spread",24,0,24);
    h_spreadY[i] = new TH1F("","Y spread",8,0,8);
    h_spreadXZ[i] = new TH1F("","XZ spread",24,0,24);
    h_spreadYZ[i] = new TH1F("","YZ spread",8,0,8);    
  }

  // Filling Histograms_______________________________________________________________________________

  int counter2 = 0;
  int counter3 = 0;
  int counter4 = 0;
  int counter5 = 0;
  //Int_t nNeutrons = 0;
  //Int_t QmaxX = 0;
  //Int_t QmaxY = 0;
  //Int_t QmaxZ = 0;
  //Int_t QmaxNum = 0;
  Int_t prevEventTime = 0;
  //loop over events
  for (size_t iev=evtIni;iev<evtFin; iev++){
    if (iev == maxEvents-1) break;
    cout<<"_Processing events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;
    
    std::vector<ND280SFGDHit*> mppc_hits;
    mppc_hits = getEventMPPChits(iev);

    //skip FEB12 triggers that are within 1.8us of the previous trigger
    //std::cout<<unpackEvent->GetFEB12hitTfromSpill()<<" "<<prevEventTime<<std::endl;
    if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 1400)
      continue;

    counter2 ++;
    //Int_t highChargeHits = 0;

    //Skip event if there are less than 3 hits in the event
    //if (mppc_hits.size()<10) continue; 

    Double_t time = 1e10;
    Double_t energy = -1;
    double timm = 1e9;
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      //time cut 
      counter3++;
      if (hit->GetDt() < -358 || hit->GetDt() > 1040) continue;
      counter4++;
      //if (hit->GetFEB() == 56 || hit->GetFEB() == 57 || hit->GetFEB() == 58 || hit->GetFEB() == 59 || hit->GetFEB() == 60 || hit->GetFEB() == 61)
      if (hit->GetFEB() == 56 && hit->GetCh() < 32)
      {
      if (hit->GetPE() > 20)
      {
	h_tspill -> Fill(hit->GetTfromSpill()); 
	h_dt->Fill(hit->GetDt());
        //cout<<"hit time "<<hit->GetDt()<<endl;
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
      else if (hit->GetFEB()< 56) {
      if (hit->GetPE() > 20){
        h_dtt->Fill(hit->GetDt());
	counter5++;
      }
      }
      if (hit->GetFEB() == 56 && hit->GetCh() > 31 && hit->GetCh()<64 ){  h_dt2->Fill(hit->GetDt()); }
      if (hit->GetFEB() == 56 && hit->GetCh() > 63 ){  h_dt3->Fill(hit->GetDt()); }
      //if (hit->GetFEB() == 59){  h_dt4->Fill(hit->GetDt()); }
      //if (hit->GetFEB() == 60){  h_dt5->Fill(hit->GetDt()); }
      //if (hit->GetFEB() == 61){  h_dt6->Fill(hit->GetDt()); }
    }
    //cout<<"time and energy : "<<time<<" "<<energy<<endl;
    //h_dt->Fill(timm);
    h_energy ->Fill(energy);

    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      //if (hit->GetDt() < -326 || hit->GetDt() > 358)
      if (hit->GetDt() < -358 || hit->GetDt() > 1058)
        continue;
      if (hit->GetPE() < 20) continue;
      if (hit->GetFEB() == 56 || hit->GetFEB() == 57 || hit->GetFEB() == 58 || hit->GetFEB() == 59 || hit->GetFEB() == 60 || hit->GetFEB() == 61){ 
        if (hit->GetView() == 0){
          EventsMapXYUSJ->Fill(hit->GetX(),hit->GetY(),1);
        // if ( hit->GetFEB() == 57 || hit->GetFEB() == 60 || hit->GetFEB() == 61){ 
          // if ((hit->GetFEB() == 61 && hit->GetCh()<32) || (hit->GetFEB() == 60 && hit->GetCh()>63)) 
            // EventsMapXZUSJ->Fill(hit->GetX(),hit->GetZ(),1);
          // else
            // EventsMapZYUSJ->Fill(hit->GetX(),hit->GetZ(),1);
        }
        else if (hit->GetView() == 1){
        // else if (hit->GetFEB() == 58 || hit->GetFEB() == 59){
          EventsMapXZUSJ->Fill(hit->GetX(),hit->GetZ(),1);
        }
        else if (hit->GetView() == 2){
        // else{
          // if (hit->GetCh()>63)
            EventsMapZYUSJ->Fill(hit->GetZ(),hit->GetY(),1);
          // else
            // EventsMapXYUSJ->Fill(hit->GetX(),hit->GetZ(),1);
        }
      }
      else{
        if (hit->GetView() == 0){
          if (hit->GetX()<24 && hit->GetX()>-1 && hit->GetY()<8 && hit->GetY()>-1){
            EventsMapXY->Fill(hit->GetX(),hit->GetY(),1);
	    for(int dtloop = 0;dtloop< 10; dtloop++){
	      if (hit->GetDt()> -358 + dtloop * 70 ){
		h_spreadX[dtloop] ->Fill(hit->GetX());
	        h_spreadY[dtloop] ->Fill(hit->GetY());
	      }
	    }	    
          }
        }
        else if (hit->GetView() == 1){
          if (hit->GetX()<24 && hit->GetX()>-1 && hit->GetZ()<48 && hit->GetZ()>-1){
            EventsMapXZ->Fill(hit->GetX(),hit->GetZ(),1);
            for(int dtloop = 0;dtloop< 10; dtloop++){
              if (hit->GetDt()> -358 + dtloop * 70 && hit->GetZ()> 46 ){
                h_spreadXZ[dtloop] ->Fill(hit->GetX());
              }
            }
          }
        }
        else if (hit->GetView() == 2){
          if (hit->GetZ()<48 && hit->GetZ()>-1 && hit->GetY()<8 && hit->GetY()>-1){
            EventsMapZY->Fill(hit->GetZ(),hit->GetY(),1);
            for(int dtloop = 0;dtloop< 10; dtloop++){
              if (hit->GetDt()> -358 + dtloop * 70 && hit->GetZ()>46 ){
                h_spreadYZ[dtloop] ->Fill(hit->GetY());
              }
            }
          }
        }
      }
      prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
      //cout << hit->GetSpillTime()*2.5 << endl;
    }    //Draw histograms___________________________________________________________________________________

    for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
  }
  cout<<"_Processing events...100% done"<<endl;
  cout<<"event number summary: "<<evtFin<<" "<<counter2<<" "<<counter3<<" "<<counter4<<" "<<counter5<<endl;

  TGraph* gg1 = new TGraph(10);
  TGraph* gg2 = new TGraph(10);
  TGraph* gg3 = new TGraph(10);
  TGraph* gg4 = new TGraph(10);

  for(int i=0;i<10;i++){
    gg1->SetPoint(i, -358+i*70, h_spreadX[i]->GetRMS());
    gg2->SetPoint(i, -358+i*70, h_spreadY[i]->GetRMS());
    gg3->SetPoint(i, -358+i*70, h_spreadXZ[i]->GetRMS());
    gg4->SetPoint(i, -358+i*70, h_spreadYZ[i]->GetRMS());
  }

  TCanvas* c100 = new TCanvas();
  c100->Divide(2,2);
  c100->cd(1);
  gg1->Draw("ALP");
  gg1->SetTitle("X spread on XY plane");
  c100->cd(2);
  gg2->Draw("ALP");
  gg2->SetTitle("Y spread on XY plane");
  c100->cd(3);
  gg3->Draw("ALP");  
  gg3->SetTitle("X spread on XZ plane");
  c100->cd(4);
  gg4->Draw("ALP");
  gg4->SetTitle("Y spread on YZ plane");

  //cout << "Number of neutrons: " << nNeutrons << endl;
  //cout << "Average XY position: (" << QmaxX/QmaxNum << ", " << QmaxY/QmaxNum << ")" << endl;

  // Plotting //Events__________________________________________________//______________________________

  //new TCanvas();
  //h_dt->Draw();

  new TCanvas();
  h_dtt->Draw();

  h_dt2->Scale(h_dt->Integral()/h_dt2->Integral());
  h_dt3->Scale(h_dt->Integral()/h_dt3->Integral());
/*
  h_dt4->Scale(h_dt->Integral()/h_dt4->Integral());
  h_dt5->Scale(h_dt->Integral()/h_dt5->Integral());
  h_dt6->Scale(h_dt->Integral()/h_dt6->Integral());
*/
  new TCanvas();
  h_dt->SetLineColor(1);
  h_dt->SetLineWidth(1);
  h_dt2->SetLineColor(2);
  h_dt2->SetLineWidth(1);  
  h_dt3->SetLineColor(3);
  h_dt3->SetLineWidth(1);  
  //h_dt4->SetLineColor(4);
  //h_dt4->SetLineWidth(1);
  //h_dt5->SetLineColor(5);
  //h_dt5->SetLineWidth(1);
  //h_dt6->SetLineColor(6);
  //h_dt6->SetLineWidth(1);  
  h_dt2->Draw("");
  h_dt->Draw("same");
  h_dt3->Draw("same");
/*
  new TCanvas();
  h_dt->Draw("hist");
  h_dt3->Draw("same");

  new TCanvas();
  h_dt->Draw("hist");
  h_dt4->Draw("same");
  
  new TCanvas();
  h_dt->Draw("hist");
  h_dt5->Draw("same");

  new TCanvas();
  h_dt->Draw("hist");
  h_dt6->Draw("same");
*/  
  new TCanvas();
  h_tspill->Draw();

  new TCanvas();
  h_energy->Draw();


  gStyle->SetOptStat(0);
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
  EventsMapXYUSJ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapXYUSJ->GetZaxis()->SetTitleOffset(1.5);
  cc->cd(5);
  EventsMapXZUSJ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapXZUSJ->GetZaxis()->SetTitleOffset(1.5);
  cc->cd(6);
  EventsMapZYUSJ->Draw("COLZ");
  gPad->SetRightMargin(0.15);
  EventsMapZYUSJ->GetZaxis()->SetTitleOffset(1.5);

  cc->Update();
  cc->WaitPrimitive();

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
