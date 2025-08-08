#define THIS_NAME TimeOffset
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS

#include "../src/tools/global_header.hh"

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool    batch          = false;
bool    IsMC           = true; 
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "~/Documents/SFGD/sfgd_framework/analysis/bin/wilfs/TimeOffset.root";
TString fileIn         = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
bool IsPrototype        = false;
bool IsReversed         = false;
bool IsTwenty          = false;
int SFGD_X             = 204;
int SFGD_Y             = 56;
int SFGD_Z             = 188;
bool useClean          = false;
bool useNN             = false;
int dTimeIni           = -100;
int dTimeFin           = -100;
int evtIni             = 0;
int evtFin             = 0;
Int_t chargeMin = 25;
Double_t fibreSpeed = 1.7e8; //m/s
Double_t attenLen = 4.5;

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
#include "TF1.h"
#include "TFitResult.h"

void TimeOffset() 
{
  //standard way to handle arguments for SFGD executables
  parseArguments();
  //set up TFiles FileInput and FileOutput and link to TTrees dataOut and data
  linkFilesToTTrees();

  // define FEB numbers without gaps
  Int_t FEBs[28] = {0,1,2,3,4,0,0,0,5,6,7,8,0,0,0,0,9,10,11,12,13,0,0,0,14,15,16,17};
 
  //initialise histograms to store amplitude and hit time for each channel
  Int_t FEB[18] = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};
  TH2F *AmpTime[18*96];
  for (Int_t i=0; i<18; i++) {
    for (Int_t j=0; j<96; j++) {
      AmpTime[96*i+j] = new TH2F("","", 300,0,300, 40,150,250);
      AmpTime[96*i+j]->SetTitle(("FEB"+to_string(FEB[i])+"_channel"+to_string(j)+"; Hit Charge [p.e.]; Hit Time before Trigger [ns]").c_str());
    }
  }

  std::vector<ND280SFGDHit*> mppc_hits;

  for (int iev=0; iev<nEvents; iev++){
    if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;
    cout << " Processing events..." << (int)(iev*100/(Double_t)nEvents) << "% done\r" << flush;
    //returns a vector of ND280SFGDHit objects that contain hit properties.
    //Skips first 20 events and ignores all hits prior to one with a negative
    //charge
    mppc_hits = getEventMPPChits(iev); 
    
    Double_t previousTime = 0;
    Double_t triggerTime = 0;
    size_t nHits = mppc_hits.size();
    //cout << nHits << endl;
    // iterate through hits of a single event
    for (size_t ihit=0; ihit<nHits; ihit++){
      Double_t timeSinceTrigger = mppc_hits[ihit]->GetDt()*2.5; //convert to ns
      Double_t hitTime = mppc_hits[ihit]->GetTfromSpill()*2.5;
      triggerTime = hitTime - timeSinceTrigger;
      Int_t hitCharge = mppc_hits[ihit]->GetPE();
      Int_t HGCharge = mppc_hits[ihit]->GetHG_pe();
      Int_t hitFEB = mppc_hits[ihit]->GetFEB();
      Int_t hitView = mppc_hits[ihit]->GetView();
      Int_t hitZ = mppc_hits[ihit]->GetZ();
      Int_t chNum = mppc_hits[ihit]->GetCh();
      Int_t ToT = mppc_hits[ihit]->GetToT();
      Double_t travelDist = 0;
      Double_t calibTime = 0;
      Int_t pairCharge = 0;
      bool keepHit = false;
      //cout << "test " << ihit << endl;

      // find distance light travelled to reach MPPC
      if (hitView == 2){ // side channels
        for (size_t ihit2=0; ihit2<nHits; ihit2++){ //find hit in other view with same z value
          if (mppc_hits[ihit2]->GetView() == 2 || mppc_hits[ihit2]->GetView() == 0) continue;
          if (hitZ == mppc_hits[ihit2]->GetZ()){
            if (mppc_hits[ihit2]->GetPE() > pairCharge){
              pairCharge = mppc_hits[ihit2]->GetPE();
              keepHit=true;
              if (hitFEB == 1 || hitFEB == 24) // right channels
                travelDist = (24.-mppc_hits[ihit2]->GetX()-0.5)*0.01; // meters
              else // left channels
                travelDist = (mppc_hits[ihit2]->GetX()+0.5)*0.01; // meters
            }
          }
        }
      }
      else if (hitView == 1) { // top/bottom channels
        for (size_t ihit2=0; ihit2<nHits; ihit2++){
          if (mppc_hits[ihit2]->GetView() == 1 || mppc_hits[ihit2]->GetView() == 0) continue;
          if (hitZ == mppc_hits[ihit2]->GetZ()){
            if (mppc_hits[ihit2]->GetPE() > pairCharge){
              pairCharge = mppc_hits[ihit2]->GetPE();
              keepHit=true;
              if (hitFEB == 3 || hitFEB == 4 || hitFEB == 8 || hitFEB == 18 || hitFEB == 19 || hitFEB == 20 ) // top channels
                travelDist = (8.-mppc_hits[ihit2]->GetY()+0.5)*0.01; // meters
              else // bottom channels
                travelDist = (mppc_hits[ihit2]->GetY()+0.5)*0.01; // meters
            }
          }
        }
      }
      else keepHit=true;
      //cout << "dist: " << travelDist << endl;
      // account for light attenuation with hit charge
      Int_t hitChargeCalib = hitCharge*exp(-travelDist/attenLen);
      // calibrate hit time to account for fibre travel time
      calibTime = hitTime-travelDist/fibreSpeed*1e9;
      //cout << "time " << triggerTime << endl;

      // ensure time since previous trigger, or the trigger before that, isn't short
      if (abs(triggerTime - previousTime) > 100){
        //cout << "Passed" << endl;
        // ensure time length of signal isn't too long
        if (ToT < 60){
          if (keepHit){
            if (hitCharge > 0 && hitCharge < 10000){
              //store amplitude and corrected hit time in histogram for the
              //channel
              AmpTime[96*FEBs[hitFEB]+chNum]->Fill(hitCharge, triggerTime-calibTime, 1);
            }
          }
        }
      }
    }
    previousTime = triggerTime; //trigger time from spill for the last processed event

    for (int it=0; it<mppc_hits.size(); it++)
      delete mppc_hits[it];
  } // iterate through events

  // initiate canvas to plot individual charge amp plots
  TCanvas *c1 = new TCanvas("c1","c1",4000,3000);

  //----Fit the charge amplitude histograms and save the fitted parameters in a
  //root file----
  TDirectory *AmpVSTime = FileOutput->mkdir("AmpVSTime");
  // Function to use in time vs amp fit
  Float_t TimeOffsets[5*1728];
  dataOut->Branch("TimeOffsets",&TimeOffsets,"TimeOffsets[8640]/F");
  AmpVSTime->cd();
  TF1 *f1 = new TF1("f1", "[0]-[1]*exp(-[2]*x)-[3]/(x+[4])",0,300);
  TF1 *f2 = new TF1("f2", "[0]-[1]/(-x+[2])",25,300);
  f1->SetParLimits(0,190,205);
  f2->SetParLimits(0,190,205);
  f1->SetParLimits(1,-100,50);
  f2->SetParLimits(1,-1000,0);
  f2->SetParLimits(2,-300,0);
  //f1->SetParLimits(2,0,1);
  //f1->SetParLimits(3,0,1000);
  //f1->SetParLimits(3,190,205);
  //TF1 *f1 = new TF1("f1", "[0]-[1]*exp(-[2]*x)",0,300);
  for (size_t i = 0; i < 1728; i++){
    c1->Clear();
    if ((i > 479 && i < 576) || (i > 767 && i < 864)){ // different fit for type 3 MPPCs
      if (AmpTime[i]->GetEntries() > 100){
        //cout << AmpTime[i]->GetEntries() << endl;
        TFitResultPtr r = AmpTime[i]->Fit("f2","QSR");
        TimeOffsets[5*i+0] = r->Parameter(0);
        TimeOffsets[5*i+1] = 0;
        TimeOffsets[5*i+2] = 0;
        TimeOffsets[5*i+3] = r->Parameter(1);
        TimeOffsets[5*i+4] = r->Parameter(2);
      }
      else {
        for (int j = 0; j < 5; j++){
          TimeOffsets[5*i+j] = 0;
          //cout << TimeOffsets[3*i+j] << endl;
        }
      }
    }
    else {
      if (AmpTime[i]->GetEntries() > 100){
        //cout << "Passed" << endl;
        //cout << AmpTime[i]->GetEntries() << endl;
        TFitResultPtr r = AmpTime[i]->Fit("f1","QSR");
        for (int j = 0; j < 5; j++){
          TimeOffsets[5*i+j] = r->Parameter(j);
          //cout << TimeOffsets[3*i+j] << endl;
        }
      }
      else {
        for (int j = 0; j < 5; j++){
          TimeOffsets[5*i+j] = 0;
          //cout << TimeOffsets[3*i+j] << endl;
        }
      }
    }
    //AmpTime[i]->Draw();
    c1->Write();
  }
  FileOutput->cd();
  dataOut->Fill();
  dataOut->Print();
  dataOut->Write("", TObject::kOverwrite);
  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
