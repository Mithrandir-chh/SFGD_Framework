#define THIS_NAME nEnergySpectrum
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

double offsetZ = 50;

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

void nEnergySpectrum()
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

  Double_t m = 939.5654133; //neutron mass in MeV
  Double_t c = 299792458; //speed of light in m/s
  Int_t prevEventTime = -999;
  Double_t offset = 0; //352; //time between FEB12 signal and the center of the gamma flash peak in units of 2.5 ns
  Double_t l = 90; //distance travelled by neutron
  // Int_t nbins = 10000;
  //neutron time window at 90m
  Int_t lowT = -330;
  Int_t highT = 340;
  //change variables if data is from 20m location (must be specified in arguments)
  if (IsTwenty){
    lowT = -281;
    highT = 418;
    l = 20;
    offset = 285;
  }
  Double_t gammaToF = l/c; //photon travel time in secs

  //calculate variable bins  
  int nbins = 10000;
  Double_t xbins[nbins+1];
  for (size_t ibin=nbins+1; ibin+1>0; ibin--)
  {
    Double_t x = l/(ibin*2.5*1e-9+gammaToF);
    xbins[nbins+1-ibin] = m/sqrt(1-pow(x/c,2)) - m;
  }
  // Int_t nbins = highT-lowT;
  // Double_t xbins[nbins+1];
  // Double_t x = 0;
  // Double_t xprev = 1e5;
  // for (Int_t ibin=0; ibin<nbins+1; ibin++)
  // {
    // Double_t s = l/((ibin+lowT+offset)*2.5*1e-9+gammaToF);
    // x = m/sqrt(1-pow(s/c,2)) - m;
    // if (xprev-x < 5)
      // x = xprev-5;
    // xbins[nbins-ibin] = x;
    // cout << x << endl;
    // xprev = x;
  // }

  // xprev = -1e5;
  // for (Int_t ibin=0; ibin<nbins; ibin++){
    // if (xbins[ibin] < xprev){
      // cout << xbins[ibin] << endl;
    // }
    // xprev = xbins[ibin];
  // }
  
  TH1 *energySpec = new TH1D("Detected Neutron Energy Spectrum", "Detected Neutron Energy Spectrum;Kinetic Energy [MeV];Number of Entries/MeV", nbins, xbins);
  TH1D *timeSpec = new TH1D("Detected Neutron Time Spectrum", "Detected Neutron Time Distribution;Event Time [2.5 ns];Number of Entries", 200, -400,1500);
  TH1D *timeSpec2 = new TH1D("time from spill", "time from spill ;Event Time from Spill [2.5 ns];Number of Entries", 260e3, 0,260e3);

  TH1D* h_dt = new TH1D("","dt;dt [ns]; Counts",100,0,1000);

  TH2D* h_nhit = new TH2D("","nhit; neutron energy (MeV); layer", 10,0,800,48,0,48);
  TH2D* h_edep = new TH2D("","avg. edep; neutron energy (MeV); layer", 10,0,800,48,0,48);
  TH2D* h_spread = new TH2D("","spread; neutron energy (MeV); layer", 10,0,800,48,0,48);  

  TH2D* h_XZ_1D = new TH2D("","XZ; X; Z", 24,-12,12,100,-50,50);

  TH1D* h_spread_1D = new TH1D("","1D spread; layer; spread",50,0,50);

  TH1D* h_dtd = new TH1D("","dt diff.; dt (tick) ;Counts ",10,0,10);

  TH1D* h_timeLayer_0 = new TH1D("","time view 0; Counts",2000,0,2000);
  TH1D* h_timeLayer_1 = new TH1D("","time view 1; Counts",2000,0,2000);
  TH1D* h_timeLayer_2 = new TH1D("","time view 2; Counts",2000,0,2000);

  TH1F* h_XZ[50][10];
  for (int i=0;i<50;i++){
    for (int j=0;j<10;j++){
      h_XZ[i][j] = new TH1F("","",100,-50,50);
    }
  }
  double nEventEnergy[20]={};
  double nonEmpty = 0;

  gStyle->SetOptStat(0);
  //loop over events

  for (size_t iev = evtIni; iev < evtFin/10 ; iev++)
  {
    if (iev == maxEvents - 1)
      break;
    cout << "_Processing events..." << (int)(iev * 100 / (evtFin - evtIni)) << "% done         \r" << flush;

    //cout<<"_Plotting events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;

    Int_t highChargeHits = 0;
    Int_t highChargeHitsXZ = 0;
    Double_t energy = 0;
    Double_t time = 1e10;

    mppc_hits = getEventMPPChits(iev);

    //cout<<"true energy and flux weight : "<<unpackEvent->GetXStdDev()<<" "<<unpackEvent->GetYStdDev()<<endl;

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

      if (hit->GetFEB() == 56 || hit->GetFEB() == 57 || hit->GetFEB() == 58 || hit->GetFEB() == 59 || hit->GetFEB() == 60 || hit->GetFEB() == 61) continue;

      //time window for neutrons. Different for 90 m and 20 m locations
      //if (hit->GetDt() < lowT || hit->GetDt() > highT)
      //  continue;
      //fiducial cuts
      /*
       if (hit->GetView() == 0) {
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
      */
      //cut on charge reduces crosstalk hits
      if (hit->GetPE() > 20)
      {
        highChargeHits += 1;
        // calculate neutron ToF by adding time diff of neutron hit and gamma peak to gamma ToF. 
        //Double_t ToF = (hit->GetDt()+offset)*2.5*1e-9+gammaToF; //neutron time of flight in seconds

	timeSpec2->Fill(hit->GetTfromSpill());

	if (hit->GetView() == 0)
	  h_timeLayer_0-> Fill(hit->GetDt());
        if (hit->GetView() == 1)
          h_timeLayer_1-> Fill(hit->GetDt());
        if (hit->GetView() == 2)
          h_timeLayer_2-> Fill(hit->GetDt());	

        Double_t ToF = hit->GetDt()*1.;
	h_dt->Fill(ToF);
	Double_t v = l/(ToF*1e-9); //neutron velocity using time of flight and travel distance
        // use the ToF and energy of the earliest high charge hit
        if (time > ToF)
        {
          energy = m/sqrt(1-pow(v/c,2)) - m;
          time = ToF;
	  //cout<<hit->GetDt()<<" "<<time<<" "<<v<<" "<<m<<" "<<c<<" "<<energy<<endl;
        }
	//if(hit->GetView() == 1){
	  //int temppp = hit->GetZ()+50;
	  //cout<<hit->GetZ()<<" "<<hit->GetX()<<endl;
	  //highChargeHitsXZ += 1;
	  //if(hit->GetZ()<-10) {cout<<"exit"<<endl; exit(1);}
	  //cout<<"check 1 "<<temppp<<" "<<hit->GetX()<<endl;
          //if (temppp<50 && temppp>=0 && hit->GetX()<24 && hit->GetX()>-24){
	    //h_XZ[temppp][]->Fill(hit->GetX());
	  //}
	  //cout<<"filled"<<endl;
	//}
      }
      //cout<<"check 1.2"<<endl;
      //prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
    }
    //cout<<"check 1.5"<<endl;
    if (highChargeHits>2){
      energySpec->Fill(energy,1);
      timeSpec->Fill(time,1);
    }

    for (size_t i=0;i<mppc_hits.size();i++) delete mppc_hits[i];
  }
  cout << "_Processing events...100% done" << endl;

  cout << "Event per micropulse : "<<nonEmpty/evtFin<<endl;

  // Plotting //Events__________________________________________________//______________________________
  // Int_t loop = 0;
  // while (loop < 100){
    // for (int i=0; i<418+281; i++){
      // Int_t flatTime = -281+i;
      // Double_t ToF = (flatTime+offset)*2.5*1e-9+gammaToF; //neutron time of flight in seconds
      // Double_t v = l/ToF; //neutron velocity using time of flight and travel distance
      // Double_t energy = m/sqrt(1-pow(v/c,2)) - m;
      // energySpec->Fill(energy);
      // loop++;
    // }
  // }

  new TCanvas();
  timeSpec ->Draw();

  new TCanvas();
  timeSpec2 ->Draw();

  TCanvas *cc = new TCanvas("cc","cc",0,0,1000,600);
  gPad->SetLogy();
  //cc->cd();
  cc->SetRightMargin(0.09);
  cc->SetLeftMargin(0.15);
  energySpec->Scale(1,"width");
  energySpec->Draw("hist");
  energySpec->GetXaxis()->SetRangeUser(0,1000);

  cc->Update();
  cc->SaveAs("energySpectrum.C");
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
