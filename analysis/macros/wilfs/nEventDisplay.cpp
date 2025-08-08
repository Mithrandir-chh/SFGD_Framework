#define THIS_NAME nEventDisplay
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS

#include "../../src/tools/global_header.hh"
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
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "/media/disk_b/standard_software/sfgd_framework/cosmic_data_2025/che_special/005/00/MCR0_Run_005_SubRun_00_2025_06_26_11_32_05__display.root";
TString fileIn         = "/media/disk_b/standard_software/sfgd_framework/cosmic_data_2025/che_special/005/00/MCR0_Run_005_SubRun_00_2025_06_26_11_32_05__events.root";
bool IsPrototype       = true;
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
TTree*          dataIn;
ND280SFGDEvent* inputEvent;
ND280SFGDEvent* recoEvent;
TBranch*        recoBranch; 
TBranch*        inputBranch;
int             nEvents;
Event*          unpackEvent;
Hit*            sampleHit;

#include "../src/tools/global_tools.cc" 
#include "../src/tools/reconstruction_methods.hh" 

void parseArgument(){
    for (int iarg=0; iarg<gApplication->Argc(); iarg++){

        if (string( gApplication->Argv(iarg))=="-h" || string( gApplication->Argv(iarg))=="--help" ){
            cout << "**************************************" << endl;
            cout << "Macros run options:" << endl;
            cout << "   -h || --help      Print help info." << endl;
            cout << "   -d || --data      Allows to read beamtest data input." << endl;
            cout << "   -r || --reverse   Swaps the front and back of the detector, as though it had been rotated 180 degrees about." << endl;
            cout << "   -i || --input     Input file, including the path(*.root)." << endl;
            cout << "   -o || --output    Output file, including the path." << endl;
            cout << "   -b || --batch     Do not shows canvases and closes program at end of execution." << endl;
            cout << "   -f || --flag      Use crosstalk flag to reconstruct voxels." << endl;
            cout << "   -a || --evtIni    Initial event" << endl;
            cout << "   -z || --evtFin    Final event" << endl;
            cout << "   -t || --twenty    Beam data from 20 m location" << endl;
            cout << "Display options: " << endl;
            cout << "   -strue            Shows true event displays." << endl;
            cout << "   -sreco            Shows reco event displays." << endl;
            cout << "**************************************" << endl;
            exit(0);
        }
        else if (string( gApplication->Argv(iarg))=="-d" || string( gApplication->Argv(iarg))=="--data" ){
            IsMC = false;
        }
        else if (string(gApplication->Argv(iarg)) == "-r" || string(gApplication->Argv(iarg)) == "--reverse")
        {
            IsReversed = true;
        }
        else if (string( gApplication->Argv(iarg))=="-i" || string( gApplication->Argv(iarg))=="--input" ){
            iarg++;
            fileIn = gApplication->Argv(iarg);
            if (fileIn.Contains(".root")) {
                cout << "Input file: " << fileIn << endl;
            }
            else {
                cerr << "input file must be a .root file!" << endl;
                exit(1);
            }
        }
        else if (string( gApplication->Argv(iarg))=="-o" || string( gApplication->Argv(iarg))=="--output" ){
            iarg++;
            fileOut = gApplication->Argv(iarg);
        }
        else if (string( gApplication->Argv(iarg))=="-b" || string( gApplication->Argv(iarg))=="--batch" ){
            batch = true;
            gROOT->SetBatch(kTRUE);
        }
        else if (string( gApplication->Argv(iarg))=="-a" || string( gApplication->Argv(iarg))=="--evtIni" ){
            iarg++;
            evtIni = atoi(gApplication->Argv(iarg));
        }
        else if (string( gApplication->Argv(iarg))=="-z" || string( gApplication->Argv(iarg))=="--evtFin" ){
            iarg++;
            evtFin = atoi(gApplication->Argv(iarg));
        }
        else if (string( gApplication->Argv(iarg))=="-t" || string( gApplication->Argv(iarg))=="--twenty" ){
            IsTwenty = true;
        }
        else if (string( gApplication->Argv(iarg))=="-strue"){
            SHOW_TRUE = true;
        }
        else if (string( gApplication->Argv(iarg))=="-sreco"){
            SHOW_RECO = true;
        }
  }
  if(IsPrototype){
    SFGD_X             = 24;
    SFGD_Y             = 8;
    SFGD_Z             = 48;
  }
}


void linkFilesToTTree(){

    FileOutput = new TFile(fileOut.Data(),"RECREATE");
    cout << "s1 works" << endl; exit(1);
    if(!FileOutput) {cout << "no output file! exit." << endl; exit(1);}
    
    dataOut  = new TTree("Events","Events");

    recoEvent = new ND280SFGDEvent();
    recoBranch = dataOut->Branch("recoEvent", "recoEvent", recoEvent);
    recoBranch->Reset();
    recoBranch->SetAddress(&recoEvent);

    FileInput  = new TFile(fileIn.Data(),"update");
    if(!FileInput)  {cout << "no input file! exit." << endl; exit(1);}

    dataIn = (TTree*) FileInput->Get("AllEvents");
    nEvents = dataIn->GetEntries();

    inputBranch = dataIn->GetBranch("Event");
    inputEvent = new ND280SFGDEvent();
    unpackEvent = new Event();

    if(IsMC){
        //cout << "MC FILE" << endl;
        inputBranch->SetAddress(&inputEvent);
    }
    else{
        //cout << "DATA FILE" << endl;
        inputBranch->SetAddress(&unpackEvent);
    }

  if(!evtFin) evtFin = nEvents;
  if(evtFin > nEvents) evtFin = nEvents;
  if(!nEvents){
    cerr << "Input file is empty!" << endl;
    exit(1);
  }
  if(evtIni < 0) {
    cerr << "evtIni can not be negative!" << endl;
    exit(1);
  }
  if(evtFin <= evtIni){
    cerr << "evtFin must be larger than event Ini!" << endl;
  }

}

void nEventDisplay() {

  // read in file names etc. defined on command line
  parseArgument();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent, 
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTree();

  std::vector<ND280SFGDHit*> mppc_hits;

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

  dataIn->Print();
  cout<<"total events = "<<nEvents<<endl;
  cout<<"displaying "<<evtFin-evtIni<<" events"<<endl;

  // Filling Histograms_______________________________________________________________________________

  TH2F *hXY = new TH2F("XY view", "XY view; X position [cm];Y position [cm]; Hit charge [pe]", 24, 0, 24, 8, 0, 8);
  TH2F *hXZ = new TH2F("XZ view", "XZ view; X position [cm];Z position [cm]; Hit charge [pe]", 24, 0, 24, 48, 0, 48);
  TH2F *hZY = new TH2F("ZY view", "ZY view; Z position [cm];Y position [cm]; Hit charge [pe]", 48, 0, 48, 8, 0, 8);

/*
  TH2F *hXY = new TH2F("XY view", "XY view; X position [cm];Y position [cm]; Hit charge [pe]", 24, -12, 12, 24, -12, 12);
  TH2F *hXZ = new TH2F("XZ view", "XZ view; X position [cm];Z position [cm]; Hit charge [pe]", 24, -12, 12, 100, -50, 50);
  TH2F *hZY = new TH2F("ZY view", "ZY view; Z position [cm];Y position [cm]; Hit charge [pe]", 100, -50, 50, 24, -12, 12);
*/ 
/*
  TH2F *hXY = new TH2F("XY view", "XY view; X position [cm];Y position [cm]; Hit charge [pe]", 8, 0, 8, 8, 0, 8);
  TH2F *hXZ = new TH2F("XZ view", "XZ view; X position [cm];Z position [cm]; Hit charge [pe]", 8, 0, 8, 32, 0, 32);
  TH2F *hZY = new TH2F("ZY view", "ZY view; Z position [cm];Y position [cm]; Hit charge [pe]", 32, 0, 32, 8, 0, 8);
*/
  TH1F* h_1D1 = new TH1F("Z dist","Z dist; Z position [cm]; Counts",48,0,48);

  Int_t nNeutrons = 0;
  Int_t QmaxX = 0;
  Int_t QmaxY = 0;
  Int_t QmaxZ = 0;
  Int_t QmaxNum = 0;
  Int_t eventsMissed = 0;
  Int_t prevEventTime = 0;
  gStyle->SetOptStat(0);

  cout<<"start the event loop with event number "<<evtFin<<endl;

  //loop over events
  for (int iev=evtIni; iev<evtFin; iev++){
    if (iev == maxEvents-1) break;

    //cout<<"_Plotting events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;
    cout<<"getting event "<<iev<<endl;

    //mppc_hits = getEventMPPChits(iev);


    dataIn->GetEntry(iev);
    bool jumpNextEvent = false;

    cout<<"got entry "<<iev<<endl;
    vector <ND280SFGDHit*> listOfHits;

        TClonesArray * unpackHits = unpackEvent->GetHits();

	cout<<"unpackHits "<<endl;
        for(Int_t ihit=0; ihit<unpackEvent->GetNHits(); ihit++)
        {
            Hit *hit = (Hit*) unpackHits->At(ihit);
            //if(hit->GetDt() < -100 || hit->GetDt() > -60) continue;
            //if(hit->GetPE() <= 0 ) {jumpNextEvent = true; break;}
            ND280SFGDHit* sfgdHit = new ND280SFGDHit();
            sfgdHit->SetX(hit->GetX());
            sfgdHit->SetY(hit->GetY());
            sfgdHit->SetZ(hit->GetZ());
	    cout<< "location "<< ihit<<" x, y, z "<<hit->GetX()<<" "<<hit->GetY()<<" "<<hit->GetZ()<<" "<<endl;
            if (IsReversed)
            {
                {
                if (hit->GetView() == 0)
                    sfgdHit->SetX(SFGD_X-1 - hit->GetX());
                }
                if (hit->GetView() == 1)
                {
                    sfgdHit->SetX(SFGD_X-1 - hit->GetX());
                    sfgdHit->SetZ(SFGD_Z-1 - hit->GetZ());
                }
                if (hit->GetView() == 2)
                {
                    sfgdHit->SetZ(SFGD_Z-1 - hit->GetZ());
                }
            }
            sfgdHit->SetDt(hit->GetDt());
            sfgdHit->SetView(hit->GetView());
            sfgdHit->SetMultiplicity(0);
            sfgdHit->SetHG_pe(hit->GetHG_pe());
            sfgdHit->SetLG_pe(hit->GetLG_pe());
            sfgdHit->SetToT_pe(hit->GetToT_pe());
            sfgdHit->SetPE(hit->GetPE());
            sfgdHit->SetTfromSpill(hit->GetTfromSpill());
            sfgdHit->SetRE(hit->GetRE());
            sfgdHit->SetFE(hit->GetFE());
            sfgdHit->SetSpillTag(hit->GetSpillTag());
            sfgdHit->SetGTrigTime(hit->GetGTrigTime());
            sfgdHit->SetGTrigTag(hit->GetGTrigTag());
            sfgdHit->SetSpillTime(hit->GetSpillTime());
            sfgdHit->SetSpillTrailTime(hit->GetSpillTrailTime());
            sfgdHit->SetTfromSpill(hit->GetTfromSpill());
            sfgdHit->SetHitAmpl(hit->GetHG_ADC());
            sfgdHit->SetHitLGAmpl(hit->GetLG_ADC());
            sfgdHit->SetFEB(hit->GetFEB());
            sfgdHit->SetCh(hit->GetCh());
            sfgdHit->SetToT(hit->GetToT());
            sfgdHit->SetTrueXTalk(hit->GetTrueXTalk());
            sfgdHit->SetRecoHit(hit->GetRecoHit());
            listOfHits.push_back(sfgdHit);
            //delete sfgdHit;
            //delete hit;
        }



    cout<<"MPPC hits done"<<endl;

    //cout << unpackEvent->GetFEB12LeadTime() << endl;
    //skip FEB12 triggers that are within 1.8us of the previous trigger
    //if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 700)
    //  continue;

    Int_t Qmax0 = 0;
    Int_t Qmax1 = 0;
    Int_t Qmax2 = 0;
    Int_t QmaxEl0[2] = {0, 0};
    Int_t QmaxEl1[2] = {0, 0};
    Int_t QmaxEl2[2] = {0, 0};
    Int_t highChargeHits = 0;

    //Skip event if there are less than 5 hits in the event
    if (mppc_hits.size()<3) 
    {
      eventsMissed += 1;
      continue;
    }
    //cout << "mppc hits = " << mppc_hits.size() << endl; 


    Double_t time = 1e10;
    Double_t energy = -1;
    Double_t l = 90;
    Double_t m = 939.5654133;
    Double_t c = 299792458; //speed of light in m/s
    Double_t gammaToF = l/c;
    Double_t offset = 352;

    int countTest1 = 0;
    int countTest2 = 0;
    int countTest3 = 0;
    

    cout<<"parameter setting done"<<endl;
    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){

      cout<<"mppc hit number "<<ihit<<endl;
      ND280SFGDHit *hit = mppc_hits[ihit];
      //time cut 
      //if (hit->GetDt() < -326 || hit->GetDt() > 360) continue;
      if (hit->GetPE() > 5)
      {

        Double_t ToF = (hit->GetDt()+offset)*2.5*1e-9+gammaToF;
        Double_t v = l/(ToF); //neutron velocity using time of flight and travel distance
        // use the ToF and energy of the earliest high charge hit
        if (time > ToF)
        {
          energy = m/sqrt(1-pow(v/c,2)) - m;
          time = ToF;
        }
	cout<<"ToF, v, energy : "<<ToF<<" "<<v<<" "<<energy<<endl;
      }
    }  

    cout<<"fRange of event "<<iev<<" is : "<<unpackEvent->GetRange()<<endl;
    //if (unpackEvent->GetRange()!=1) continue;
    //loop over hits
    bool inXY = false;
    bool inXZ = false;
    bool inYZ = false;
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      ND280SFGDHit *hit = mppc_hits[ihit];
      //time cut 
      //if (hit->GetDt() < -326 || hit->GetDt() > 360) continue;

      //cut on charge reduces crosstalk hits
      //if (hit->GetPE() > 15 || hit->GetPE() == 0)
      //if (sampleHit->GetTrueXTalk() == false)
      //cout<<sampleHit->GetTrueXTalk() <<endl;
      if(hit->GetPE() > 20)
      {   
        highChargeHits++;
        if (hit->GetView() == 0){
	  if(hit->GetPE() == 0) hXY->Fill(hit->GetX(),hit->GetY(),1);	
	  //else if(hit->GetPE() > 15) 
	  cout<<"XY hit : "<<hit->GetX()<<" "<<hit->GetY()<<endl;
	  inXY = true;
	  hXY->Fill(hit->GetX(),hit->GetY(),hit->GetPE());
          if (hit->GetPE() > Qmax0){
             Qmax0 = hit->GetPE();
             QmaxEl0[0] = hit->GetX();
             QmaxEl0[1] = hit->GetY();
          }
        }
        if (hit->GetView() == 1){
          if(hit->GetPE() == 0) hXY->Fill(hit->GetX(),hit->GetZ(),1);		
	  //else if(hit->GetPE() > 15) 
	  cout<<"XZ hit : "<<hit->GetX()<<" "<<hit->GetZ()<<endl;
          hXZ->Fill(hit->GetX(), hit->GetZ(), hit->GetPE());
	  inXZ = true;
          if (hit->GetPE() > Qmax1){
            Qmax1 = hit->GetPE();
            QmaxEl1[0] = hit->GetX();
            QmaxEl1[1] = hit->GetZ();
          }
        }
        else if(hit->GetView() == 2) {
          if(hit->GetPE() == 0) hZY->Fill(hit->GetZ(),hit->GetY(),1);		
	  //else if (hit->GetPE() > 15) 
	  cout<<"ZY hit : "<<hit->GetZ()<<" "<<hit->GetY()<<endl;
	  hZY->Fill(hit->GetZ(), hit->GetY(), hit->GetPE());
	  inYZ = true;
	  if(hit-> GetY()<6 && hit->GetY()>3)
	    h_1D1->Fill(hit->GetZ(), hit->GetPE());
          if (hit->GetPE() > Qmax2)
          {
            Qmax2 = hit->GetPE();
            QmaxEl2[0] = hit->GetZ();
            QmaxEl2[1] = hit->GetY();
          }
        }
	countTest1 ++;
	if(hit->GetTrueXTalk() ) countTest2 ++;
        if(hit->GetRecoHit() ) countTest3 ++;
        //cout << hit->GetSpillTime()*2.5 << endl;
      }
      prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
    }

    cout<<"hey totoal hit "<<countTest1<<" xtalk "<<countTest2<<" .. and RecoHit "<<countTest3<<endl;

    if (Qmax0 > 0) {
      QmaxX += QmaxEl0[0];
      QmaxY += QmaxEl0[1];
      QmaxNum += 1;
    }

    // Plotting //Events__________________________________________________//______________________________

    if (highChargeHits < 3) {
      eventsMissed += 1;
      continue;
    }
    cout << "Micropulses missed: " << eventsMissed << endl;
    eventsMissed = 0;

    //new TCanvas();
    //h_1D1->Draw();

    cout<<"event energy : "<<energy<<endl;
    cout<<"event momentum : "<<unpackEvent->GetdEdz(0)<<" "<<unpackEvent->GetdEdz(1)<<" "<<unpackEvent->GetdEdz(2)<<endl;

    if(!inXY || !inYZ || !inXZ) continue;
    TCanvas *cc = new TCanvas("cc", "cc", 800, 800);
    cc->Divide(2,2);

    cc->cd(1);
    hXY->Draw("COLZ");
    gPad->SetRightMargin(0.15);
    hXY->GetZaxis()->SetTitleOffset(1.5);
    cc->cd(2);
    hXZ->Draw("COLZ");
    hXZ->GetZaxis()->SetTitleOffset(1.5);
    gPad->SetRightMargin(0.15);
    cc->cd(3);
    hZY->Draw("COLZ");
    hZY->GetZaxis()->SetTitleOffset(1.5);
    gPad->SetRightMargin(0.15);
    cc->Update();
    cc->WaitPrimitive();

    hXY->Reset();
    hXZ->Reset();
    hZY->Reset();
    cout<<"saving event "<<iev<<endl;
    hXY->Write(Form("event_XY_%s",iev));
    hXZ->Write(Form("event_XZ_%s",iev));
    hZY->Write(Form("event_ZY_%s",iev));
    cc->SaveAs(("event_"+to_string(iev)+".png").c_str());
    for (size_t i=0; i<mppc_hits.size();i++) delete mppc_hits[i];
  }
  cout<<"_Plotting events...100% done"<<endl;

  cout << "Number of neutrons: " << nNeutrons << endl;
  cout << "Average XY position: (" << QmaxX/QmaxNum << ", " << QmaxY/QmaxNum << ")" << endl;

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
