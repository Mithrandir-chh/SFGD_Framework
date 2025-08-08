#define THIS_NAME TimeResMuon
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS

#include "../src/tools/global_header.hh"

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool    batch          = false;
bool    IsMC           = false; 
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "";
TString fileIn         = "";
bool IsPrototype       = false;
bool IsReversed        = false;
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
Double_t sideTimeRef = 0;
Double_t topTimeRef = 0;
Double_t TimeRef = 0;
vector<Double_t> sideTimeRefVec;
vector<Int_t> sideTimeRefZ;
vector<Double_t> topTimeRefVec;
vector<Int_t> topTimeRefZ;
Int_t sTRcount = 0;
Int_t tTRcount = 0;
Int_t NumberEvDis = 100;
Int_t eventsPassed = 0;
string foutPNG;
Int_t chargeMin = 35;
Double_t fibreSpeed = 1.7e8; // meters per sec
Double_t attenLen = 4.5; // meters (only short attenuation used as the short fibre lengths mean the long regime is not entered. Value from arxiv.org/pdf/1511.06225.pdf)
Int_t diffCut = 0;
Int_t peakCut = 0;
Int_t devCut = 0;
Int_t entryCut = 0;
Int_t layersCut = 0;
Int_t hitsCut = 0;
Int_t cloneCut = 0;
Int_t missedLayers = 0;
bool minPeak = 0;
bool largeHitTimeDif = 0;
bool clonedEvent = 0;
Int_t totalHits = 0;

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
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TImage.h"
#include "TLeaf.h"
#include "TColor.h"
#include "TPaletteAxis.h"

// function to retrieve the first part of the input filename to add to the output file
string GetLocation(string str)
{
  int i = str.rfind("_events.root");
  string way = str.substr(0,i);  
  return way;
}

// cuts applied to select muon events. Uses event info from XZ and ZY faces. Returns 0 if event passes cuts
bool muon_cut(TH2D *eventXZ, TH2D *eventZY)
{
  // cut event if the time over threshold was unusually long for one of the hits
  if (largeHitTimeDif == 0){
    // remove duplicates of events that had more than one trigger
    if (clonedEvent == 0){
      // only save event if enough channels fired
      if (eventXZ->GetEntries()>10 && eventZY->GetEntries()>10){
        // cut event if peak energy deposit less than 50 pe
        if (minPeak == 1){
          // only save event if 2 or less z layers have no hit
          if (missedLayers < 3){
            // only save event if std dev in x and y direction for xy plot is less than 1, to ensure the event is just a particle passing straight through and not recoiling
            if (eventXZ->GetStdDev(1)<1 && eventZY->GetStdDev(2)<1){
              return 1;
            }
            else devCut++;
          }
          else layersCut++;
        }
        else peakCut++;
      }
      else entryCut++;
    }
    else cloneCut++;
  }
  else diffCut++;
  return 0;
}

// cuts applied to select stopping proton events. Uses event info from XZ and ZY faces. Returns 0 if event passes cuts
bool proton_cut(TH2D *histXZ, TH2D *histZY)
{
  if (largeHitTimeDif == 0){
    // ignore event if previous trigger was within 100 ticks (250 ns). This corresponds to two triggers for the same event, so we remove the clone.
    if (clonedEvent == 0){
      // only save event if enough channels fired
      if (histXZ->GetEntries()>10 && histZY->GetEntries()>10){
        // cut event if peak energy deposit less than 50 pe
        if (minPeak == 1){
          // only save event if 2 or less z layers have no hit
          if (missedLayers < 3){
            // only save event if std dev in x and y direction for xy plot is less than 1
            if (histXZ->GetStdDev(1)<1 && histZY->GetStdDev(2)<1){
              return 0;
            }
          }
        }
      }
    }
  }
  return 1;
}

struct combo{
  Double_t result;
  Int_t maxFreq;
  vector<Int_t> elem;
};

// Function for calculating median 
combo findMedian(vector<Double_t> a) 
{ 
  combo out;
  out.maxFreq = 0;
  size_t n = a.size();
  if (n == 0){
    cout << "empty vector" << endl;
    out.result = 0;
    out.elem.push_back(-2);
    return out;
  }
  sort(a.begin(), a.end());
  
  // check for odd case 
  if (n % 2 != 0){
     out.result = (double)a[(n-1)/2];
     out.elem.push_back((n-1)/2);
     return out;
  }

  out.result = (double)(a[(n-2)/2] + a[n/2])/2.0; 
  out.elem.push_back(n/2);
  return out;
} 

// Function for calculating median 
combo findMedian(vector<Int_t> a) 
{ 
  combo out;
  out.maxFreq = 0;
  size_t n = a.size();
  if (n == 0){
    cout << "empty vector" << endl;
    out.result = 0;
    out.elem.push_back(-2);
    return out;
  }
  sort(a.begin(), a.end());
  
  // check for odd case 
  if (n % 2 != 0){
     out.result = (double)a[(n-1)/2];
     out.elem.push_back((n-1)/2);
     return out;
  }

  out.result = (double)(a[(n-2)/2] + a[n/2])/2.0; 
  out.elem.push_back(n/2);
  return out;
} 

// Function for finding mode
combo findMode(vector<Double_t> a){
  Int_t localFreq = 1;
  combo out;
  out.maxFreq = 1;
  Int_t maxIter = 0;
  bool multi = 0;
  vector<Int_t> vec = {0};
  size_t n = a.size();
  if (n == 0)
    return (combo){0,0,{-2}};
  sort(a.begin(), a.end());

  for (Int_t i=0; i<(int)n-1; i++){
    if (a[i] == a[i+1]){
      localFreq += 1;
      vec.push_back(i+1);
      if (i==(int)n-2){
        if (localFreq > out.maxFreq){
          out.maxFreq = localFreq;
          out.elem.clear();
          out.elem = vec;
          maxIter = i;
          multi = 0;
        }
        else if (localFreq == out.maxFreq){
          multi = 1;
        }
      }
    }
    else {
      if (localFreq > out.maxFreq){
        out.maxFreq = localFreq;
        out.elem.clear();
        out.elem = vec;
        maxIter = i;
        multi = 0;
      }
      else if (localFreq == out.maxFreq)
        multi = 1;
      vec.clear();
      vec.push_back(i+1);
      localFreq = 1;
    }
  }
  out.result = a[maxIter];
  if (multi){
    cout << "No unique mode" << endl;
    //cout << "Freq = " << out.maxFreq << endl;
    //cout << "size = " << n << endl;
    return (combo){0,out.maxFreq,{-2}};
  }
  return out;
}

void TimeResMuon() 
{
  if (gApplication->Argc() != 5){
    cerr << "\nMust give four arguments.\n" << endl;
    cerr << "'avg': use the average hit time for each event as the time reference." << endl;
    cerr << "'front': use the earliest hit at the front of the detector as a time reference." << endl;
    cerr << "'median': use the median hit time for each event as a time reference." << endl;
    cerr << "'trigger': use the trigger time for the event as a time reference." << endl;
    cerr << "\n'fit': use a Gaussian fit to find time resolution." << endl;
    cerr << "'stat': use the statistical standard deviation to find the time resolution.\n" << endl;
    cerr << "'test': save the results in a test file so other analyses aren't overwritten." << endl;
    cerr << "'notest': create standard analysis file, overwriting previous ones." << endl;
    cerr << "\n'<path to data file>': give the path to the input data file.\n" << endl;
    exit(1);
  }
  string argv1 = gApplication->Argv(1);
  string argv2 = gApplication->Argv(2);
  string argv3 = gApplication->Argv(3);
  string argv4 = gApplication->Argv(4);
  if ((argv1 != "avg" && argv1 != "front" && argv1 != "median" && argv1 != "trigger" && argv1 != "mode") || (argv2 != "fit" && argv2 != "stat") || (argv3 != "test" && argv3 != "notest")){
    cerr << "\nIncorrect argument(s).\n" << endl;
    cerr << "'avg': use the average hit time as the time reference." << endl;
    cerr << "'front': use the earliest hit at the front of the detector as a time reference." << endl;
    cerr << "'median': use the median hit time for each event as a time reference." << endl;
    cerr << "'trigger': use the trigger time for the event as a time reference." << endl;
    cerr << "'mode': use the mode time for each event as a time reference." << endl;
    cerr << "\n'fit': use a Gaussian fit to find time resolution." << endl;
    cerr << "'stat': use the statistical standard deviation to find the time resolution\n" << endl;
    cerr << "'test': save the results in a test file so other analyses aren't overwritten." << endl;
    cerr << "'notest': create standard analysis file, overwriting previous ones." << endl;
    cerr << "\n'<path to data file>': give the path to the input data file.\n" << endl;

    exit(1);
  }

  // define FEB numbers 
  Int_t FEBs[28] = {0,1,2,3,4,0,0,0,5,6,7,8,0,0,0,0,9,10,11,12,13,0,0,0,14,15,16,17};

  // get input and output files
  fileIn = argv4;
  fileOut = GetLocation(argv4);
  if (argv3 == "test")
    fileOut += "_TimeRes_test.root";
  else
    fileOut += "_TimeRes_"+argv1+"_"+argv2+"_new.root";
  cout << "Output: " << fileOut << endl;

  linkFilesToTTrees();

  cout << "Number of events: " << nEvents << endl;

  // initiate plots for channel time res measurements. One plot per channel (18 FEBs with 96 channels each)
  Int_t FEB[18] = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};
  TH1F *Time_res[18*96];
  // retrieve mapping info and put into 3D matrix
  for (int iFEB = 0; iFEB<18; iFEB++) {
    for (int jCh = 0; jCh<96; jCh++) {
      if (argv1 == "trigger")
        Time_res[96*iFEB+jCh] = new TH1F("","",40,230,330);
      else
        Time_res[96*iFEB+jCh] = new TH1F("","",40,-50,50);
      Time_res[96*iFEB+jCh]->SetTitle(("FEB"+to_string(FEB[iFEB])+"_channel"+to_string(jCh)+"_TimeRes; Time Difference [ns]; Counts").c_str());
    }
  }
  
  // initiate our event displays for all events with relevant titles and dimensions
  //TH2F *EventsMap_XY = new TH2F("All_events_map_XY","All_events_map_XY",  24,0,24, 8,0,8);
  //EventsMap_XY->SetTitle("XY View; x position [cm]; y position [cm]");
  //TH2F *EventsMap_ZY = new TH2F("All_events_map_ZY","All_events_map_ZY",  48,0,48, 8,0,8);
  //EventsMap_ZY->SetTitle("ZY View; z position [cm]; y position [cm]");
  //TH2F *EventsMap_XZ = new TH2F("All_events_map_XZ","All_events_map_XZ",  24,0,24, 48,0,48);
  //EventsMap_XZ->SetTitle("XZ View; x position [cm]; z position [cm]");
  TH2F *TimeResMap_ZY = new TH2F("Time_Res_map_ZY","; z position [cm]; y position [cm];",  48,0,48, 8,0,8);
  TimeResMap_ZY->SetMinimum(0.5);
  TimeResMap_ZY->SetMaximum(1.5);
  TH2F *TimeResMap_XZ = new TH2F("Time_Res_map_XZ","; z position [cm]; x position [cm]; Time Resolution [ns]",  48,0,48, 24,0,24);
  TimeResMap_XZ->GetZaxis()->CenterTitle(true);
  TimeResMap_XZ->SetMinimum(0.5);
  TimeResMap_XZ->SetMaximum(1.5);
  
  // create directories in our output file to store our event displays and MIP plots
  //TDirectory *events2D = FileOutput->mkdir("events2D");
  TDirectory *timeRes = FileOutput->mkdir("timeRes");
 
  // initiate canvases to plot individual event displays, MIP peaks and time resolution plots
  if (batch) gROOT->SetBatch(kTRUE);

  // initiate event displays for individual events
  //TH2D *event_XY = new TH2D("","", 24,0,24, 8,0,8);
  //event_XY->SetTitle("; x position [cm]; y position [cm]; Light Yield [p.e.]");
  //event_XY->GetXaxis()->SetTitleSize(0.05);
  //event_XY->GetYaxis()->SetTitleSize(0.05);
  //event_XY->GetZaxis()->SetTitleSize(0.05);
  //event_XY->GetXaxis()->SetLabelSize(0.05);
  //event_XY->GetYaxis()->SetLabelSize(0.05);
  //event_XY->GetZaxis()->SetLabelSize(0.05);
  
  TH2D *event_ZY = new TH2D("","", 48,0,48, 8,0,8);
  event_ZY->SetTitle("; x position [cm]; y position [cm]; Light Yield [p.e.]");
  event_ZY->GetXaxis()->SetTitleSize(0.05);
  event_ZY->GetYaxis()->SetTitleSize(0.05);
  event_ZY->GetZaxis()->SetTitleSize(0.05);
  event_ZY->GetXaxis()->SetLabelSize(0.05);
  event_ZY->GetYaxis()->SetLabelSize(0.05);
  event_ZY->GetZaxis()->SetLabelSize(0.05);
  
  TH2D *event_XZ = new TH2D("","", 24,0,24, 48,0,48);
  event_ZY->SetTitle("; x position [cm]; z position [cm]; Light Yield [p.e.]");
  event_XZ->GetXaxis()->SetTitleSize(0.05);
  event_XZ->GetYaxis()->SetTitleSize(0.05);
  event_XZ->GetZaxis()->SetTitleSize(0.05);
  event_XZ->GetXaxis()->SetLabelSize(0.05);
  event_XZ->GetYaxis()->SetLabelSize(0.05);
  event_XZ->GetZaxis()->SetLabelSize(0.05);

  Double_t energyDep[48];
  Double_t energyDepXZ[48];
  Double_t energyDepZY[48];
  Int_t eventNum = 0;
  Double_t timeSinceTrigger;
  Double_t hitTime;
  Double_t hitToT;
  Int_t hitCharge;
  Double_t hitChargeCalib;
  Int_t hitFEB;
  Int_t hitView;
  Int_t hitZ;
  Int_t hitX;
  Int_t hitY;
  Int_t chNum;
  Double_t travelDist;
  Double_t calibTime;
  Int_t pairCharge;
  bool keepHit;
  Double_t previousTime;

  // retrieve time walk offset parameters for each channel
  TFile *f = new TFile("/home/ws1216/Documents/SFGD/Unpacking_august_v2.2.1/bin/TimeOffset.root", "READ");
  TTree *Offtree = (TTree*)f->Get("tree");
  Float_t TimeOffsets[5*1728];
  Offtree->SetBranchAddress("TimeOffsets",&TimeOffsets[0]);
  Offtree->GetEvent();
  f->Close();
  // for (int j=0; j<1728; j++){
    // cout << endl;
    // for (int i=0; i<5; i++)
      // cout << TimeOffsets[j*5+i] << ",";
  // }

  size_t maxEv = nEvents;
  if (argv3 == "test") maxEv = 1;

  //TCanvas *c1 = new TCanvas("c1","c1",800,600);

//=========Start iterating through events=================

  for (size_t iev = 0; iev<maxEv; iev++) {
    //if(selEvents >= maxSelEvents) break;
    cout << "Processing events..." << (int)(iev*100/(Double_t)maxEv) << "% done   \r" << flush;
    //cout << "Processing events..." << iev << "/" << maxEv << " done        \r" << flush;

    clonedEvent = 0;
    largeHitTimeDif = 0;

    for (size_t i=0; i<48; i++){
      energyDep[i] = 0;
      energyDepXZ[i] = 0;
      energyDepZY[i] = 0;
    }

    std::vector<ND280SFGDHit*> mppc_hits;

    //returns a vector of ND280SFGDHit objects that contain hit properties.
    //Skips first 20 events and ignores events that contain hits with negative
    //charges
    mppc_hits = getEventMPPChits(iev);

    sideTimeRef = 0;
    topTimeRef = 0;
    TimeRef = 0;
    sideTimeRefVec.clear();
    sideTimeRefZ.clear();
    topTimeRefVec.clear();
    topTimeRefZ.clear();
    sTRcount = 0;
    tTRcount = 0;
    Int_t element = 0;

    Double_t triggerTime = 0;
    size_t nHits = mppc_hits.size();
    //if (nHits == 0) {hitsCut++; continue;}
    // iterate through hits of a single event
    for (size_t ihit=0; ihit<nHits; ihit++){
      totalHits++;
      timeSinceTrigger = mppc_hits[ihit]->GetDt()*2.5; //convert to ns
      hitTime = mppc_hits[ihit]->GetTfromSpill()*2.5;
      hitToT = mppc_hits[ihit]->GetToT();
      triggerTime = hitTime - timeSinceTrigger;
      hitCharge = mppc_hits[ihit]->GetPE();
      hitFEB = mppc_hits[ihit]->GetFEB();
      hitView = mppc_hits[ihit]->GetView();
      hitZ = mppc_hits[ihit]->GetZ();
      hitX = mppc_hits[ihit]->GetX();
      hitY = mppc_hits[ihit]->GetY();
      chNum = mppc_hits[ihit]->GetCh();
      travelDist = 0;
      calibTime = 0;
      pairCharge = 0;
      keepHit = false;

      // find distance light travelled to reach MPPC
      if (hitView == 2){ // side channels
        for (size_t ihit2=0; ihit2<nHits; ihit2++){ //find hit in other view with same z value
          if (mppc_hits[ihit2]->GetView() == 1){
            if (hitZ == mppc_hits[ihit2]->GetZ() && mppc_hits[ihit2]->GetPE() > pairCharge){
              pairCharge = mppc_hits[ihit2]->GetPE();
              keepHit = true; // keep hit only if there is a hit in the same z layer of the top/bottom channels
              if (hitFEB == 1 || hitFEB == 24) // right channels
                travelDist = (24.-mppc_hits[ihit2]->GetX()+0.5)*0.01; // meters
              else // left channels
                travelDist = (mppc_hits[ihit2]->GetX()+0.5)*0.01; // meters
            }
          }
        }
      }
      else if (hitView == 1) { // top/bottom channels
        for (size_t ihit2=0; ihit2<nHits; ihit2++){
          if (mppc_hits[ihit2]->GetView() == 2){
            if (hitZ == mppc_hits[ihit2]->GetZ() && mppc_hits[ihit2]->GetPE() > pairCharge){
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
      else keepHit = true;
      // account for light attenuation with hit charge
      hitChargeCalib = hitCharge*exp(-travelDist/attenLen);
      // calibrate hit time to account for time-walk effect, fibre travel time and light attenuation
      element = FEBs[hitFEB]*96*5+chNum*5;
      if (hitFEB == 8 || hitFEB == 11) //type 3 MPPCs have a different fit
        calibTime = hitTime-TimeOffsets[element+1]/(TimeOffsets[element+2]+hitCharge)-(travelDist/fibreSpeed)*1e9;
      else
        calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
      // cut on time window for muons (and electrons and pions)
      if (triggerTime - hitTime < 300){
          //&& triggerTime - calibTime > 250){
        if (hitCharge > chargeMin){
          if (keepHit){
            //skip event if it's trigger time is similar to the previous event's, meaning there was probably two triggers for a single event
            if (abs(triggerTime - previousTime) < 250) clonedEvent = 1;
            // record if time length of signal is long
            if (hitToT > 60) largeHitTimeDif = 1;
            if (hitCharge < 10000){
              // front and back: plot event 
              //if (hitView == 0){
                //event_XY->Fill(hitX, hitY, hitChargeCalib);
              //} 
              // sides: plot event
              if (hitView == 2){
                event_ZY->Fill(hitZ, hitY, hitCharge);
                if (argv1 == "avg"){
                  sideTimeRef += calibTime;
                  sTRcount += 1;
                }
                if (argv1 == "median" || argv1 == "mode"){
                  sideTimeRefVec.push_back(calibTime); // round time to 2 dp so we can find a mode
                  sideTimeRefZ.push_back(hitZ);
                }
                energyDepZY[(int)hitZ] += hitCharge;
              } 
              // top and bottom: plot event
              else if (hitView == 1){
                event_XZ->Fill(hitX, hitZ, hitCharge);

                if (argv1 == "avg"){
                  topTimeRef += calibTime;
                  tTRcount += 1;
                }
                if (argv1 == "median" || argv1 == "mode"){
                  topTimeRefVec.push_back(hitTime); // round time to 2 dp so we can find a mode
                  topTimeRefZ.push_back(hitZ);
                }
                energyDepXZ[(int)hitZ] += hitChargeCalib;
              } 
            } //max charge cut
          } //keep hit
        } //min charge cut
      } //time window
    } //loop through hits
    // fill TH1Fs with summed light yield from z layers
    for (int ij = 0; ij < 48; ij++ ) energyDep[ij] = energyDepZY[ij] + energyDepXZ[ij];
    missedLayers = 0;
    minPeak = 0;
    for (int ik = 0; ik < 48; ik++ ){
      // check if hits recorded in all z layers
      if (energyDep[ik] == 0)
        missedLayers += 1;
      // check if any bins are over 50 pe
      if (energyDep[ik] > 50)
        minPeak = 1;
    }
    // only take events that go straight through the detector with no interactions
    if (muon_cut(event_XZ, event_ZY) == 1){
      eventsPassed++;
      // use the collected times from the processed hits to find the time
      // reference using the method specified by the user
      if (argv1 == "avg"){
        sideTimeRef = sideTimeRef/sTRcount; // time reference is the average hit time of the event
        topTimeRef = topTimeRef/tTRcount;
      }
      if (argv1 == "median"){
        sort(sideTimeRefZ.begin(), sideTimeRefZ.end());
        sort(topTimeRefZ.begin(), topTimeRefZ.end());
        combo median = findMedian(sideTimeRefVec);
        sideTimeRef = median.result;
        //cout << "median:" << endl;
        //cout << median.result << endl;
        //cout << "Element:" << endl;
        //cout << median.elem[0] << endl;
        median = findMedian(topTimeRefVec);
        topTimeRef = median.result;
        //topTimeRefVec.insert(end(topTimeRefVec),begin(sideTimeRefVec),end(sideTimeRefVec));
        //TimeRef = findMedian(topTimeRefVec);
      }
      if (argv1 == "mode"){
        sort(sideTimeRefZ.begin(), sideTimeRefZ.end());
        sort(topTimeRefZ.begin(), topTimeRefZ.end());
        //cout.precision(12);
        //for (int i=0; i<sideTimeRefVec.size(); i++)
        // cout << sideTimeRefVec[i] << ", ";
        //cout << endl;
        //cout << findMode(sideTimeRefVec).maxFreq << endl;
        combo mode = findMode(sideTimeRefVec);
        sideTimeRef = mode.result;
        vector<Int_t> modeLayer;
        for (Int_t k=0; k<(int)mode.elem.size(); k++)
          modeLayer.push_back(sideTimeRefZ[mode.elem[k]]);
        if (mode.elem[0] != -2){
          combo midmode = findMedian(modeLayer);
        }
        mode = findMode(topTimeRefVec);
        topTimeRef = mode.result;
        modeLayer.clear();
        for (Int_t k=0; k<(int)mode.elem.size(); k++)
          modeLayer.push_back(topTimeRefZ[mode.elem[k]]);
        if (mode.elem[0] != -2){
          combo midmode = findMedian(modeLayer);
        }
      }

      //timeRef = FEB[12].hitTimefromSpill->at(TOFtrigger);
   
    // report if no reference time was found for selected event
    //if (timeRef == 0){
    //  cout << "NO REF" << endl;
    //  cout << subSpill+1 << " " << numberEvents+1 << endl;
    //}
      // iterate through hits again and build histograms using the time
      // references 
      for (size_t ihits=0; ihits<nHits; ihits++){
        timeSinceTrigger = mppc_hits[ihits]->GetDt()*2.5; //convert to ns
        hitTime = mppc_hits[ihits]->GetTfromSpill()*2.5;
        triggerTime = hitTime - timeSinceTrigger; 
        hitCharge = mppc_hits[ihits]->GetPE();
        hitFEB = mppc_hits[ihits]->GetFEB();
        hitView = mppc_hits[ihits]->GetView();
        hitZ = mppc_hits[ihits]->GetZ();
        hitX = mppc_hits[ihits]->GetX();
        hitY = mppc_hits[ihits]->GetY();
        chNum = mppc_hits[ihits]->GetCh();
        travelDist = 0;
        calibTime = 0;
        pairCharge = 0;
        keepHit = false;

        // find distance light travelled to reach MPPC
        if (hitView == 2){ // side channels
          for (size_t ihit2=0; ihit2<nHits; ihit2++){
            if (mppc_hits[ihit2]->GetView() == 1){
              if (mppc_hits[ihit2]->GetZ() == hitZ && mppc_hits[ihit2]->GetPE() > pairCharge){
                pairCharge = mppc_hits[ihit2]->GetPE();
                keepHit = true;
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
            if (mppc_hits[ihit2]->GetView() == 2){
              if (hitZ == mppc_hits[ihit2]->GetZ() && mppc_hits[ihit2]->GetPE() > pairCharge){
                pairCharge = mppc_hits[ihit2]->GetPE();
                keepHit = true;
                if (hitFEB == 3 || hitFEB == 4 || hitFEB == 8 || hitFEB == 18 || hitFEB == 19 || hitFEB == 20 ) // top channels
                  travelDist = (8.-mppc_hits[ihit2]->GetY()-0.5)*0.01; // meters
                else // bottom channels
                  travelDist = (mppc_hits[ihit2]->GetY()+0.5)*0.01; // meters
              }
            }
          }
        }
        else keepHit = true;
        // calibrate hit time to account for time-walk effect, fibre travel time and light attenuation
        element = FEBs[hitFEB]*96*5+chNum*5;
        // if (hitFEB == 8 || hitFEB == 11) //type 3 MPPCs have a different fit
          // calibTime = hitTime-TimeOffsets[element+1]/(TimeOffsets[element+2]-hitCharge)-(travelDist/fibreSpeed)*1e9;
        // else
          calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
        //calibTime = hitTime;

        // cut on time window for muons (and electrons and pions)
        if (triggerTime - hitTime < 300
            && triggerTime - hitTime > 230){
          if (keepHit) {
            if (hitCharge > chargeMin && hitCharge < 10000){
              // front and back: plot event 
              //if (hitView == 0){
                // fill map showing all events
                //EventsMap_XY->Fill(hitX, hitY, 1);
              //} 
              // sides: plot event
              if (hitView == 2){
                //EventsMap_ZY->Fill(hitZ, hitY, 1);
                // fill time res plot for relevant channels, but only if a time reference was able to be found and if the hit does not look like cross-talk
                if (argv1 == "trigger")
                  Time_res[FEBs[hitFEB]*96+chNum]->Fill(-calibTime + triggerTime);
                else {
                  if (sideTimeRef != 0)
                    Time_res[FEBs[hitFEB]*96+chNum]->Fill(calibTime - sideTimeRef);
                }
              } 
              // top and bottom: plot event
              else if (hitView == 1){
                //EventsMap_XZ->Fill(hitX, hitZ, 1);
                // fill time res plot for relevant channels, but only if a time reference was able to be found and if the hit does not look like cross-talk
                if (hitFEB == 18) continue;
                if (topTimeRef != 0){
                  if (argv1 == "trigger")
                    Time_res[FEBs[hitFEB]*96+chNum]->Fill(-calibTime + triggerTime);
                  else {
                    if (topTimeRef != 0)
                      Time_res[FEBs[hitFEB]*96+chNum]->Fill(calibTime - topTimeRef);
                  }
                }
              }
            }
          }
        }
      } //iterate through hits
      
//============create event displays================
  
      //const Int_t Number = 3;
      //Double_t Length[Number] = {0., 0.5, 0.};
      //Double_t DBlue[Number] = {0., 0., 0.};
      //Double_t Blue[Number] = {0., 0., 0.};
      //Double_t Orange[Number] = {0., 0., 1.};
      //Int_t nb = 50;
      //if (eventNum<NumberEvDis){       
      //  c1->Clear();
      //  // customise stats box
      //  gROOT->ForceStyle();
      //  gStyle->SetStatX(0.3);
      //  gStyle->SetStatY(1);
      //  gStyle->SetStatH(0.1155555);
      //  gStyle->SetStatW(0.2);
      //  gStyle->SetOptStat("rme");
      //  gStyle->SetTitleY(0.95);
      //  gStyle->SetPadRightMargin(0.2);
      //  gStyle->SetPadTopMargin(0.15);
      //  gStyle->SetFrameLineWidth(2);
      //  TColor::CreateGradientColorTable(Number, Length, DBlue, Blue, Orange, nb);
      //  c1->Update();
      //  event_ZY->SetContour(nb);
      //  event_ZY->Draw("colorz");
      //  c1->Update();
      //  // save the canvas in events2D folder in output file
      //  c1->Write();
      //}
    } // muon cut
    
    //##########
    //##########

    previousTime = triggerTime;
    eventNum++;
    //event_XY->Reset();
    event_ZY->Reset();
    event_XZ->Reset();
    for (int it=0;it<mppc_hits.size();it++) delete mppc_hits[it];
  } // iterate through events
  cout << "Processing events...100% done" << endl;
  // iterate through side channels to find minimum entry number for time res
  cout << "Events passed: " << eventsPassed << endl;
  cout << "hits cut: " << hitsCut << endl;
  cout << "large ToT cut: " << diffCut << endl;
  cout << "clone cut: " << cloneCut << endl;
  cout << "low entries cut: " << entryCut << endl;
  cout << "missed layers cut: " << layersCut << endl;
  cout << "low peak cut: " << peakCut << endl;
  cout << "std dev cut: " << devCut << endl;
  cout << "\ntotal hits: " << totalHits << endl;
  int n = 0;
  int SminEn = 1e6;
  for (int i = 0; i < 18; i++){
    if (FEB[i] == 1 || FEB[i] == 24  || FEB[i] == 2  || FEB[i] == 17){
      for (int ch = 0; ch < 96; ch++){
        n = Time_res[i*96+ch]->GetEntries();
        if (n < SminEn && n != 0) SminEn = n;
      }
    }
  }
  if (SminEn<150)
    SminEn = 150;
  FileOutput->cd();
  //EventsMap_XY->Write();
  //EventsMap_ZY->Write();
  //EventsMap_XZ->Write();
  // create function to use in Time Res fit
  //TF1 *f1 = new TF1("f1", "gaus", -20, 15);
  if (argv1 == "trigger")
    TF1 *f1 = new TF1("f1", "gaus", 260, 330);
  else
    TF1 *f1 = new TF1("f1", "gaus", -20, 15);
  vector<Double_t> sigmas1;
  vector<Double_t> channels1;
  vector<Double_t> xe1;
  vector<Double_t> ye1;
  vector<Double_t> sigmas2;
  vector<Double_t> channels2;
  vector<Double_t> xe2;
  vector<Double_t> ye2;
  vector<Double_t> sigmas3;
  vector<Double_t> channels3;
  vector<Double_t> xe3;
  vector<Double_t> ye3;
  TH1F *TimeResDist = new TH1F("TimeResDist","TimeResDist",200,0,4);
  TimeResDist->SetTitle("; Channel Time Resolution [ns]; Channels");
  TH1F *TimeResDistL = new TH1F("TimeResDistL","TimeResDistL",200,0,4);
  TimeResDistL->SetTitle("; Channel Time Resolution [ns]; Channels");
  TH1F *TimeResDistS = new TH1F("TimeResDistS","TimeResDistS",200,0,4);
  TimeResDistS->SetTitle("; Channel Time Resolution [ns]; Channels");
  TH1F *TimeResDist1 = new TH1F("TimeResDist1","TimeResDist1",200,0,4);
  TimeResDistS->SetTitle("; Channel Time Resolution [ns]; Channels");
  TH1F *TimeResDist2 = new TH1F("TimeResDist2","TimeResDist2",200,0,4);
  TimeResDistS->SetTitle("; Channel Time Resolution [ns]; Channels");
  TH1F *TimeResDist3 = new TH1F("TimeResDist3","TimeResDist3",200,0,4);
  TimeResDistS->SetTitle("; Channel Time Resolution [ns]; Channels");
  //cout << "SminEn = " << SminEn << endl;
  
  // create a 3D matrix to store mapping data for our event displays
  int MapCon[28][2][96];
  for (int iFEB=0; iFEB<18; iFEB++){
      // retrieve a txt file with mapping information for the FEBs
      // txt file has 3 columns: channel no., x position, y position
      string sFEB = "../../../data_preprocessing/mapping/" + to_string(FEB[iFEB]) + ".txt";
      ifstream fmap(sFEB.c_str());
      //cout <<endl<< "FEB "<< FEB[iFEB]<< " mapping"<<endl;
      // read mapping info into 3D matrix
      int temp=0;
      while (!fmap.eof()) {
        fmap >> temp >> MapCon[FEB[iFEB]][0][temp] >> MapCon[FEB[iFEB]][1][temp]; // COLUMNS: temp=(ch. no.), FEB=(FEB no.), 1/0=(vertical/horizontal position) 
        // so MapCon[2][0][26] will be the horizontal position of the cube in FEB_2 at channel 26
      }
      fmap.close();
  }

  TCanvas *c2 = new TCanvas("c2","c2",800,630);
  cout << "Plotting results, this will take a few minutes..." << endl;
  // set t2k style for plots
  TString localStyleName = "T2K";
  Int_t localWhichStyle = 3;
  TStyle *t2kstyle = SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->ForceStyle(t2kstyle);

  for (Int_t ih=96; ih<1728; ih++){
    c2->Clear();
    if (Time_res[ih]->GetEntries() >= 750){ // cut out plots with low stats
      TFitResultPtr r = Time_res[ih]->Fit("f1","QSR");
      if (r == -1) continue;
      int FEBind = (int)(ih/96);
      if (FEB[FEBind] == 3 || FEB[FEBind] == 4 || FEB[FEBind] == 9 || FEB[FEBind] == 10){ //type2 mppc
        if (argv2 == "fit"){
          sigmas2.push_back(r->Parameter(2)/sqrt(2));
          ye2.push_back(r->ParError(2)/sqrt(2));
          xe2.push_back(0.5);
          channels2.push_back(ih);
          TimeResDist->Fill(sigmas2.back());
          TimeResDist2->Fill(sigmas2.back());
          TimeResMap_XZ->Fill(
            MapCon[FEB[FEBind]][1][(int)ih%96], // z position
            MapCon[FEB[FEBind]][0][(int)ih%96], // x position
            sigmas2.back());
        }
        else {
          sigmas2.push_back(Time_res[ih]->GetStdDev()/sqrt(2));
          ye2.push_back(Time_res[ih]->GetStdDevError()/sqrt(2));
          xe2.push_back(0.5);
          channels2.push_back(ih);
          TimeResDist->Fill(sigmas2.back());
          TimeResDist2->Fill(sigmas2.back());
          TimeResMap_XZ->Fill(
            MapCon[FEB[FEBind]][1][(int)ih%96], // z position
            MapCon[FEB[FEBind]][0][(int)ih%96], // x position
            sigmas2.back());
        }
      }
      else if (FEB[FEBind] == 8 || FEB[FEBind] == 11){ //type3 mppc
        if (argv2 == "fit"){
          sigmas3.push_back(r->Parameter(2)/sqrt(2));
          ye3.push_back(r->ParError(2)/sqrt(2));
          xe3.push_back(0.5);
          channels3.push_back(ih);
          TimeResDist->Fill(sigmas3.back());
          TimeResDist3->Fill(sigmas3.back());
          TimeResMap_XZ->Fill(
            MapCon[FEB[FEBind]][1][(int)ih%96], // z position
            MapCon[FEB[FEBind]][0][(int)ih%96], // x position
            sigmas3.back());
        }
        else {
          sigmas3.push_back(Time_res[ih]->GetStdDev()/sqrt(2));
          ye3.push_back(Time_res[ih]->GetStdDevError()/sqrt(2));
          xe3.push_back(0.5);
          channels3.push_back(ih);
          TimeResDist->Fill(sigmas3.back());
          TimeResDist3->Fill(sigmas3.back());
          TimeResMap_XZ->Fill(
            MapCon[FEB[FEBind]][1][(int)ih%96], // z position
            MapCon[FEB[FEBind]][0][(int)ih%96], // x position
            sigmas3.back());
        }
      }
      else { //type1 mppc
        if (argv2 == "fit"){
          sigmas1.push_back(r->Parameter(2)/sqrt(2));
          ye1.push_back(r->ParError(2)/sqrt(2));
          xe1.push_back(0.5);
          channels1.push_back(ih);
          TimeResDist->Fill(sigmas1.back());
          if ((FEB[FEBind] > 17 && FEB[FEBind] < 21) || (FEB[FEBind] > 24)){
            TimeResDistS->Fill(sigmas1.back());
            TimeResDist1->Fill(sigmas1.back());
            TimeResMap_XZ->Fill(
              MapCon[FEB[FEBind]][1][(int)ih%96], // z position
              MapCon[FEB[FEBind]][0][(int)ih%96], // x position
              sigmas1.back());
          }
          else {
            TimeResDistL->Fill(sigmas1.back());
            TimeResMap_ZY->Fill(
              MapCon[FEB[FEBind]][0][(int)ih%96], // z position
              MapCon[FEB[FEBind]][1][(int)ih%96], // y position
              sigmas1.back());
          }
        }
        else{
          sigmas1.push_back(Time_res[ih]->GetStdDev()/sqrt(2));
          ye1.push_back(Time_res[ih]->GetStdDevError()/sqrt(2));
          xe1.push_back(0.5);
          channels1.push_back(ih);
          TimeResDist->Fill(sigmas1.back()/sqrt(2));
          TimeResDist1->Fill(sigmas1.back()/sqrt(2));
          if ((FEB[FEBind] > 17 && FEB[FEBind] < 21) || (FEB[FEBind] > 24)){
            TimeResDistS->Fill(sigmas1.back());
            TimeResMap_XZ->Fill(
              MapCon[FEB[FEBind]][1][(int)ih%96], // z position
              MapCon[FEB[FEBind]][0][(int)ih%96], // x position
              sigmas1.back());
          }
          else{
            TimeResDistL->Fill(sigmas1.back());
            TimeResMap_ZY->Fill(
              MapCon[FEB[FEBind]][0][(int)ih%96], // z position
              MapCon[FEB[FEBind]][1][(int)ih%96], // y position
              sigmas1.back());
          }
        }
      }
      timeRes->cd();
      gStyle->SetOptFit(0111);
      gStyle->SetStatX(0.95);
      gStyle->SetStatH(0.1155555);
      gStyle->SetStatW(0.1);
      c2->Modified(); c2->Update();
      c2->Write();
    }
  }
  c2->Close();

  // plot graph of channel vs Time Res
  FileOutput->cd();
  TCanvas *c3 = new TCanvas("c3","c3",2000,630);
  c3->cd();
  gStyle->SetOptFit(0111);
  gStyle->SetStatH(0.1155555);
  gStyle->SetStatW(0.1);
  //gStyle->SetPalette(0);
  c3->Update(); 
  TGraphErrors *ChannelsTimeRes1 = new TGraphErrors(sigmas1.size(), &channels1[0], &sigmas1[0], &xe1[0], &ye1[0]);
  ChannelsTimeRes1->Draw("AP");
  ChannelsTimeRes1->SetTitle("; Channel; Time Resolution [ns]");
  ChannelsTimeRes1->GetYaxis()->SetRangeUser(0.5,3.0);
  ChannelsTimeRes1->GetXaxis()->SetRangeUser(-0.5,1727.5);
  ChannelsTimeRes1->SetMarkerStyle(kDot);
  TGraphErrors *ChannelsTimeRes2 = new TGraphErrors(sigmas2.size(), &channels2[0], &sigmas2[0], &xe2[0], &ye2[0]);
  ChannelsTimeRes2->SetMarkerStyle(26);
  TGraphErrors *ChannelsTimeRes3 = new TGraphErrors(sigmas3.size(), &channels3[0], &sigmas3[0], &xe3[0], &ye3[0]);
  ChannelsTimeRes3->SetMarkerStyle(24);
  c3->Update();
  TBox *longFibre1 = new TBox(95.5,0.5,287.5,3.0);
  longFibre1->SetFillColorAlpha(kBlue, 0.35);
  TBox *longFibre2 = new TBox(1343.5,0.5,1439.5,3.0);
  longFibre2->SetFillColorAlpha(kBlue, 0.35);
  TBox *longFibre3 = new TBox(959.5,0.5,1055.5,3.0);
  longFibre3->SetFillColorAlpha(kBlue, 0.35);
  TBox *shortFibre1 = new TBox(287.5,0.5,863.5,3.0);
  shortFibre1->SetFillColorAlpha(kYellow, 0.35);
  TBox *shortFibre2 = new TBox(1055.5,0.5,1343.5,3.0);
  shortFibre2->SetFillColorAlpha(kYellow, 0.35);
  TBox *shortFibre3 = new TBox(1439.5,0.5,1727.5,3.0);
  shortFibre3->SetFillColorAlpha(kYellow, 0.35);
  TLegend *legend = new TLegend(0.91,0.7,0.99,0.95);
  legend->AddEntry(longFibre1, "24 cm Fibres", "f");
  legend->AddEntry(shortFibre1, "8 cm Fibres", "f");
  legend->AddEntry(ChannelsTimeRes1, "MPPC Type 1", "p");
  legend->AddEntry(ChannelsTimeRes2, "MPPC Type 2", "p");
  legend->AddEntry(ChannelsTimeRes3, "MPPC Type 3", "p");
  longFibre1->Draw("same");
  longFibre2->Draw("same");
  longFibre3->Draw("same");
  shortFibre1->Draw("same");
  shortFibre2->Draw("same");
  shortFibre3->Draw("same");
  ChannelsTimeRes1->Draw("P");
  ChannelsTimeRes2->Draw("P");
  ChannelsTimeRes3->Draw("P");
  TLine *line = new TLine(0,0,0,0);
  for (Int_t l=1; l<18; l++){
    line = new TLine(96*l-0.5, 0.5, 96*l-0.5, 3.0);
    line->SetLineStyle(9);
    line->Draw("same");
  }
  TGaxis *axis = new TGaxis(-0.5,3.0,1727.5,3.0,0,27,18,"-NMS");
  for (Int_t i=0; i<18; i++)
    axis->ChangeLabel(i+1,-1,-1,-1,-1,-1, to_string(FEB[i]));
  axis->SetTitle("FEB Number");
  axis->SetTitleOffset(1.2);
  axis->SetLabelOffset(0.);
  axis->SetTickSize(-0.03);
  axis->Draw();
  legend->Draw();
  c3->Modified(); c3->Update();
  c3->SaveAs("ChannelTimeRes.C");
  c3->Close();

  TCanvas *c4 = new TCanvas("c4","c4",800,630);
  c4->cd();
  //gStyle->SetStatH(0.1155555);
  //gStyle->SetStatW(0.25);
  c4->Update();

  TimeResDistL->SetFillColorAlpha(kBlue, 0.35);
  TimeResDistL->SetLineColor(kBlack);
  TimeResDistL->SetTitle("");
  TimeResDistL->GetXaxis()->SetRangeUser(0.7, 1.3);
  TimeResDistL->GetYaxis()->SetRangeUser(0, 110);
  TimeResDistL->Draw();
  TimeResDistS->SetFillColorAlpha(kYellow, 0.35);
  TimeResDistS->SetLineColor(kBlack);
  TimeResDistS->Draw("same");
  TimeResDistL->SetStats(0);
  TimeResDistS->SetStats(0);
  TLegend *DistLegend = new TLegend(0.695,0.82,0.954,0.99);
  DistLegend->AddEntry(TimeResDistS, "8 cm Fibres", "f");
  DistLegend->AddEntry(TimeResDistL, "24 cm Fibres", "f");
  DistLegend->Draw();
  TPaveText *pave = new TPaveText(0.694,.503,.955,.819);
  Double_t entriesS = TimeResDistS->GetEntries();
  Double_t meanS = TimeResDistS->GetMean();
  Double_t stdS = TimeResDistS->GetStdDev();
  pave->AddText("8 cm Fibres");
  pave->AddText(Form("Entries %g",entriesS));
  pave->AddText(Form("Mean %.2f ns",meanS));
  pave->AddText(Form("Std Dev %.2f ns",stdS));
  pave->AddLine(0,0.5,1,0.5);
  Double_t entriesL = TimeResDistL->GetEntries();
  Double_t meanL = TimeResDistL->GetMean();
  Double_t stdL = TimeResDistL->GetStdDev();
  pave->AddText("24 cm Fibres");
  pave->AddText(Form("Entries %g",entriesL));
  pave->AddText(Form("Mean %.2f ns",meanL));
  pave->AddText(Form("Std Dev %.2f ns",stdL));
  pave->Paint("NDC");
  pave->Draw();
  c4->Modified(); c4->Update();
  c4->SaveAs("TimeResDistFibres.C");
  c4->Close();

  TCanvas *c5 = new TCanvas("c5","c5",800,630);
  c5->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetMargin(0.059,0.141,0.1,0.06);
  pad1->SetFixedAspectRatio();
  pad1->Draw();
  float label1 = 22/(pad1->GetWh()*pad1->GetAbsHNDC());
  pad1->cd();
  gStyle->SetOptStat(10);
  TimeResMap_XZ->SetTitle("");
  TimeResMap_XZ->SetTitleFont(43, "XYZ"); // title size now in pixels
  TimeResMap_XZ->SetTitleSize(20, "XYZ");
  TimeResMap_XZ->SetTitleOffset(2, "XZ");
  TimeResMap_XZ->SetTitleOffset(1.5, "Y");
  TimeResMap_XZ->GetXaxis()->SetLabelSize(label1);
  TimeResMap_XZ->GetYaxis()->SetLabelSize(label1);
  TimeResMap_XZ->GetZaxis()->SetLabelSize(label1);
  TimeResMap_XZ->Draw();
  c5->Update();
  TPaveStats *st = (TPaveStats*)TimeResMap_XZ->FindObject("stats");
  st->SetY1NDC(0.97);
  st->SetY2NDC(0.999);
  st->SetX1NDC(0.87);
  st->SetX2NDC(0.999);
  TimeResMap_XZ->Draw("colorz");
  c5->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
  pad2->SetMargin(0.059,0.141,0.18,0.02);
  pad2->SetFixedAspectRatio();
  pad2->Draw();
  float label2 = 22/(pad2->GetWh()*pad2->GetAbsHNDC());
  pad2->cd();
  TimeResMap_ZY->SetTitle("");
  TimeResMap_ZY->SetTitleFont(43, "XYZ");
  TimeResMap_ZY->SetTitleSize(20, "XYZ");
  TimeResMap_ZY->SetTitleOffset(1.5, "Y");
  TimeResMap_ZY->SetTitleOffset(4, "X");
  TimeResMap_ZY->GetXaxis()->SetLabelSize(label2);
  TimeResMap_ZY->GetYaxis()->SetLabelSize(label2);
  TimeResMap_ZY->GetZaxis()->SetLabelSize(label2);
  TimeResMap_ZY->Draw();
  c5->Update();
  TPaveStats *st1 = (TPaveStats*)TimeResMap_ZY->FindObject("stats");
  st1->SetY1NDC(0.91);
  st1->SetY2NDC(0.999);
  st1->SetX1NDC(0.87);
  st1->SetX2NDC(0.999);
  TLatex *txt = new TLatex();
  txt->SetNDC();
  txt->SetTextAngle(90);
  txt->SetTextFont(43);
  txt->SetTextSize(20);
  TimeResMap_ZY->Draw("color");
  //txt->DrawLatexNDC(0.92,0.19,"#splitline{Time}{Resolution [ns]}");
  c5->Modified(); c5->Update();
  c5->SaveAs("TimeResMap.C");
  c5->Close();

  TCanvas *c6 = new TCanvas("c6","c6",800,630);
  c6->cd();
  TimeResDist1->SetFillColorAlpha(kYellow, 0.35);
  TimeResDist1->SetLineColor(kBlack);
  TimeResDist1->SetTitle("");
  TimeResDist1->GetXaxis()->SetRangeUser(0.7, 1.3);
  TimeResDist1->GetYaxis()->SetRangeUser(0, 60);
  TimeResDist1->Draw();
  TimeResDist2->SetFillColorAlpha(kBlue, 0.35);
  TimeResDist2->SetLineColor(kBlack);
  TimeResDist2->Draw("same");
  TimeResDist3->SetFillColorAlpha(kRed, 0.35);
  TimeResDist3->SetLineColor(kBlack);
  TimeResDist3->Draw("same");
  TimeResDist1->SetStats(0);
  TimeResDist2->SetStats(0);
  TimeResDist3->SetStats(0);
  TLegend *DistLegend1 = new TLegend(0.695,0.82,0.954,0.99);
  DistLegend1->AddEntry(TimeResDist1, "Type1", "f");
  DistLegend1->AddEntry(TimeResDist2, "Type2", "f");
  DistLegend1->AddEntry(TimeResDist3, "Type3", "f");
  DistLegend1->Draw();
  TPaveText *pave1 = new TPaveText(0.694,.503,.955,0.819);
  Double_t entries1 = TimeResDist1->GetEntries();
  Double_t mean1 = TimeResDist1->GetMean();
  Double_t std1 = TimeResDist1->GetStdDev();
  pave1->AddText("Type1");
  pave1->AddText(Form("Entries %g",entries1));
  pave1->AddText(Form("Mean %.2f ns",mean1));
  pave1->AddText(Form("Std Dev %.2f ns",std1));
  pave1->AddLine(0,0.33,1,0.33);
  Double_t entries2 = TimeResDist2->GetEntries();
  Double_t mean2 = TimeResDist2->GetMean();
  Double_t std2 = TimeResDist2->GetStdDev();
  pave1->AddText("Type2");
  pave1->AddText(Form("Entries %g",entries2));
  pave1->AddText(Form("Mean %.2f ns",mean2));
  pave1->AddText(Form("Std Dev %.2f ns",std2));
  pave1->AddLine(0,0.66,1,0.66);
  Double_t entries3 = TimeResDist3->GetEntries();
  Double_t mean3 = TimeResDist3->GetMean();
  Double_t std3 = TimeResDist3->GetStdDev();
  pave1->AddText("Type3");
  pave1->AddText(Form("Entries %g",entries3));
  pave1->AddText(Form("Mean %.2f ns",mean3));
  pave1->AddText(Form("Std Dev %.2f ns",std3));
  pave1->Paint("NDC");
  pave1->Draw();
  c6->Modified(); c6->Update();
  c6->SaveAs("TimeResDistTypes.C");
  c6->Close();

  TCanvas *c7 = new TCanvas("c7","c7",800,630);
  c7->cd();
  TimeResDist->GetXaxis()->SetRangeUser(0.7, 1.3);
  TimeResDist->GetYaxis()->SetRangeUser(0, 160);
  TimeResDist->Draw();
  c7->Modified(); c7->Update();
  c7->SaveAs("TimeResDist.C");
  c7->Close();

  writeOutput();
  cout << "Execution successful" << endl;
  exit(0);
}
