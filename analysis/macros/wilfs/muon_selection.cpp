#define THIS_NAME MuonSelection
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
TString fileOut        = "";
TString fileIn         = "";
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
Double_t sideTimeRef = 0;
Double_t topTimeRef = 0;
Double_t TimeRef = 0;
vector<Double_t> sideTimeRefVec;
vector<Int_t> sideTimeRefZ;
vector<Double_t> topTimeRefVec;
vector<Int_t> topTimeRefZ;
Int_t sTRcount = 0;
Int_t tTRcount = 0;
Int_t earlyTime = 0;
Int_t lateTime = 0;
Int_t NumberEvDis = 1000;
Int_t eventsPassed = 0;
string foutPNG;
Int_t chargeMin = 25;
Double_t fibreSpeed = 1.7e8; // meters per sec
Double_t attenLen = 4.5; // meters (only short attenuation used as the short fibre lengths mean the long regime is not entered. Value from arxiv.org/pdf/1511.06225.pdf)
string tRef = "median";

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

// function to retrieve the first part of the input filename to add to the output file
string GetLocation(string str)
{
  int i = str.rfind("_events.root");
  string way = str.substr(0,i);  
  return way;
}

// cuts applied to select muon events. Uses event info from XZ and ZY faces. Returns 0 if event passes cuts
bool muon_cut(Event* event, TH2D *eventXZ, TH2D *eventZY, std::vector<ND280SFGDHit*> Hits)
{
  Int_t missedLayers = 0;
  bool minPeak = 0;
  bool largeHitTimeDif = 0;
  for (size_t i=0; i<48; i++){
    Int_t eventDep = event->GetdEdz(i);
    if (eventDep == 0) missedLayers += 1;
    if (eventDep > 50) minPeak = 1;
  }
  size_t nHits = Hits.size();
  for (size_t ihit=0; ihit<nHits; ihit++){
    if (Hits[ihit]->GetToT() > 60) {largeHitTimeDif = 1; break;}
  }
  // cut event if the time over threshold was unusually long
  if (! largeHitTimeDif){ 
    // only save event if enough channels fired
    if (eventXZ->GetEntries()>10 && eventZY->GetEntries()>10){
      // cut event if peak energy deposit less than 50 pe
      if (minPeak){
        // only save event if 2 or less z layers have no hit
        if (missedLayers < 3){
          // only save event if std dev in x and y direction for xy plot is less than 1, to ensure the event is just a particle passing straight through and not recoiling
          if (eventXZ->GetStdDev(1)<1 && eventZY->GetStdDev(2)<1){
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

// cuts applied to select stopping proton events. Uses event info from XZ and ZY faces. Returns 0 if event passes cuts
bool proton_cut(Event *event, TH2D *histXZ, TH2D *histZY, Double_t previousTime)
{
  Int_t missedLayers = 0;
  bool minPeak = 0;
  bool largeHitTimeDif = 0;
  for (size_t i=0; i<48; i++){
    Int_t eventDep = event->GetdEdz(i);
    if (eventDep == 0) missedLayers += 1;
    if (eventDep > 50) minPeak = 1;
  }
  TClonesArray *Hits = 0;
  Hits = event->GetHits();
  Hits->Clear();
  size_t nHits = event->GetNHits();
  for (size_t iev=0; iev<nHits; iev++){
    Hit *hit = (Hit*)Hits->At(iev);
    if (hit->GetToT() > 60) {largeHitTimeDif = 1; break;}
  }
  if (! largeHitTimeDif){
    // ignore event if previous trigger was within 100 ticks (250 ns). This corresponds to two triggers for the same event, so we remove the clone.
    if (abs(event->GetFEB12LeadTime()*2.5 - previousTime) < 250){
      // only save event if enough channels fired
      if (histXZ->GetEntries()>10){
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

void MuonSelection() 
{
  for (int iarg=0; iarg<gApplication->Argc(); iarg++){
    if (string( gApplication->Argv(iarg))=="-h" || string( gApplication->Argv(iarg))=="--help" ){
      cout << "**************************************" << endl;
      cout << "Macros run options:" << endl;
      cout << "   -h || --help      print help info." << endl;
      cout << "   -d || --data      Allows to read beamtest data input." << endl;
      cout << "   -i || --input     input file, including the path(*.root)." << endl;
      cout << "   -o || --output    output file, including the path." << endl;
      cout << "   -b || --batch     Do not shows canvases and closes program at end of execution." << endl;
      cout << "   -f || --flag      Use crosstalk flag to reconstruct voxels." << endl;
      cout << "   -r <method> || --ref <method> Specify the method to extract the time reference (avg, median, mode, front, trigger).\n" << endl;      
      cout << "Display options: " << endl;
      cout << "   -strue            Shows true event displays." << endl;
      cout << "   -sreco            Shows reco event displays." << endl;
      cout << "**************************************" << endl;
    }
    else if (string( gApplication->Argv(iarg))=="-d" || string( gApplication->Argv(iarg))=="--data" ){
      IsMC = false;
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
    }
    else if (string( gApplication->Argv(iarg))=="-r" || string( gApplication->Argv(iarg))=="--ref" ){
      tRef = gApplication->Argv(iarg+1);
    }
    else if (string( gApplication->Argv(iarg))=="-strue"){
      SHOW_TRUE = true;
    }
    else if (string( gApplication->Argv(iarg))=="-sreco"){
      SHOW_RECO = true;
    }
  }
  if(IsMC && IsPrototype){
    SFGD_X             = 24;
    SFGD_Y             = 8;
    SFGD_Z             = 48;
  }

  // define FEB numbers 
  Int_t FEBs[28] = {0,1,2,3,4,0,0,0,5,6,7,8,0,0,0,0,9,10,11,12,13,0,0,0,14,15,16,17};

  // get input and output files
  if (fileOut=""){
    fileOut=GetLocation(fileIn.Data());
    fileOut+="_muons.root";
  }
  cout << "Output: " << fileOut << endl << endl;
  
  //set up TFiles FileInput and FileOutput and link to TTrees dataOut and data
  linkFilesToTTrees();

  // TDirectory *SubFolder = FileInput->GetDirectory("TimeRes");
  // if (!SubFolder) SubFolder = FileInput->mkdir("TimeRes","TimeRes");

  TTree *SelectedEvents = new TTree("SelectedEvents", "SelectedEvents");
  Event* SelEvent = new Event();
  SelectedEvents->Branch("Event", "Event", SelEvent);

  cout << "Number of events: " << nEvents << endl;

  gStyle->SetPalette(0);

  Int_t eventNum = 0;
  Double_t previousTime = 0;

  //event_XY->SetTitle("; x position [cm]; y position [cm]; Light Yield [p.e.]");
  //event_XY->GetXaxis()->SetTitleSize(0.05);
  //event_XY->GetYaxis()->SetTitleSize(0.05);
  //event_XY->GetZaxis()->SetTitleSize(0.05);
  //event_XY->GetXaxis()->SetLabelSize(0.05);
  //event_XY->GetYaxis()->SetLabelSize(0.05);
  //event_XY->GetZaxis()->SetLabelSize(0.05);
  //event_XY->GetZaxis()->SetRangeUser(0.,140.);
  
  TH2D *event_ZY = new TH2D("","", 48,0,48, 8,0,8);
  event_ZY->SetTitle("; x position [cm]; y position [cm]; Light Yield [p.e.]");
  event_ZY->GetXaxis()->SetTitleSize(0.05);
  event_ZY->GetYaxis()->SetTitleSize(0.05);
  event_ZY->GetZaxis()->SetTitleSize(0.05);
  event_ZY->GetXaxis()->SetLabelSize(0.05);
  event_ZY->GetYaxis()->SetLabelSize(0.05);
  event_ZY->GetZaxis()->SetLabelSize(0.05);
  //event_ZY->GetZaxis()->SetRangeUser(0.,140.);
  
  TH2D *event_XZ = new TH2D("","", 24,0,24, 48,0,48);
  event_ZY->SetTitle("; x position [cm]; z position [cm]; Light Yield [p.e.]");
  event_XZ->GetXaxis()->SetTitleSize(0.05);
  event_XZ->GetYaxis()->SetTitleSize(0.05);
  event_XZ->GetZaxis()->SetTitleSize(0.05);
  event_XZ->GetXaxis()->SetLabelSize(0.05);
  event_XZ->GetYaxis()->SetLabelSize(0.05);
  event_XZ->GetZaxis()->SetLabelSize(0.05);
  //event_XZ->GetZaxis()->SetRangeUser(0.,140.);

  // retrieve time walk offset parameters for each channel
  TFile *f = new TFile("TimeOffset.root", "READ");
  TTree *Offtree = (TTree*)f->Get("Events");
  Float_t TimeOffsets[5*1728];
  Offtree->SetBranchAddress("TimeOffsets",&TimeOffsets[0]);
  Offtree->GetEvent();
  f->Close();
  //for (int i=0; i<5; i++)
    //cout << TimeOffsets[i] << endl;

  size_t maxEv = nEvents;

  for (size_t iev = 0; iev<maxEv; iev++) {
    if(selEvents >= maxSelEvents) break;
    cout << "Selecting events..." << (int)(iev*100/(Double_t)maxEv) << "% done        \r";
    std::cout.flush();

    std::vector<ND280SFGDHit*> mppc_hits;

    //returns a vector of ND280SFGDHit objects that contain hit properties.
    //Skips first 20 events and ignores all hits prior to one with a negative
    //charge
    mppc_hits = getEventMPPChits(iev);

    sideTimeRef = 0;
    topTimeRef = 0;
    TimeRef = 0;
    earlyTime = 0;
    lateTime = 0;
    sideTimeRefVec.clear();
    sideTimeRefZ.clear();
    topTimeRefVec.clear();
    topTimeRefZ.clear();
    sTRcount = 0;
    tTRcount = 0;
    Int_t element = 0;

    Double_t triggerTime = 0;
    size_t nHits = mppc_hits.size();
    if (nHits < 90) continue;
    // iterate through hits of a single event
    for (size_t ihit=0; ihit<nHits; ihit++){
      Double_t timeSinceTrigger = mppc_hits[ihit]->GetDt()*2.5; //convert to ns
      Double_t hitTime = mppc_hits[ihit]->GetTfromSpill()*2.5;
      triggerTime = hitTime - timeSinceTrigger;
      if (abs(triggerTime - previousTime) < 250) break; //skip event if it's trigger time is similar to the previous event's, meaning there was probably two triggers for a single event
      Int_t hitCharge = mppc_hits[ihit]->GetPE();
      Int_t hitFEB = mppc_hits[ihit]->GetFEB();
      Int_t hitView = mppc_hits[ihit]->GetView();
      Int_t hitZ = mppc_hits[ihit]->GetZ();
      Int_t hitX = mppc_hits[ihit]->GetX();
      Int_t hitY = mppc_hits[ihit]->GetY();
      Int_t chNum = mppc_hits[ihit]->GetCh();
      Double_t travelDist = 0;
      Double_t calibTime = 0;
      Int_t pairCharge = 0;
      bool keepHit = false;

      // cut on time window for muons (and electrons and pions)
      if (timeSinceTrigger < 300 && hitCharge > chargeMin){
        //&& timeSinceTrigger > 250){
        //Hit* SelHit = SelEvent->AddHit();
        //SelHit->SetAll(hit);
        // find distance light travelled to reach MPPC
        if (hitView == 2){ // side channels
          for (size_t ihit2=0; ihit2<nHits; ihit2++){ //find hit in other view with same z value
            if (mppc_hits[ihit2]->GetView() == 2 || mppc_hits[ihit2]->GetView() == 0) continue;
            if (hitZ == mppc_hits[ihit2]->GetZ()){
              if (mppc_hits[ihit2]->GetPE() > pairCharge){
                pairCharge = mppc_hits[ihit2]->GetPE();
                keepHit = true; // keep hit only if there is a hit in the same z layer of the top/bottom channels
                if (hitZ % 2 == 0) // right channels
                  travelDist = (24.-mppc_hits[ihit2]->GetX()+0.5)*0.01; // meters
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
                if (hitZ % 2 == 0) // bottom channels
                  travelDist = (mppc_hits[ihit2]->GetY()+0.5)*0.01; // meters
                else // top channels
                  travelDist = (8.-mppc_hits[ihit2]->GetY()+0.5)*0.01; // meters
              }
            }
          }
        }
        else keepHit = true;
        // account for light attenuation with hit charge
        Int_t hitChargeCalib = hitCharge*exp(-travelDist/attenLen);
        // calibrate hit time to account for time-walk effect, fibre travel time and light attenuation
        element = FEBs[hitFEB]*96*5+chNum*5;
        if (hitFEB == 8 || hitFEB == 11) //type 3 MPPCs have a different fit
          calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(-hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
        else
          calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
        //calibTime = hitTime;
        if (keepHit){
          if (hitCharge < 10000 && hitCharge < chargeMin){ // leave out hits with very large or small amplitudes
            // fill event time histo
            //event_Time->Fill(triggerTime - calibTime);
            // front and back: plot event using mapping info 
            //if (hitView == 0){
            //  event_XY->Fill(hitX, hitY, hitChargeCalib);
            //} 
            // sides: plot event
            if (hitView == 2){
              event_ZY->Fill(hitZ, hitY, hitChargeCalib);
              //event_Time_ZY->Fill(triggerTime - calibTime);
              if (tRef == "avg"){
                sideTimeRef += calibTime;
                sTRcount += 1;
              }
              if (tRef == "median" || tRef == "mode"){
                sideTimeRefVec.push_back(calibTime); // round time to 2 dp so we can find a mode
                //cout << sideTimeRefVec.back() << endl;
                sideTimeRefZ.push_back(hitZ);
              }
              // fill MIP plot with deposited charge (dE/dx) data. Only use data from sides/top so we can see the change in energy along the detector
              //event_MIP->Fill(FEB[FEBs[i]].hitCharge_pe->at(check));
            } 
            // top and bottom: plot event
            else if (hitView == 1){
              event_XZ->Fill(hitX, hitZ, hitCharge);
              //event_Time_XZ->Fill(triggerTime - calibTime);

              if (tRef == "avg"){
                topTimeRef += calibTime;
                tTRcount += 1;
              }
              if (tRef == "median" || tRef == "mode"){
                topTimeRefVec.push_back(calibTime); // round time to 2 dp so we can find a mode
                topTimeRefZ.push_back(hitZ);
              }
              // fill MIP plot with deposited charge (dE/dx) data
              //event_MIP->Fill(FEB[FEBs[i]].hitCharge_pe->at(check));
            }
          } //charge window
        } //keepHit
      } //time window cut
    } //loop through hits

    // only take events that look go straight through the detector with no interactions
    if (muon_cut(unpackEvent, event_XZ, event_ZY, mppc_hits) == 0){
      if (tRef == "avg"){
        sideTimeRef = sideTimeRef/sTRcount; // time tReference is the average hit time of the event
        topTimeRef = topTimeRef/tTRcount;
      }
      if (tRef == "median"){
        sort(sideTimeRefZ.begin(), sideTimeRefZ.end());
        sort(topTimeRefZ.begin(), topTimeRefZ.end());
        combo median = findMedian(sideTimeRefVec);
        sideTimeRef = median.result;
        //cout << "median:" << endl;
        //cout << median.result << endl;
        //cout << "Element:" << endl;
        //cout << median.elem[0] << endl;
        //refMap_ZY->Fill(sideTimeRefZ[median.elem[0]]);
        median = findMedian(topTimeRefVec);
        topTimeRef = median.result;
        //refMap_XZ->Fill(topTimeRefZ[median.elem[0]]);
        //topTimeRefVec.insert(end(topTimeRefVec),begin(sideTimeRefVec),end(sideTimeRefVec));
        //TimeRef = findMedian(topTimeRefVec);
      }
      if (tRef == "mode"){
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
          //refMap_ZY->Fill(sideTimeRefZ[midmode.result]);
        }
        mode = findMode(topTimeRefVec);
        topTimeRef = mode.result;
        modeLayer.clear();
        for (Int_t k=0; k<(int)mode.elem.size(); k++)
          modeLayer.push_back(topTimeRefZ[mode.elem[k]]);
        if (mode.elem[0] != -2){
          combo midmode = findMedian(modeLayer);
          //refMap_XZ->Fill(topTimeRefZ[midmode.result]);
        }
      }

      //timeRef = FEB[12].hitTimefromSpill->at(TOFtrigger);
   
      SelEvent->SetEventID(iev); 
      SelEvent->SetFEB12ch(unpackEvent->GetFEB12ch());
      SelEvent->SetFEB12LeadTime(topTimeRef); //store event time references haphazardly
      SelEvent->SetRange(sideTimeRef);
      SelectedEvents->Fill();
      selEvents++;
    }
    //event_XY->Reset();
    event_XZ->Reset();
    event_ZY->Reset();
    SelEvent->Clear();
    for (int it=0; it<mppc_hits.size(); it++) delete mppc_hits[it];
  }
  cout << "\nSelecting events...100% done" << endl;
  FileOutput->cd();
  SelectedEvents->Write("",TObject::kOverwrite);
  //dataOut->Fill();
  //dataOut->Print();
  //dataOut->Write("", TObject::kOverwrite);
  f->Close();
  FileOutput->Close();
  FileInput->Close();
  exit(0);
}
