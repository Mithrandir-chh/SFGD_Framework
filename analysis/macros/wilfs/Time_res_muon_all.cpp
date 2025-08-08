#define THIS_NAME TimeResMuonAll
#define NOINTERACTIVE_INPUT
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
int IsPrototype        = false;
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
Int_t chargeMin = 25;
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
#include "TImage.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TBox.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLeaf.h"

using namespace std;

// function to retrieve the first part of the input filename to add to the output file
string GetLocation(string str)
{
  int i = str.rfind("_Slot_");
  string way = str.substr(0,i);  
  return way;
}

// function to locate the directory of a given file
string GetDir(string str)
{
    int i = str.rfind("/");
    string way = str.substr(0,i);
    return way;
}

// cuts applied to select muon events. Uses event info from XZ and ZY faces. Returns 0 if event passes cuts
bool muon_cut(TH2F *histXZ, TH2F *histZY, TH2F *histXY)
{
  if (LargehitTimeDif == 0){
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
              return 1;
            }
          }
        }
      }
    }
  }
  return 0;
}

// cuts applied to select stopping proton events. Uses event info from XZ and ZY faces. Returns 0 if event passes cuts
bool proton_cut(TH2F *histXZ, TH2F *histZY, TH2F *histXY)
{
  if (LargehitTimeDif == 0){
    // ignore event if previous trigger was within 100 ticks (250 ns). This corresponds to two triggers for the same event, so we remove the clone.
    if (clonedEvent == 0){
      // only save event if enough channels fired
      if (histXZ->GetEntries()>10){
        // cut event if peak energy deposit less than 50 pe
        if (minPeak == 1){
          // only save event if 2 or less z layers have no hit
          if (missedLayers < 3){
            // only save event if std dev in x and y direction for xy plot is less than 1
            if (histXZ->GetStdDev(1)<1 && histZY->GetStdDev(2)<1 && histXY->GetStdDev(1)<1 && histXY->GetStdDev(2)<1){
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
  Int_t maxIter;
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

void TimeResMuonAll()
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
    cerr << "'notest': create standard analysis file, overwriting previous ones.\n" << endl;
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
    cerr << "'notest': create standard analysis file, overwriting previous ones.\n" << endl;
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
    fileOut += "_TimeRes_"+argv1+"_"+argv2+"_allChannels_noFEB18.root";
  cout << "Output: " << fileOutput << endl;

  // open a new TFile with our generated name. We get the data from here
  TFile *FileInput = new TFile(rootFileInput.c_str());
  
  // initiate our custom struct to store TTree data
  vectorsTree FEB[NumberOfEB];
  
  // give each element of the struct a zero value
  for (Int_t i=0;i<NumberOfEB;i++){
    FEB[i].FEBSN=0;
    FEB[i].SpillTag=0;
    FEB[i].hitsChannel=0;
    FEB[i].hitAmpl=0;
    FEB[i].hitAmplRec=0;
    FEB[i].hitLeadTime=0;
    FEB[i].GTrigTag=0;
    FEB[i].GTrigTime=0;
    FEB[i].hitLGAmpl=0;
    FEB[i].hitLGAmplRec=0;
    FEB[i].hitTrailTime=0;
    FEB[i].hitTimeDif=0;
    FEB[i].hitTimefromSpill=0;
    
    FEB[i].hitHG_pe=0;
    FEB[i].hitLG_pe=0;
    FEB[i].hitToT_pe=0;
    FEB[i].hitCharge_pe=0;      
  }
  
  // initialise an array of TTrees
  TTree *FEBtree[NumberOfEB];
  Long64_t nentries[NumberOfEB];
  std::fill(nentries, nentries + NumberOfEB, 0);
  
  ostringstream sFEBnum;
  string sFEB;
  ostringstream sChannelnum;
  string sChannel;
  ostringstream sEventnum;
  string sEvent;
  ostringstream spillNum;
  
  // set the branches of our TTrees using the branch addresses of the TTrees in the input file
  for (Int_t ih=0; ih<NumberOfEB; ih++) {
    sFEBnum.str("");
    sFEBnum << ih;
    sFEB = "FEB_"+sFEBnum.str();
    FEBtree[ih] = (TTree*)FileInput->Get(sFEB.c_str());
    if ((TTree*)FileInput->Get(sFEB.c_str())){
      FEBtree[ih]->SetBranchAddress((sFEB+"_SN").c_str(),&FEB[ih].FEBSN);
      FEBtree[ih]->SetBranchAddress((sFEB+"_SpillTag").c_str(),&FEB[ih].SpillTag);
      FEBtree[ih]->SetBranchAddress((sFEB+"_GTrigTag").c_str(),&FEB[ih].GTrigTag);
      FEBtree[ih]->SetBranchAddress((sFEB+"_GTrigTime").c_str(),&FEB[ih].GTrigTime);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitsChannel").c_str(),&FEB[ih].hitsChannel);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitAmpl").c_str(),&FEB[ih].hitAmpl);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitLGAmpl").c_str(),&FEB[ih].hitLGAmpl);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitLeadTime").c_str(),&FEB[ih].hitLeadTime);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitTrailTime").c_str(),&FEB[ih].hitTrailTime);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitTimeDif").c_str(),&FEB[ih].hitTimeDif);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitTimefromSpill").c_str(),&FEB[ih].hitTimefromSpill);
      
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitAmplRecon").c_str(), &FEB[ih].hitAmplRec);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitLGAmplRecon").c_str(), &FEB[ih].hitLGAmplRec);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitHG_pe").c_str(), &FEB[ih].hitHG_pe);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitLG_pe").c_str(), &FEB[ih].hitLG_pe);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitToT_pe").c_str(), &FEB[ih].hitToT_pe);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitCharge_pe").c_str(), &FEB[ih].hitCharge_pe);
      
      nentries[ih] = FEBtree[ih]->GetEntries();
      FEBtree[ih]->GetEntry(0);
      //std::cout << "Number of events = " <<FEB[ih].FEBSN->size()<<std::endl;
    }
  }

  // find the minimum number of entries, i.e. the number of spills
  double minEn = nentries[0];
  for (int i = 0; i < NumberOfEB; i++){
    if (nentries[i] < minEn && nentries[i] > 0)
      minEn = nentries[i];
  }
  cout << "Number of spills " << minEn << endl;

  // open the output file to write to
  TFile wfile(rootFileOutput.c_str(), "recreate");
  cout<<rootFileOutput<<endl;
  //TFile wfile2("AmpTime.root", "recreate");
  
  // create a 3D matrix to store mapping data for our event displays
  int MapCon[28][2][96];

  // initiate plots for channel time res measurements. One plot per channel (18 FEBs with 96 channels each)
  TH1F *Time_res;
  //TH2F *AmpTime[1728];
  // retrieve mapping info and put into 3D matrix
  for (int iFEB = 0; iFEB<19; iFEB++) {
    if (FEBs[iFEB] != 12){ // FEB_12 was used for TOFtrigger, not reading data from the prototype
      sFEBnum.str("");
      sFEBnum << FEBs[iFEB];
      // retrieve a txt file with mapping information for the FEBs
      // txt file has 3 columns: channel no., x position, y position
      sFEB = "../mapping/" + sFEBnum.str() + ".txt";
      ifstream fmap(sFEB.c_str());
      //cout <<endl<< "FEB "<< FEBs[iFEB]<< " mapping"<<endl;
      // read mapping info into 3D matrix
      int temp=0;
      while (!fmap.eof()) {
        fmap >> temp >> MapCon[FEBs[iFEB]][0][temp] >> MapCon[FEBs[iFEB]][1][temp]; // COLUMNS: temp=(ch. no.), FEBs=(FEB no.), 1/0=(vertical/horizontal position) 
        // so MapCon[2][0][26] will be the horizontal position of the cube in FEB_2 at channel 26
      }
      fmap.close();
    }
  }
  // name and define time res plot
  if (argv1 == "trigger")
    Time_res = new TH1F("AllChannels_TimeResolution","AllChannels_TimeResolution",40,230,330);
  else
    Time_res = new TH1F("AllChannels_TimeResolution","AllChannels_TimeDifference",40,-50,50);
  Time_res->GetXaxis()->SetTitle("Time Difference [ns]");
  Time_res->GetYaxis()->SetTitle("Counts");
  //sChannel = "FEB"+sFEBnum.str()+"_channel"+sChannelnum.str()+"_AmpTime";
  //AmpTime[96*iFEB+ih] = new TH2F(sChannel.c_str(),sChannel.c_str(), 300,0,300, 200,0,500);
  //AmpTime[96*iFEB+ih]->GetXaxis()->SetTitle("Hit Charge [p.e.]");
  //AmpTime[96*iFEB+ih]->GetYaxis()->SetTitle("Hit Time before Trigger [ns]");
  
  // initiate our event displays for all events with relevant titles and dimensions
  TH2F *EventsMap_XY = new TH2F("All_events_map_XY","All_events_map_XY",  24,0,24, 8,0,8);
  EventsMap_XY->GetXaxis()->SetTitle("x position [cm]");
  EventsMap_XY->GetYaxis()->SetTitle("y position [cm]");
  TH2F *EventsMap_ZY = new TH2F("All_events_map_ZY","All_events_map_ZY",  48,0,48, 8,0,8);
  EventsMap_ZY->GetXaxis()->SetTitle("z position [cm]");
  EventsMap_ZY->GetYaxis()->SetTitle("y position [cm]");
  TH2F *EventsMap_XZ = new TH2F("All_events_map_XZ","All_events_map_XZ",  24,0,24, 48,0,48);
  EventsMap_XZ->GetXaxis()->SetTitle("x position [cm]");
  EventsMap_XZ->GetYaxis()->SetTitle("z position [cm]");



  // create directories in our output file to store our event displays and MIP plots
  TDirectory *events2D = wfile.mkdir("events2D");
 
  // initiate canvases to plot individual event displays, MIP peaks and time resolution plots
  TCanvas *c1 = new TCanvas("c1","c1",800,3000);
  TCanvas *c3 = new TCanvas("c3","c3",800,600);

  bool SpillMissed = false;

  Double_t energyDep[48];
  Double_t energyDepXZ[48];
  Double_t energyDepZY[48];
  Int_t eventNum = 0;
  int previousTime = 0;

  // retrieve walk time offset parameters for each channel
  TFile *f = new TFile("TimeOffset.root", "READ");
  TTree *tree = (TTree*)f->Get("tree");
  Float_t TimeOffsets[5*1728];
  tree->SetBranchAddress("TimeOffsets",&TimeOffsets[0]);
  tree->GetEvent();
  for (int i=0; i<5; i++)
    cout << TimeOffsets[i] << endl;

  if (argv3 == "test")
    minEn = 1;
  // iterate through spills
  for (Int_t subSpill = 0; subSpill<minEn; subSpill++) { 
  //for (Int_t subSpill = 0; subSpill<10; subSpill++) {
    // determine if the spill has any data. Some spills don't have (much) data due to problems with the data transfer rate and were skipped in the data unpacking 
    cout << "Getting Spill Number " << subSpill + 1 << endl;
    for (int ik = 0; ik < 19; ik++){
      FEBtree[FEBs[ik]]->GetEntry(subSpill);
      if (FEB[FEBs[ik]].SpillTag->back() != subSpill + 1)
        cout << "wtf" << endl;
      if (FEB[FEBs[ik]].SpillTag->size() < 2){
        cout << "NULL" << endl;
        SpillMissed = true;
        break;
      } 
      else{ 
        SpillMissed = false;
        //cout << "FEB_" << FEBs[ik] << " "<< FEB[FEBs[ik]].SpillTag->size() << endl; //
        //cout << FEB[12].FEBSN->size();
      }
    }

    // process the data, if present, and fill event displays
    if (!SpillMissed){
      cout << FEB[12].FEBSN->size() << " triggers in spill" << endl;
      int numberEvents = 0;
      // iterate through TOF triggers (i.e. events)
      for ( int TOFtrigger = 0; TOFtrigger < FEB[12].FEBSN->size(); TOFtrigger++){ // remember FEB_12 contains the TOF/trigger info
        //if (FEB[12].hitTimeDif->at(TOFtrigger) > 0 && NumberEvDis > eventNum){ // require time diff info from TOF to make sense and don't plot more than "NumberEvDis" events
        if (FEB[12].hitTimeDif->at(TOFtrigger) > 0){
          for (int ij = 0; ij < 48; ij++ ){
            energyDep[ij] = 0;
            energyDepXZ[ij] = 0;
            energyDepZY[ij] = 0;
          }
          // initiate event displays for individual events
          TH2F *event_XY;
          sEventnum.str("");
          spillNum.str("");
          sEventnum << numberEvents+1;
          spillNum << subSpill+1;
          sEvent = "event_S"+spillNum.str()+"_E"+sEventnum.str()+"_XY(front)";
          event_XY = new TH2F(sEvent.c_str(),sEvent.c_str(), 24,0,24, 8,0,8);
          event_XY->GetXaxis()->SetTitle("x position [cm]");
          event_XY->GetYaxis()->SetTitle("y position [cm]");
          event_XY->GetZaxis()->SetTitle("Light Yield [p.e.]");
          event_XY->GetZaxis()->SetRangeUser(0.,140.);
          
          TH2F *event_ZY;
          sEvent = "event_S"+spillNum.str()+"_E"+sEventnum.str()+"_ZY(side)";
          event_ZY = new TH2F(sEvent.c_str(),sEvent.c_str(), 48,0,48, 8,0,8);
          event_ZY->GetXaxis()->SetTitle("z position [cm]");
          event_ZY->GetYaxis()->SetTitle("y position [cm]");
          event_ZY->GetZaxis()->SetTitle("Light Yield [p.e.]");
          event_ZY->GetZaxis()->SetRangeUser(0.,140.);
          
          TH2F *event_XZ;
          sEvent = "event_S"+spillNum.str()+"_E"+sEventnum.str()+"_XZ(top)";
          event_XZ = new TH2F(sEvent.c_str(),sEvent.c_str(), 24,0,24, 48,0,48);
          event_XZ->GetXaxis()->SetTitle("x position [cm]");
          event_XZ->GetYaxis()->SetTitle("z position [cm]");
          event_XZ->GetZaxis()->SetTitle("Light Yield [p.e.]");
          event_XZ->GetZaxis()->SetRangeUser(0.,140.);
         
          // initiate TH1F to store light yield data
          TH1F *event_LY;
          sEvent = "event_LY_S"+spillNum.str()+"_E"+sEventnum.str();
          event_LY = new TH1F(sEvent.c_str(),sEvent.c_str(),48,0,48);
          event_LY->GetXaxis()->SetTitle("z position [cm]");
          event_LY->GetYaxis()->SetTitle("Light Yield [p.e.]");

          //TH1F *event_MIP;
          //sEvent = "event_S"+spillNum.str()+"_E"+sEventnum.str()+"_MIP";
          //event_MIP = new TH1F(sEvent.c_str(),sEvent.c_str(),200,0,200);
          //event_MIP->GetXaxis()->SetTitle("hitCharge [p.e.]");
          //event_MIP->GetYaxis()->SetTitle("Counts");

          TH1F *event_LY_XZ;
          sEvent = "event_LY_S"+spillNum.str()+"_E"+sEventnum.str()+"_XZ";
          event_LY_XZ = new TH1F(sEvent.c_str(),sEvent.c_str(),48,0,48);
          event_LY_XZ->SetXTitle ("z position [cm]");
          event_LY_XZ->SetYTitle ("Light Yield [p.e.]");


          TH1F *event_LY_ZY;
          sEvent = "event_LY_S"+spillNum.str()+"_E"+sEventnum.str()+"_ZY";
          event_LY_ZY = new TH1F(sEvent.c_str(),sEvent.c_str(),48,0,48);
          event_LY_ZY->SetXTitle ("z position [cm]");
          event_LY_ZY->SetYTitle ("Light Yield [p.e.]");


          TH1F *event_Time;
          sEvent = "event_Time"+sEventnum.str();
          sEvent = "event_Time_S"+spillNum.str()+"_E"+sEventnum.str();
          event_Time = new TH1F(sEvent.c_str(),sEvent.c_str(), 180,0,450);
          event_Time->SetXTitle("Hit time before trigger [ns]");
          event_Time->SetYTitle("Count");


          TH1F *event_Time_XZ;
          sEvent = "event_Time_S"+spillNum.str()+"_E"+sEventnum.str()+"_XZ";
          event_Time_XZ = new TH1F(sEvent.c_str(),sEvent.c_str(), 180,0,450);
          event_Time_XZ->SetXTitle("Hit time before trigger [ns]");
          event_Time_XZ->SetYTitle("Count");


          TH1F *event_Time_ZY;
          sEvent = "event_Time_S"+spillNum.str()+"_E"+sEventnum.str()+"_ZY";
          event_Time_ZY = new TH1F(sEvent.c_str(),sEvent.c_str(), 180,0,450);
          event_Time_ZY->SetXTitle("Hit time before trigger [ns]");
          event_Time_ZY->SetYTitle("Count");

          LargehitTimeDif = 0;
          clonedEvent = 0;
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
          Int_t GTindex[2] = {0,0};
          Int_t GTindex2[2] = {0,0};
          Int_t element;
          // iterate through FEBs, missing out FEB_12
          for (int i = 0; i < 19; i++){
            if (FEBs[i]!=12){
              // get lower and upper positions of GTrigTags that are equal to the GTrigTag of the TOF trigger, i.e. get signals from a single event
              auto bounds=std::equal_range(FEB[FEBs[i]].GTrigTag->begin(), FEB[FEBs[i]].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
              GTindex[0] = bounds.first - FEB[FEBs[i]].GTrigTag->begin(); // element no. of first GTrigTag that equals trigger's
              GTindex[1] = bounds.second - FEB[FEBs[i]].GTrigTag->begin(); // element no. of last GTrigTag that equals trigger's
              // iterate through elements of a single event
              for (int check = GTindex[0]; check < GTindex[1]; check++){
                Double_t triggerTime = FEB[12].hitTimefromSpill->at(TOFtrigger)*2.5; //convert to ns
                Int_t hitCharge = FEB[FEBs[i]].hitCharge_pe->at(check);
                Int_t HGCharge = FEB[FEBs[i]].hitHG_pe->at(check);
                Double_t hitTime = FEB[FEBs[i]].hitTimefromSpill->at(check)*2.5;
                Int_t hitChannel = FEB[FEBs[i]].hitsChannel->at(check);
                Double_t travelDist = 0;
                Double_t calibTime = 0;
                Int_t pairCharge = 0;
                bool keepHit = false;

                // find distance light travelled to reach MPPC
                if (FEBs[i] == 1 || FEBs[i] == 2 || FEBs[i] == 17 || FEBs[i] == 24){ // side channels
                  for (int j=0; j<28; j++){
                    if (FEBpairs[j] == FEBs[i]){
                      Int_t pair = j;
                      auto bounds2=std::equal_range(FEB[pair].GTrigTag->begin(), FEB[pair].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
                      GTindex2[0] = bounds2.first - FEB[pair].GTrigTag->begin(); 
                      GTindex2[1] = bounds2.second - FEB[pair].GTrigTag->begin();
                      for (int check2 = GTindex2[0]; check2 < GTindex2[1]; check2++){
                        if (MapCon[FEBs[i]][0][hitChannel] == MapCon[pair][1][(int)FEB[pair].hitsChannel->at(check2)]){
                          if (FEB[pair].hitCharge_pe->at(check2) > pairCharge){
                            pairCharge = FEB[pair].hitCharge_pe->at(check2);
                            keepHit = true; // keep hit only if there is a hit in the same z layer of the top/bottom channels
                            if (FEBs[i] == 1 || FEBs[i] == 24) // right channels
                              travelDist = (24.-MapCon[pair][0][(int)FEB[pair].hitsChannel->at(check2)]-0.5)*0.01; // meters
                            else // left channels
                              travelDist = (MapCon[pair][0][(int)FEB[pair].hitsChannel->at(check2)]+0.5)*0.01; // meters
                          }
                        }
                      }
                    }
                  }
                }
                else if (FEBs[i] != 0 && FEBs[i] != 16) { // top/bottom channels
                  Int_t pair = FEBpairs[FEBs[i]];
                  auto bounds2=std::equal_range(FEB[pair].GTrigTag->begin(), FEB[pair].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
                  GTindex2[0] = bounds2.first - FEB[pair].GTrigTag->begin(); 
                  GTindex2[1] = bounds2.second - FEB[pair].GTrigTag->begin();
                  for (int check2 = GTindex2[0]; check2 < GTindex2[1]; check2++){
                    if (MapCon[FEBs[i]][1][hitChannel] == MapCon[pair][0][(int)FEB[pair].hitsChannel->at(check2)]){
                      if (FEB[pair].hitCharge_pe->at(check2) > pairCharge){
                        pairCharge = FEB[pair].hitCharge_pe->at(check2);
                        keepHit = true;
                        if (FEBs[i] == 3 || FEBs[i] == 4 || FEBs[i] == 8 || FEBs[i] == 18 || FEBs[i] == 19 || FEBs[i] == 20) // top channels
                          travelDist = (8.-MapCon[pair][1][(int)FEB[pair].hitsChannel->at(check2)]-0.5)*0.01; // meters
                        else // bottom channels
                          travelDist = (MapCon[pair][1][(int)FEB[pair].hitsChannel->at(check2)]+0.5)*0.01; // meters
                      }
                    }
                  }
                }
                else keepHit = true;
                // account for light attenuation with hit charge
                Int_t hitChargeCalib = hitCharge*exp(-travelDist/attenLen);
                // calibrate hit time to account for time-walk effect, fibre travel time and light attenuation
                if (FEBs[i] < 12)
                  element = i*96*5+hitChannel*5;
                else
                  element = (i-1)*96*5+hitChannel*5;
                if (FEBs[i] == 8 || FEBs[i] == 11)
                  calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(-hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
                else
                  calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
                //calibTime = hitTime;
                // cut on time window for muons (and electrons and pions)
                if (triggerTime - calibTime < 300){
                    //&& triggerTime - calibTime > 250){
                  //if (hitCharge > chargeMin){
                    if (keepHit){
                      // fill event time histo
                      event_Time->Fill(triggerTime - calibTime);
                      // ensure time since previous trigger, or the trigger before that, isn't short
                      if (abs(triggerTime - previousTime) < 250)
                        clonedEvent = 1;
                      // record if time length of signal is long
                      if (FEB[FEBs[i]].hitTimeDif->at(check) > 60)
                        LargehitTimeDif = 1;
                      // front and back: plot event using mapping info 
                      if (FEBs[i] == 0 || FEBs[i] == 16){
                        if (hitCharge < 10000) // leave out hits with very large or negative amplitudes
                          event_XY->Fill( 
                              MapCon[FEBs[i]][0][hitChannel], // x position
                              MapCon[FEBs[i]][1][hitChannel], // y position
                              hitChargeCalib // signal amplitude in pe
                              );
                      } 
                      // sides: plot event
                      else if (FEBs[i] == 1 || FEBs[i] == 2 || FEBs[i] == 17 || FEBs[i] == 24){
                        if (hitCharge < 10000){
                          event_ZY->Fill(
                              MapCon[FEBs[i]][0][hitChannel], // z position
                              MapCon[FEBs[i]][1][hitChannel], // y position
                              hitChargeCalib
                              );
                          event_Time_ZY->Fill(triggerTime - calibTime);
                          if (argv1 == "avg"){
                            if (hitCharge > chargeMin){
                              sideTimeRef += calibTime;
                              sTRcount += 1;
                            }
                          }
                          if (argv1 == "median" || argv1 == "mode"){
                            if (hitCharge > chargeMin){
                              sideTimeRefVec.push_back(calibTime); // round time to 2 dp so we can find a mode
                              //cout << sideTimeRefVec.back() << endl;
                              sideTimeRefZ.push_back(MapCon[FEBs[i]][0][hitChannel]);
                            }
                          }
                          // fill MIP plot with deposited charge (dE/dx) data. Only use data from sides/top so we can see the change in energy along the detector
                          //event_MIP->Fill(FEB[FEBs[i]].hitCharge_pe->at(check));
                          // record sum of measured charge for each layer in z direction for this event
                          energyDepZY[MapCon[FEBs[i]][0][hitChannel]] += hitChargeCalib;
                        }
                      } 
                      // top and bottom: plot event
                      else {
                        if (hitCharge < 10000){
                          event_XZ->Fill(
                              MapCon[FEBs[i]][0][hitChannel], // x pos
                              MapCon[FEBs[i]][1][hitChannel], // z pos
                              hitCharge
                              );
                          event_Time_XZ->Fill(triggerTime - calibTime);

                          if (argv1 == "avg"){
                            if (hitCharge > chargeMin){
                              topTimeRef += calibTime;
                              tTRcount += 1;
                            }
                          }
                          if (argv1 == "median" || argv1 == "mode"){
                            if (hitCharge > chargeMin){
                              topTimeRefVec.push_back(calibTime); // round time to 2 dp so we can find a mode
                              topTimeRefZ.push_back(MapCon[FEBs[i]][1][hitChannel]);
                            }
                          }
                          // fill MIP plot with deposited charge (dE/dx) data
                          //event_MIP->Fill(FEB[FEBs[i]].hitCharge_pe->at(check));
                          // record sum of measured charge for each layer in z direction for this event
                          energyDepXZ[MapCon[FEBs[i]][1][hitChannel]] += hitChargeCalib;
                        }
                      }
                      // get hit times from front channels and record earliest time
                      //cout << timeRef << endl;
                      //if (FEBs[i] == 3 || FEBs[i] == 9){ // top/bottom channels at z=0
                      //  if (calibTime > 0){
                      //    if (MapCon[FEBs[i]][1][hitChannel] == 0){
                      //      if (argv1 == "front"){
                      //        if (FEBs[i] != 0 || FEBs[i] != 16){
                      //          if (FEBs[i] == 1 || FEBs[i] == 2 || FEBs[i] == 17 || FEBs[i] == 24){
                      //            if (sideTimeRef > calibTime || sideTimeRef == 0)
                      //              sideTimeRef = calibTime;
                      //          }
                      //          else{
                      //            if (topTimeRef > calibTime || topTimeRef == 0)
                      //              topTimeRef = calibTime;
                      //          }
                      //        }
                      //      }
                      //      if (earlyTime > calibTime || earlyTime == 0){
                      //        if (hitCharge > chargeMin){
                      //          earlyTime = calibTime;
                      //        }
                      //      }
                      //    }
                      //  }
                      //}

                      //// get hit times from back channels and record latest time
                      ////cout << timeRef << endl;
                      //if (FEBs[i] == 20 || FEBs[i] == 27){ // top/bottom channels at z=47
                      //  if (calibTime > 0){
                      //    if (MapCon[FEBs[i]][1][hitChannel] == 47){
                      //      if (lateTime < calibTime || lateTime == 0){
                      //        if (hitCharge > chargeMin)
                      //          lateTime = calibTime;
                      //      }
                      //    }
                      //  }
                      //}
                    }
                  //}
                }
              }
            }
          }
          // fill TH1Fs with summed light yield from z layers
          for (int ij = 0; ij < 48; ij++ ){
            energyDep[ij] = energyDepZY[ij] + energyDepXZ[ij];
            event_LY->Fill(ij,energyDep[ij]);
            event_LY_XZ->Fill(ij,energyDepXZ[ij]);
            event_LY_ZY->Fill(ij,energyDepZY[ij]);
          }
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

          // only take events that look go straight through the detector with no interactions
          if (muon_cut(event_XZ, event_ZY, event_XY) == 0){
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
              //refMap_ZY->Fill(sideTimeRefZ[median.elem[0]]);
              median = findMedian(topTimeRefVec);
              topTimeRef = median.result;
              //refMap_XZ->Fill(topTimeRefZ[median.elem[0]]);
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
                //refMap_ZY->Fill(sideTimeRefZ[midmode.result]);
              }
              //cout << "test" << endl;
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
         
            // iterate through FEBs again to fill display showing all selected events
            for (int i = 0; i < 19; i++){
              if (FEBs[i]!=12){
                // report if no reference time was found for selected event
                //if (timeRef == 0){
                //  cout << "NO REF" << endl;
                //  cout << subSpill+1 << " " << numberEvents+1 << endl;
                //}
                // get lower and upper positions of GTrigTags that are equal to the GTrigTag of the TOF trigger, i.e. get signals from a single event
                auto bounds=std::equal_range (FEB[FEBs[i]].GTrigTag->begin(), FEB[FEBs[i]].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
                GTindex[0] = bounds.first - FEB[FEBs[i]].GTrigTag->begin(); // element no. of first GTrigTag that equals trigger's
                GTindex[1] = bounds.second - FEB[FEBs[i]].GTrigTag->begin(); // element no. of last GTrigTag that equals trigger's
                // iterate through elements of a single event
                for (int check = GTindex[0]; check < GTindex[1]; check++){
                  Double_t triggerTime = FEB[12].hitTimefromSpill->at(TOFtrigger)*2.5; //convert to ns
                  Int_t hitCharge = FEB[FEBs[i]].hitCharge_pe->at(check);
                  Int_t HGCharge = FEB[FEBs[i]].hitHG_pe->at(check);
                  Double_t hitTime = FEB[FEBs[i]].hitTimefromSpill->at(check)*2.5;
                  Int_t hitChannel = FEB[FEBs[i]].hitsChannel->at(check);
                  Double_t travelDist = 0;
                  Double_t calibTime = 0;
                  Int_t pairCharge = 0;
                  bool keepHit = false;

                  // find distance light travelled to reach MPPC
                  if (FEBs[i] == 1 || FEBs[i] == 2 || FEBs[i] == 17 || FEBs[i] == 24){ // side channels
                    for (int j=0; j<28; j++){
                      if (FEBpairs[j] == FEBs[i]){
                        Int_t pair = j;
                        auto bounds2=std::equal_range(FEB[pair].GTrigTag->begin(), FEB[pair].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
                        GTindex2[0] = bounds2.first - FEB[pair].GTrigTag->begin(); 
                        GTindex2[1] = bounds2.second - FEB[pair].GTrigTag->begin();
                        for (int check2 = GTindex2[0]; check2 < GTindex2[1]; check2++){
                          if (MapCon[FEBs[i]][0][hitChannel] == MapCon[pair][1][(int)FEB[pair].hitsChannel->at(check2)]){
                            if (FEB[pair].hitCharge_pe->at(check2) > pairCharge){
                              pairCharge = FEB[pair].hitCharge_pe->at(check2);
                              keepHit=true;
                              if (FEBs[i] == 1 || FEBs[i] == 24) // right channels
                                travelDist = (24.-MapCon[pair][0][(int)FEB[pair].hitsChannel->at(check2)]-0.5)*0.01; // meters
                              else // left channels
                                travelDist = (MapCon[pair][0][(int)FEB[pair].hitsChannel->at(check2)]+0.5)*0.01; // meters
                            }
                          }
                        }
                      }
                    }
                  }
                  else if (FEBs[i] != 0 && FEBs[i] != 16) { // top/bottom channels
                    Int_t pair = FEBpairs[FEBs[i]];
                    auto bounds2=std::equal_range(FEB[pair].GTrigTag->begin(), FEB[pair].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
                    GTindex2[0] = bounds2.first - FEB[pair].GTrigTag->begin(); 
                    GTindex2[1] = bounds2.second - FEB[pair].GTrigTag->begin();
                    for (int check2 = GTindex2[0]; check2 < GTindex2[1]; check2++){
                      if (MapCon[FEBs[i]][1][hitChannel] == MapCon[pair][0][(int)FEB[pair].hitsChannel->at(check2)]){
                        if (FEB[pair].hitCharge_pe->at(check2) > pairCharge){
                          pairCharge = FEB[pair].hitCharge_pe->at(check2);
                          keepHit=true;
                          if (FEBs[i] == 3 || FEBs[i] == 4 || FEBs[i] == 8 || FEBs[i] == 18 || FEBs[i] == 19 || FEBs[i] == 20) // top channels
                            travelDist = (8.-MapCon[pair][1][(int)FEB[pair].hitsChannel->at(check2)]-0.5)*0.01; // meters
                          else // bottom channels
                            travelDist = (MapCon[pair][1][(int)FEB[pair].hitsChannel->at(check2)]+0.5)*0.01; // meters
                        }
                      }
                    }
                  }
                  else keepHit=true;
                  // account for light attenuation with hit charge
                  Int_t hitChargeCalib = hitCharge*exp(-travelDist/attenLen);
                  // calibrate hit time to account for time-walk effect, fibre travel time and light attenuation
                  if (FEBs[i] < 12)
                    element = i*96*5+hitChannel*5;
                  else
                    element = (i-1)*96*5+hitChannel*5;
                  if (FEBs[i] == 8 || FEBs[i] == 11)
                    calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(-hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
                  else
                    calibTime = hitTime-TimeOffsets[element+1]*exp(-TimeOffsets[element+2]*hitCharge)-TimeOffsets[element+3]/(hitCharge+TimeOffsets[element+4])-(travelDist/fibreSpeed)*1e9;
                  //calibTime = hitTime;

                  // cut on time window for muons (and electrons and pions)
                  if (triggerTime - calibTime < 300
                        && triggerTime - calibTime > 230){
                    //if (hitCharge > chargeMin){
                      if (keepHit) {
                        //if (FEBs[i]<12)
                        //  AmpTime[i*96+hitChannel]->Fill(hitChargeCalib,triggerTime-calibTime,1);
                        //else
                        //  AmpTime[(i-1)*96+hitChannel]->Fill(hitChargeCalib,triggerTime-calibTime,1);
                        // front and back: plot event using mapping info 
                        if (FEBs[i] == 0 || FEBs[i] == 16){
                          if (hitCharge > chargeMin && hitCharge < 10000){
                            // fill map showing all events
                            EventsMap_XY->Fill(
                                MapCon[FEBs[i]][0][hitChannel], // x position
                                MapCon[FEBs[i]][1][hitChannel], // y position
                                1);
                          }
                        } 
                        // sides: plot event
                        else if (FEBs[i] == 1 || FEBs[i] == 2 || FEBs[i] == 17 || FEBs[i] == 24){
                          if (hitCharge > chargeMin && hitCharge < 10000){
                            EventsMap_ZY->Fill(
                                MapCon[FEBs[i]][0][hitChannel], // z pos
                                MapCon[FEBs[i]][1][hitChannel], // y pos
                                1);
                            // fill all events MIP plot with deposited charge (dE/dx) data
                            //EventsMIP_ZY->Fill(FEB[FEBs[i]].hitCharge_pe->at(check));
                            // fill time res plot for relevant channels, but only if a time reference was able to be found and if the hit does not look like cross-talk
                            if (FEBs[i]<12){
                              if (argv1 == "trigger")
                                Time_res->Fill(-calibTime + triggerTime);
                              else {
                                if (sideTimeRef != 0)
                                  Time_res->Fill(calibTime - sideTimeRef);
                              }
                            }
                            else {
                              if (argv1 == "trigger")
                                Time_res->Fill(-calibTime + triggerTime);                   
                              else {
                                if (sideTimeRef != 0)
                                  Time_res->Fill(calibTime - sideTimeRef);
                              }
                            }
                          }
                        } 
                        // top and bottom: plot event
                        else if ( FEBs[i] != 18) {
                        //else {
                          if (hitCharge > chargeMin && hitCharge < 10000){
                            EventsMap_XZ->Fill(
                                MapCon[FEBs[i]][0][hitChannel], // x pos
                                MapCon[FEBs[i]][1][hitChannel], // z pos
                                1);
                            // fill time res plot for relevant channels, but only if a time reference was able to be found and if the hit does not look like cross-talk
                            if (topTimeRef != 0){
                              if (FEBs[i]<12){
                                if (argv1 == "trigger")
                                  Time_res->Fill(-calibTime + triggerTime);
                                else {
                                  if (topTimeRef != 0)
                                    Time_res->Fill(calibTime - topTimeRef);
                                }
                              }
                              else {
                                if (argv1 == "trigger")
                                  Time_res->Fill(-calibTime + triggerTime);                   
                                else {
                                  if (topTimeRef != 0)
                                    Time_res->Fill(calibTime - topTimeRef);
                                }
                              }
                            }
                            // fill MIP plot with deposited charge (dE/dx) data
                            //EventsMIP_XZ->Fill(FEB[FEBs[i]].hitCharge_pe->at(check));
                          }
                        }
                      } // lone hit cut
                    //} // charge cut
                  } // time cut
                } // iterate through events
              } // no FEB12
            } // iterate through FEBs
            //if (earlyTime!=0 && lateTime!=0)
            //  TOF->Fill((lateTime-earlyTime)*2.5);
            
            //###########
            // create event displays
            //###########
            eventsPassed += 1; 
            if (eventNum<NumberEvDis){       
 
              numberEvents += 1;
              c1->Clear();
     
              // plot the 4 event displays on one canvas 
              c1->Divide(3,3);

              // customise stats box
              gStyle->SetStatX(0.3);
              gStyle->SetStatY(1);
              gStyle->SetStatH(0.1155555);
              gStyle->SetStatW(0.2);
              gStyle->SetOptStat("rme");
              gStyle->SetTitleY(0.95);
              gStyle->SetPadRightMargin(0.2);
              gStyle->SetPadTopMargin(0.15);
              c1->Update();

              c1->cd(1);
              event_XY->Draw("colorz");
      
              c1->cd(2);
              event_ZY->Draw("colorz");

              c1->cd(3);
              event_XZ->Draw("colorz");
      
              c1->cd(4);
              event_LY->Draw("HIST");

              c1->cd(5);
              event_LY_ZY->Draw("HIST");

              c1->cd(6);
              event_LY_XZ->Draw("HIST");

              c1->cd(7);
              event_Time->Draw("HIST");

              c1->cd(8);
              event_Time_ZY->Draw("HIST");

              c1->cd(9);
              event_Time_XZ->Draw("HIST");

              c1->cd(7);
              TLatex *text = new TLatex();
              text->SetNDC();
              //if (earlyTime!=0 && lateTime!=0){
              //  text->DrawLatex(0,0,("TOF = "+to_string((lateTime-earlyTime)*2.5)+" #pm 2.5 ns").c_str());
              //}
              //else
              //  text->DrawLatex(0,0,"TOF Unavailable");            

              c1->Update();
              events2D -> cd();

              // save the canvas in events2D folder in output file
              c1->Write();

              //c2->cd();
              //event_MIP->Draw();

              //c2->Update();
               //eventsMIP->cd();
              //c2->Write();
        
              previousTime = 2.5*FEB[12].hitTimefromSpill->at(TOFtrigger);
              eventNum++;
            }
          }
          //##########
          //##########

          delete event_XY;
          delete event_ZY;
          delete event_XZ;
          delete event_LY;
          //delete event_MIP;
          delete event_LY_XZ;
          delete event_LY_ZY;
          delete event_Time;
          delete event_Time_XZ;
          delete event_Time_ZY;
      
        } // if TOF timediff positive
      } // iterate through TOFtrigger
    cout << numberEvents << " events plotted from spill" << endl;
    } // if spillMissed
  } // iterate through spills
  // iterate through side channels to find minimum entry number for time res
  cout << "Events passed: " << eventsPassed << endl;
  int n = 0;
  n = Time_res->GetEntries();
  wfile.cd();
  EventsMap_XY->Write();
  EventsMap_ZY->Write();
  EventsMap_XZ->Write();
  //TOF->Write();
  //refMap_XZ->Write();
  //refMap_ZY->Write();
  //EventsMIP_XZ->Write();
  //EventsMIP_ZY->Write();
  // create function to use in Time Res fit
  //TF1 *f1 = new TF1("f1", "gaus", -20, 15);
  if (argv1 == "trigger")
    TF1 *f1 = new TF1("f1", "gaus", 260, 330);
  else
    TF1 *f1 = new TF1("f1", "gaus", -20, 15);
  c3->cd();
  // fit time res plot and record time resolution
  TFitResultPtr r = Time_res->Fit("f1","QSR");
  gStyle->SetOptFit(1111);
  Time_res->Draw();
  c3->Write();
  wfile.Close();
  FileInput->Close();
  return 0;
}
