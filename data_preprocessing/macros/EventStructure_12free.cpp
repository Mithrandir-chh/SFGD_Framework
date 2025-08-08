#include "../src/global_header.hh"
#include <cstring>

//Function to get the directory of the input file in order to create output file in same directory
string GetLocation(string str)
{
  int i = str.rfind("_calib");
  string way = str.substr(0,i);
  return way;
}

//Function to check that the input file is indeed a _all_calib.root file
bool CheckFile(string str)
{
  if (str.find("_calib.root") != string::npos) return 0;
  return 1;
}

struct vectorsTree
{
  vector<double> *FEBSN;
  vector<double> *SpillTag;
  vector<double> *GTrigTag;
  vector<double> *GTrigTime;
  vector<double> *hitsChannel;
  vector<double> *hitAmpl;
  vector<double> *hitAmplRec;
  vector<double> *hitLGAmpl;
  vector<double> *hitLGAmplRec;
  vector<double> *hitHG_pe;
  vector<double> *hitLG_pe;
  vector<double> *hitToT_pe;
  vector<double> *hitCharge_pe;
  vector<double> *hitLeadTime;
  vector<double> *hitTrailTime;
  vector<double> *hitTimeDif;
  vector<double> *SpillTime;
  vector<double> *SpillTimeGTrig;
  vector<double> *hitTimefromSpill;
  vector<double> *SpillTrailTime;
  vector<double> *AsicTemperature;
  vector<double> *FPGATemperature;
  vector<double> *GlobalHV;
  vector<double> *BoardTemperature;
  vector<double> *BoardHumidity;
};

// Structure to hold hit information for time grouping
struct HitInfo {
  double time;
  int feb;
  int hitIndex;
  bool used;
};

ClassImp(Hit);
ClassImp(Event);

// Function to add a hit to an event and update occupancy/energy arrays
bool AddHitToEvent(Event* event, vectorsTree FEB[], int febId, int hitIdx, 
                   int MapCon[62][2][96], Int_t XZocc[], Int_t ZYocc[], 
                   Int_t EDepXZ[], Int_t EDepZY[], TH2D* hitsXY, double peThreshold =20.0) {
  
  // Filter out hits below the threshold
  if (FEB[febId].hitHG_pe->at(hitIdx) < peThreshold) {
    // cout<<"Filtered out hit with hitHG_pe = "<<FEB[febId].hitHG_pe->at(hitIdx)<<endl;
    return false;
  }
  
  Hit* hit = event->AddHit();
  
  // Set hit properties
  if (FEB[febId].hitCharge_pe->at(hitIdx)<10000) 
    hit->SetPE(FEB[febId].hitCharge_pe->at(hitIdx));
  hit->SetHG_pe(FEB[febId].hitHG_pe->at(hitIdx));
  hit->SetLG_pe(FEB[febId].hitLG_pe->at(hitIdx));
  hit->SetToT_pe(FEB[febId].hitToT_pe->at(hitIdx));
  hit->SetHG_ADC(FEB[febId].hitAmpl->at(hitIdx));
  hit->SetLG_ADC(FEB[febId].hitLGAmpl->at(hitIdx));
  hit->SetRE(FEB[febId].hitLeadTime->at(hitIdx));
  hit->SetFE(FEB[febId].hitTrailTime->at(hitIdx));
  hit->SetToT(FEB[febId].hitTimeDif->at(hitIdx));
  hit->SetSpillTag(FEB[febId].SpillTag->at(hitIdx));
  hit->SetSpillTime(FEB[febId].SpillTime->at(hitIdx));
  hit->SetSpillTrailTime(FEB[febId].SpillTrailTime->at(hitIdx));
  hit->SetTfromSpill(FEB[febId].hitTimefromSpill->at(hitIdx));
  hit->SetFEB(febId);
  hit->SetCh((int)FEB[febId].hitsChannel->at(hitIdx));
  hit->SetGTrigTag(FEB[febId].GTrigTag->at(hitIdx));
  hit->SetGTrigTime(FEB[febId].GTrigTime->at(hitIdx));
  
  // Determine view and position based on FEB number
  if ( febId == 0 || febId == 16 || (febId == 56 && FEB[febId].hitsChannel->at(hitIdx) < 64)){
    hit->SetView(0);    //XY view
    hit->SetX(MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]);
    hit->SetY(MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]);
    hit->SetZ(-1);
    // cout<<"a XY hit : "<<MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]<<" "
        //<<MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]<<" "
        //<<FEB[febId].hitCharge_pe->at(hitIdx)<<endl;
    if (hit->GetPE() > 20 && hitsXY != nullptr)
      hitsXY->Fill(hit->GetX(),hit->GetY());

  } else if ( febId == 1 || febId == 2 || febId == 17 || febId == 24 || 
              febId == 56 || febId == 57 || 
              (febId == 60 && FEB[febId].hitsChannel->at(hitIdx) < 64) || 
              (febId == 61 && FEB[febId].hitsChannel->at(hitIdx) > 31)){
    hit->SetView(2);   //ZY view
    hit->SetX(-1);
    // cout<<"a ZY hit : "<<MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]<<" "
        //<<MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]<<" "
        //<<FEB[febId].hitCharge_pe->at(hitIdx)<<endl;
    hit->SetZ(MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]);
    hit->SetY(MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]);
    if (FEB[febId].hitCharge_pe->at(hitIdx)>0 && FEB[febId].hitCharge_pe->at(hitIdx)<10000) {
      ZYocc[MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]]+=1;
      EDepZY[MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]]+=FEB[febId].hitCharge_pe->at(hitIdx);
    }
  } else {
    hit->SetView(1);   //XZ view
    hit->SetX(MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]);
    hit->SetY(-1);
    hit->SetZ(MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]);
    // cout<<"a XZ hit : "<<MapCon[febId][0][(int)FEB[febId].hitsChannel->at(hitIdx)]<<" "
        //<<MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]<<" "
        //<<FEB[febId].hitCharge_pe->at(hitIdx)<<endl;
    if (FEB[febId].hitCharge_pe->at(hitIdx)>0 && FEB[febId].hitCharge_pe->at(hitIdx)<10000) {
      XZocc[MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]]+=1;
      EDepXZ[MapCon[febId][1][(int)FEB[febId].hitsChannel->at(hitIdx)]]+=FEB[febId].hitCharge_pe->at(hitIdx);
    }
  }
  
  return true;
}

// Function to finalize event properties
void FinalizeEvent(Event* event, Int_t XZocc[], Int_t ZYocc[], 
                   Int_t EDepXZ[], Int_t EDepZY[], TH2D* hitsXY, Int_t eventNum) {
  
  Int_t range = 0;
  Int_t empty = 0;
  for (int ij = 0; ij < 48; ij++ ){
    event->SetOccupancyXZ(ij,XZocc[ij]);
    event->SetOccupancyZY(ij,ZYocc[ij]);
    event->SetdEdzXZ(ij,EDepXZ[ij]);
    event->SetdEdzZY(ij,EDepZY[ij]);
    event->SetdEdz(ij,EDepXZ[ij]+EDepZY[ij]);
    if(EDepXZ[ij]+EDepZY[ij]>0) range = ij+1;
    if (EDepXZ[ij] <= 20 || EDepZY[ij] <= 20) empty++;
  }

  Double_t xSD = 0;
  Double_t ySD = 0;
  if (hitsXY != nullptr) {
    xSD = hitsXY->GetStdDev(1);
    ySD = hitsXY->GetStdDev(2);
  }
  
  event->SetRange(range);
  event->SetMaxCharge(event->FindMaxCharge());
  event->SetEventID(eventNum);
  event->SetEmptyZ(empty);
  event->SetXStdDev(xSD);
  event->SetYStdDev(ySD);
}

// Function to build spill-based events
void BuildSpillBasedEvents(TTree* FEBtree[], vectorsTree FEB[], vector<int>& FEBs, 
                          int MapCon[62][2][96], TTree& AllEvents, 
                          int TotalSpills, int len) {
  
  // cout << "Building spill-based events..." << endl;
  
  Event* event = new Event();
  AllEvents.Branch("Event", "Event", event);
  
  Int_t eventNum = 0;
  Int_t XZocc[48], ZYocc[48], EDepXZ[48], EDepZY[48], EDep[48];
  int hitCount = 0;
  
  // Loop over spills
  for (Int_t subSpill = 0; subSpill<TotalSpills-1; subSpill++) {
  // for (Int_t subSpill = 0; subSpill < 50; subSpill++) {  // Using 15 spills as in original
    
    // cout << "_Getting Spill Number " << subSpill  << " of "<<TotalSpills-1<<"...";
    
    // Clear occupancy and energy arrays
    for (int ij = 0; ij < 48; ij++ ){
      XZocc[ij] = 0;
      ZYocc[ij] = 0;
      EDepXZ[ij] = 0;
      EDepZY[ij] = 0;
      EDep[ij] = 0;
    }
    
    TH2D *hitsXY = new TH2D("hitsXY","", 24,0,24, 8,0,8);

    // Loop over FEBs
    for (int ik = 0; ik < len; ik++){
      // cout<<"check point 2."<<ik<<" "<<FEBs[ik]<<endl;

      FEBtree[FEBs[ik]]->GetEntry(subSpill);
      
      if (FEB[FEBs[ik]].SpillTag->size() == 0){
        continue;
      }
      
      // Process all hits for this FEB in this spill
      for (int ihit = 0; ihit < FEB[FEBs[ik]].hitAmpl->size(); ihit++){
        if (AddHitToEvent(event, FEB, FEBs[ik], ihit, MapCon, XZocc, ZYocc, 
                          EDepXZ, EDepZY, hitsXY)) {
          hitCount++;
          // cout<<"hit # "<<hitCount<<endl;
        }
      }
    }
    
    // Finalize event properties
    FinalizeEvent(event, XZocc, ZYocc, EDepXZ, EDepZY, hitsXY, eventNum);
    
    AllEvents.Fill();
    delete hitsXY;
    event->Clear();
    eventNum++;
  }
  
  // cout<<"Spill-based events: "<<eventNum<<endl;
  delete event;
}

// Function to build time-grouped events
void BuildTimeGroupedEvents(TTree* FEBtree[], vectorsTree FEB[], vector<int>& FEBs, 
                           int MapCon[62][2][96], TTree& TimeGroupedEvents, 
                           int TotalSpills, int len) {
  
  // cout << "Building time-grouped events..." << endl;
  // cout << "Using trigger threshold of 20pe and secondary threshold of 0pe" << endl;
  // cout << "Time window: dynamically expanded based on connected hits >20pe" << endl;
  
  Event* timeGroupedEvent = new Event();
  TimeGroupedEvents.Branch("Event", "Event", timeGroupedEvent);
  
  Int_t timeGroupedEventNum = 0;
  Int_t XZocc[48], ZYocc[48], EDepXZ[48], EDepZY[48], EDep[48];
  
  // Loop over spills
  for (Int_t subSpill = 0; subSpill<TotalSpills-1; subSpill++) {
  // for (Int_t subSpill = 0; subSpill < 50; subSpill++) {  // Using 15 spills as in original
    
    // cout<<"Starting time-based grouping for spill "<<subSpill<<endl;
    
    // Collect ALL hits from this spill (including those under 20pe)
    vector<HitInfo> allHitsInSpill;
    
    // Gather all hits from all FEBs for this spill
    for (int ik = 0; ik < len; ik++){
      FEBtree[FEBs[ik]]->GetEntry(subSpill);
      for (int ihit = 0; ihit < FEB[FEBs[ik]].hitAmpl->size(); ihit++){
        // Store ALL hits, we'll filter later based on context
        HitInfo hi;
        hi.time = FEB[FEBs[ik]].hitTimefromSpill->at(ihit);
        hi.feb = FEBs[ik];
        hi.hitIndex = ihit;
        hi.used = false;
        allHitsInSpill.push_back(hi);
      }
    }
    
    // Sort hits by time
    sort(allHitsInSpill.begin(), allHitsInSpill.end(), 
         [](const HitInfo& a, const HitInfo& b) { return a.time < b.time; });
    
    // cout<<"Total hits in spill: "<<allHitsInSpill.size()<<endl;
    
    // Group hits by time - look for seed hits > 20pe
    for (int i = 0; i < allHitsInSpill.size(); i++) {
      if (allHitsInSpill[i].used) continue; // Skip already grouped hits
      
      // Check if this hit can be a seed (>20pe)
      int febId = allHitsInSpill[i].feb;
      int hitIdx = allHitsInSpill[i].hitIndex;
      if (FEB[febId].hitHG_pe->at(hitIdx) <= 20) continue; // Not a valid seed
      
      // Found a seed hit > 20pe
      double seedTime = allHitsInSpill[i].time;
      // cout<<"New time group starting with seed hit >20pe at time "<<seedTime<<endl;
      
      // Clear arrays for new event
      for (int ij = 0; ij < 48; ij++ ){
        XZocc[ij] = 0;
        ZYocc[ij] = 0;
        EDepXZ[ij] = 0;
        EDepZY[ij] = 0;
        EDep[ij] = 0;
      }
      TH2D *hitsXY_time = new TH2D("hitsXY_time","", 24,0,24, 8,0,8);
      
      // First pass: find all connected hits > 20pe to determine the full time window
      vector<int> highEnergyHits;
      highEnergyHits.push_back(i);
      allHitsInSpill[i].used = true;
      
      double minTime = seedTime;
      double maxTime = seedTime;
      
      // Iteratively expand to find all connected high-energy hits
      bool expanded = true;
      while (expanded) {
        expanded = false;
        for (int j = 0; j < allHitsInSpill.size(); j++) {
          if (allHitsInSpill[j].used) continue;
          
          int febId_j = allHitsInSpill[j].feb;
          int hitIdx_j = allHitsInSpill[j].hitIndex;
          if (FEB[febId_j].hitHG_pe->at(hitIdx_j) <= 20) continue;
          
          double hitTime = allHitsInSpill[j].time;
          // Check if this high-energy hit is within 10 ticks of our current window
          if ((hitTime >= minTime - 10 && hitTime <= maxTime + 10)) {
            highEnergyHits.push_back(j);
            allHitsInSpill[j].used = true;
            if (hitTime < minTime) minTime = hitTime;
            if (hitTime > maxTime) maxTime = hitTime;
            expanded = true;
            // cout<<"Expanding window: found connected hit >20pe at time "<<hitTime<<endl;
          }
        }
      }
      
      // Now define the full window: 10 ticks before earliest and 10 ticks after latest
      double windowStart = minTime - 10;
      double windowEnd = maxTime + 10;
      // cout<<"Final time window: ["<<windowStart<<", "<<windowEnd<<"]"<<endl;
      
      // Add all high-energy hits to the event
      for (int idx : highEnergyHits) {
        int febId_h = allHitsInSpill[idx].feb;
        int hitIdx_h = allHitsInSpill[idx].hitIndex;
        AddHitToEvent(timeGroupedEvent, FEB, febId_h, hitIdx_h, MapCon, 
                      XZocc, ZYocc, EDepXZ, EDepZY, hitsXY_time);
      }
      
      // Second pass: add all hits > 0pe within the expanded window
      for (int j = 0; j < allHitsInSpill.size(); j++) {
        if (allHitsInSpill[j].used) continue; // Skip already used hits
        
        double hitTime = allHitsInSpill[j].time;
        if (hitTime >= windowStart && hitTime <= windowEnd) {
          int febId_j = allHitsInSpill[j].feb;
          int hitIdx_j = allHitsInSpill[j].hitIndex;
          
          // Try to add the hit with secondary threshold
          if (AddHitToEvent(timeGroupedEvent, FEB, febId_j, hitIdx_j, MapCon, 
                           XZocc, ZYocc, EDepXZ, EDepZY, hitsXY_time, 0.0)) {
            allHitsInSpill[j].used = true;
            double original_pe = FEB[febId_j].hitCharge_pe->at(hitIdx_j);
            // cout<<"Added secondary hit with "<<original_pe<<"pe at time "
                //<<allHitsInSpill[j].time<<endl;
          }
        }
      }
      
      // Finalize event properties
      FinalizeEvent(timeGroupedEvent, XZocc, ZYocc, EDepXZ, EDepZY, 
                   hitsXY_time, timeGroupedEventNum);
      
      // Only save if event has hits
      if (timeGroupedEvent->GetNHits() > 0) {
        TimeGroupedEvents.Fill();
        timeGroupedEventNum++;
        // cout<<"Time-grouped event "<<timeGroupedEventNum-1<<" has "
            //<<timeGroupedEvent->GetNHits()<<" hits"<<endl;
      }
      if (timeGroupedEvent->GetNHits()> 20){
        // cout<<"**********************************"<<endl;
        // exit(0);

      }

      delete hitsXY_time;
      timeGroupedEvent->Clear();
    }
  }
  
  // cout<<"Time-grouped events: "<<timeGroupedEventNum<<endl;
  delete timeGroupedEvent;
}

int main(int argc, char **argv) {

  // process arguments - simplified without data type requirement
  if (argc != 2) { 
    printf("Enter source _calib.root file: ./EventStructure inputfile\n"); 
    return EXIT_FAILURE; 
  }
  string sFileName(argv[1]);
  
  if (CheckFile(sFileName) == 1) { 
    printf("Source must be a _calib.root file \n"); 
    return EXIT_FAILURE; 
  }

  int NumberOfEB = 30;
  int len = 18;
  vector<int> FEBs = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};
  
  TFile *FileInput=new TFile(sFileName.c_str(),"read");
  cout<<"Reading "<<sFileName<<endl;

  string rootFileOutput=GetLocation(sFileName.c_str());
  rootFileOutput+="_events.root";

  TFile wfile(rootFileOutput.c_str(), "recreate");
  cout<<rootFileOutput<<endl;

  vectorsTree FEB[NumberOfEB];

  for (Int_t i=0;i<NumberOfEB;i++){
    FEB[i].FEBSN=0;
    FEB[i].SpillTag=0;
    FEB[i].hitsChannel=0;
    FEB[i].hitAmpl=0;
    FEB[i].hitLeadTime=0;
    FEB[i].GTrigTag=0;
    FEB[i].GTrigTime=0;
    FEB[i].hitLGAmpl=0;
    FEB[i].hitTrailTime=0;
    FEB[i].hitTimeDif=0;
    FEB[i].SpillTime=0;
    FEB[i].SpillTimeGTrig=0;
    FEB[i].hitTimefromSpill=0;
    FEB[i].SpillTrailTime=0;
    FEB[i].AsicTemperature=0;
    FEB[i].FPGATemperature=0;
    FEB[i].GlobalHV=0;
    FEB[i].BoardTemperature=0;
    FEB[i].BoardHumidity=0;
    FEB[i].hitHG_pe=0;
    FEB[i].hitLG_pe=0;
    FEB[i].hitToT_pe=0;
    FEB[i].hitCharge_pe=0;
    FEB[i].hitAmplRec=0;
    FEB[i].hitLGAmplRec=0;
  }

  TTree *FEBtree[NumberOfEB];

  int TotalSpills=0;
  int NumofSpills[len];
  int MinSpills=1e9;
  int febx=0;

  ostringstream sFEBnum;
  string sFEB;

  vector<int> FEBnumbers;
  FEBnumbers.clear();

  for (Int_t ih=0; ih<NumberOfEB; ih++) {
    sFEBnum.str("");
    sFEBnum << ih;
    sFEB = "FEB_"+sFEBnum.str();
    FEBtree[ih] = (TTree*)FileInput->Get(sFEB.c_str());
    if ((TTree*)FileInput->Get(sFEB.c_str())){
      // std::cout<<sFEB<<'\t';

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
      FEBtree[ih]->SetBranchAddress((sFEB+"_SpillTime").c_str(),&FEB[ih].SpillTime);
      FEBtree[ih]->SetBranchAddress((sFEB+"_SpillTimeGTrig").c_str(),&FEB[ih].SpillTimeGTrig);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitTimefromSpill").c_str(),&FEB[ih].hitTimefromSpill);
      FEBtree[ih]->SetBranchAddress((sFEB+"_SpillTrailTime").c_str(),&FEB[ih].SpillTrailTime);
      FEBtree[ih]->SetBranchAddress((sFEB+"_AsicTemperature").c_str(),&FEB[ih].AsicTemperature);
      FEBtree[ih]->SetBranchAddress((sFEB+"_FPGATemperature").c_str(),&FEB[ih].FPGATemperature);
      FEBtree[ih]->SetBranchAddress((sFEB+"_GlobalHV").c_str(),&FEB[ih].GlobalHV);
      FEBtree[ih]->SetBranchAddress((sFEB+"_BoardTemperature").c_str(),&FEB[ih].BoardTemperature);
      FEBtree[ih]->SetBranchAddress((sFEB+"_BoardHumidity").c_str(),&FEB[ih].BoardHumidity);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitAmplRecon").c_str(), &FEB[ih].hitAmplRec);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitLGAmplRecon").c_str(), &FEB[ih].hitLGAmplRec);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitHG_pe").c_str(), &FEB[ih].hitHG_pe);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitLG_pe").c_str(), &FEB[ih].hitLG_pe);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitToT_pe").c_str(), &FEB[ih].hitToT_pe);
      FEBtree[ih]->SetBranchAddress((sFEB+"_hitCharge_pe").c_str(), &FEB[ih].hitCharge_pe);

      TotalSpills = FEBtree[ih]->GetEntries();
      FEBtree[ih]->GetEntry(0);

      // cout<< "Number of Spills "<< TotalSpills - 1 <<endl;
      NumofSpills[febx]=TotalSpills;
      if (MinSpills > TotalSpills) MinSpills = TotalSpills;
      febx++;
    }
  }

  TotalSpills = MinSpills;
  // cout<<"Processing "<<TotalSpills-1<<" spills"<<endl;

  // cout<<"FEB number "<<len<<endl;

  int MapCon[62][2][96];
  for (int iFEB = 0; iFEB<len; iFEB++) {
    sFEBnum.str("");
    sFEBnum << FEBs[iFEB];
    sFEB = "../mapping/" + sFEBnum.str() + ".txt";
    ifstream fmap(sFEB.c_str());
    int temp=0;
    while (temp<96) {
      // cout<< FEBs[iFEB] << temp<<endl;
      fmap >> temp >> MapCon[FEBs[iFEB]][0][temp] >>MapCon[FEBs[iFEB]][1][temp];
      temp++;
    }
    fmap.close();
  }

  // std::cout<<"mapping finished .. "<<std::endl;
  
  // Create trees
  TTree AllEvents("AllEvents", "The ROOT tree of events");
  TTree TimeGroupedEvents("TimeGroupedEvents", "Time-grouped events within spills");

  // Build spill-based events
  //BuildSpillBasedEvents(FEBtree, FEB, FEBs, MapCon, AllEvents, TotalSpills, len);

  // Build time-grouped events
  BuildTimeGroupedEvents(FEBtree, FEB, FEBs, MapCon, TimeGroupedEvents, TotalSpills, len);

  // Write output
  wfile.cd();
  AllEvents.Write("",TObject::kOverwrite);
  TimeGroupedEvents.Write("",TObject::kOverwrite);
  wfile.Close();
  FileInput->Close();

  return 0;
}
