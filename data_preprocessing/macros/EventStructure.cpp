

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

ClassImp(Hit);
ClassImp(Event);

int main(int argc, char **argv) {

  // process arguments
  if (argc != 3) { printf("Enter source _calib.root file and data type   ./EventStructure inputfile datatype\n\n Datatypes that dictate event time window\n CERN: CERN beam test data\n LANL20: LANL beam test data from 20 m location\n LANL90: LANL beam test data from 90 m location\n"); return EXIT_FAILURE; }
  string sFileName(argv[1]);
  string sDataType(argv[2]);
  int timeWindowLow;
  int timeWindowHigh;
  bool USJ = false;
  if (CheckFile(sFileName) == 1) { printf("Source must be a _calib.root file \n"); return EXIT_FAILURE; }
  if (sDataType == "CERN"){
    timeWindowLow = -130;
    timeWindowHigh = -50;
  }
  else if (sDataType == "LANL20"){
    timeWindowLow = -288;
    timeWindowHigh = 418;
  }
  else if (sDataType == "LANL90"){
    timeWindowLow = -358;
    timeWindowHigh = 340;
  }  else if (sDataType == "USJ"){
    timeWindowLow = -358;
    timeWindowHigh = 1058;
    //timeWindowLow = -358;
    //timeWindowHigh = 340;
    USJ = true;
  }
  //else if (sDatatype == "cosmics"){
  //  timeWindowLow = -1e9;
  //  timeWindowHigh = 1e9;    
  //}
  else{ printf("Enter source _calib.root file and data type   ./EventStructure inputfile datatype\n\nDatatypes that dictate event time window\n CERN: CERN beam test data\n LANL20: LANL beam test data from 20 m location\n LANL90: LANL beam test data from 90 m location\n"); return EXIT_FAILURE; }

  int NumberOfEB = 30;
  int len = 18;
  vector<int> FEBs = {0,1,2,3,8,9,10,11,12,16,17,18,19,20,24,25,26,27};
  //vector<int> FEBs = {0,1,2,3,4,8,9,10,11,16,17,18,19,24,25,26,27};
  //int FEBs[len];
  //memset(FEBs,0,len*sizeof(int));
  if (USJ){
    //NumberOfEB = 65;
    //len = 7;
    //FEBs.clear();
    //FEBs = {12,56,57,58,59,60,61};

    NumberOfEB = 65;
    len = 25;
    FEBs.clear();
    FEBs = {0,1,2,3,4,8,9,10,11,12,16,17,18,19,20,24,25,26,27,56,57,58,59,60,61};
  }

  TFile *FileInput=new TFile(sFileName.c_str(),"read");
  cout<<"Reading "<<sFileName<<endl;

  string rootFileOutput=GetLocation (sFileName.c_str());
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
      std::cout<<sFEB<<'\t';

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

      cout<< "Number of Spills "<< TotalSpills - 1 <<endl;
      NumofSpills[febx]=TotalSpills;
      if (MinSpills > TotalSpills) MinSpills = TotalSpills; //if the FEBs contain different numbers of spills, use the smallest number
      febx++;
    }
  }

  TotalSpills = MinSpills;
  cout<<"Processing "<<TotalSpills-1<<" spills"<<endl;

  //Check that the number of spills is consistent across all FEBs

  int RecSpills = FEBtree[12]->GetEntries();
  for (int ih=0; ih<len; ih++){
    if (NumofSpills[ih]==0){
      cout<<"FEB "<<FEBs[ih]<<" missing!!"<<endl;
      return EXIT_FAILURE;
    }
    //else if (NumofSpills[ih]!=RecSpills){
    //  cout<<"Number of spills is not consistent!!"<<endl;
    //  return EXIT_FAILURE;
    //}
  }

  int MapCon[62][2][96];
  for (int iFEB = 0; iFEB<len; iFEB++) {
    if (FEBs[iFEB] != 12){
      sFEBnum.str("");
      sFEBnum << FEBs[iFEB];
      sFEB = "../mapping/" + sFEBnum.str() + ".txt";
      ifstream fmap(sFEB.c_str());
      int temp=0;
      while (temp<96) {
        fmap >> temp >> MapCon[FEBs[iFEB]][0][temp] >>MapCon[FEBs[iFEB]][1][temp];
        temp++;
      }
      fmap.close();
    }
  }

  std::cout<<"mapping finished .. "<<std::endl;
  //AllEvents tree
  TTree AllEvents("AllEvents", "The ROOT tree of events");
  Event* event = new Event();
  AllEvents.Branch("Event", "Event", event);

  Int_t eventNum=0;
  bool SpillMissed = false;
  Int_t XZocc[48], ZYocc[48], EDepXZ[48],EDepZY[48],EDep[48];
  int hitCount = 0;

  //loop over spills
  for (Int_t subSpill = 0; subSpill<TotalSpills-1; subSpill++) {
  //for (Int_t subSpill = 0; subSpill<3000; subSpill++) {
    Int_t micropulse=0;
    FEBtree[12]->GetEntry(subSpill);
    if(FEB[12].SpillTag->empty()){
      cout<<"EMPTY SPILL"<<endl;
      continue;
    }

    Int_t Spilltag = FEB[12].SpillTag->back();
    //cout << "_Getting Spill Number " << Spilltag  << " of "<<TotalSpills-1<<"...";

    //loop over FEBs
    for (int ik = 0; ik < len; ik++){

      // cout << FEBs[ik] << endl;
      FEBtree[FEBs[ik]]->GetEntry(subSpill);
      if (FEB[FEBs[ik]].SpillTag->size() == 0){
        continue;
      }
      // cout << FEB[FEBs[ik]].SpillTag->back() << endl;
      if (FEB[FEBs[ik]].SpillTag->back() != Spilltag)
      {
        cout << "ERROR (SpillTag != FEB12SpillTag)"<<"FEB "<<FEBs[ik]<<" spill tag = "<<FEB[FEBs[ik]].SpillTag->back()<<endl;
        SpillMissed = true;
        break;
      }
      if (FEB[FEBs[ik]].SpillTag->size() < 2 ){
        cout << "NULL (Less than 2 hits)"<<endl;
        SpillMissed = true;
        break;
      } else {
        SpillMissed = false;
      }
    }

    if (!SpillMissed)
    {
      //loop over hits
      for ( int TOFtrigger = 0; TOFtrigger < FEB[12].FEBSN->size(); TOFtrigger++){
        //cout<<"size : "<<FEB[12].FEBSN->size()<< "  now at : "<<TOFtrigger<<endl;
        //cout <<(int)(TOFtrigger*100/(Double_t)FEB[12].FEBSN->size())<<"% done\r"<<flush;
        if (sDataType == "LANL20" || sDataType == "LANL90" || sDataType == "USJ"){
          if ( FEB[12].hitsChannel->at(TOFtrigger) != 1 && FEB[12].hitsChannel->at(TOFtrigger) != 0) {
            continue;
	  }
        }
        if (FEB[12].hitTimeDif->at(TOFtrigger) > 0){
          for (int ij = 0; ij < 48; ij++ ){
            XZocc[ij] = 0;
            ZYocc[ij] = 0;
            EDepXZ[ij] = 0;
            EDepZY[ij] = 0;
            EDep[ij] = 0;
          }
          TH2D *hitsXY = new TH2D("hitsXY","", 24,0,24, 8,0,8);

          Int_t GTindex[2] = {0,0};
          for (int i = 0; i < len; i++){
            if (FEBs[i]!=12){
              auto bounds=std::equal_range (FEB[FEBs[i]].GTrigTag->begin(), FEB[FEBs[i]].GTrigTag->end(), FEB[12].GTrigTag->at(TOFtrigger));
              GTindex[0] = bounds.first - FEB[FEBs[i]].GTrigTag->begin();
              GTindex[1] = bounds.second - FEB[FEBs[i]].GTrigTag->begin();

              for (int check = GTindex[0]; check <  GTindex[1]; check++){
                //bool passMycut = false;      
		//cout<<"dt : "<<i<<" "<<check<<" "<<FEB[FEBs[i]].hitTimefromSpill->at(check) - FEB[12].hitTimefromSpill->at(TOFtrigger)<<endl;
                //if (FEB[FEBs[i]].hitTimefromSpill->at(check) - FEB[12].hitTimefromSpill->at(TOFtrigger) > timeWindowLow && FEB[FEBs[i]].hitTimefromSpill->at(check) - FEB[12].hitTimefromSpill->at(TOFtrigger) < timeWindowHigh )
		{
                  Hit* hit = event->AddHit();

                  if (FEB[FEBs[i]].hitCharge_pe->at(check)<10000) hit->SetPE(FEB[FEBs[i]].hitCharge_pe->at(check));
                  hit->SetHG_pe(FEB[FEBs[i]].hitHG_pe->at(check));
                  hit->SetLG_pe(FEB[FEBs[i]].hitLG_pe->at(check));
                  hit->SetToT_pe(FEB[FEBs[i]].hitToT_pe->at(check));
                  hit->SetHG_ADC(FEB[FEBs[i]].hitAmpl->at(check));
                  hit->SetLG_ADC(FEB[FEBs[i]].hitLGAmpl->at(check));
                  hit->SetRE(FEB[FEBs[i]].hitLeadTime->at(check));
                  hit->SetFE(FEB[FEBs[i]].hitTrailTime->at(check));
                  hit->SetToT(FEB[FEBs[i]].hitTimeDif->at(check));
                  hit->SetDt(FEB[FEBs[i]].hitTimefromSpill->at(check) - FEB[12].hitTimefromSpill->at(TOFtrigger));
                  hit->SetSpillTag(FEB[FEBs[i]].SpillTag->at(check));
                  hit->SetSpillTime(FEB[FEBs[i]].SpillTime->at(check));
                  hit->SetSpillTrailTime(FEB[FEBs[i]].SpillTrailTime->at(check));
                  hit->SetTfromSpill(FEB[FEBs[i]].hitTimefromSpill->at(check));
                  hit->SetFEB(FEBs[i]);
                  hit->SetCh((int)FEB[FEBs[i]].hitsChannel->at(check));
                  hit->SetGTrigTag(FEB[FEBs[i]].GTrigTag->at(check));
                  hit->SetGTrigTime(FEB[FEBs[i]].GTrigTime->at(check));
		  hitCount ++;
		  //cout<<"hit # "<<hitCount<<endl;

                  if ( FEBs[i] == 0 || FEBs[i] == 16 || (FEBs[i] == 56 && FEB[FEBs[i]].hitsChannel->at(check) < 64)){
                    hit->SetView(0);    //XY view
                    hit->SetX(MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]);
                    hit->SetY(MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]);
                    hit->SetZ(-1);
		    //cout<<"a XY hit : "<<MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]<<" "<<MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]<<" "<<FEB[FEBs[i]].hitCharge_pe->at(check)<<endl;
                    if (hit->GetPE() > 20)
                      hitsXY->Fill(hit->GetX(),hit->GetY());

                  } else if ( FEBs[i] == 1 || FEBs[i] == 2 || FEBs[i] == 17 || FEBs[i] == 24 || FEBs[i] == 56 || FEBs[i] == 57 || (FEBs[i] == 60 && FEB[FEBs[i]].hitsChannel->at(check) < 64) || (FEBs[i] == 61 && FEB[FEBs[i]].hitsChannel->at(check) > 31)){
                    hit->SetView(2);   //ZY view
                    hit->SetX(-1);
		    //cout<<"a ZY hit : "<<MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]<<" "<<MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]<<" "<<FEB[FEBs[i]].hitCharge_pe->at(check)<<endl;
                    hit->SetZ(MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]);
                    hit->SetY(MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]);
                    if (FEB[FEBs[i]].hitCharge_pe->at(check)>0 && FEB[FEBs[i]].hitCharge_pe->at(check)<10000) {
                      ZYocc[MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]]+=1;
                      EDepZY[MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]]+=FEB[FEBs[i]].hitCharge_pe->at(check);
                    }

                  } else {
                    hit->SetView(1);   //XZ view
                    hit->SetX(MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]);
                    hit->SetY(-1);
		    //cout<<"a XZ hit : "<<MapCon[FEBs[i]][0][(int)FEB[FEBs[i]].hitsChannel->at(check)]<<" "<<MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]<<" "<<FEB[FEBs[i]].hitCharge_pe->at(check)<<endl;
                    hit->SetZ(MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]);
                    if (FEB[FEBs[i]].hitCharge_pe->at(check)>0 && FEB[FEBs[i]].hitCharge_pe->at(check)<10000) {
                      XZocc[MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]]+=1;
                      EDepXZ[MapCon[FEBs[i]][1][(int)FEB[FEBs[i]].hitsChannel->at(check)]]+=FEB[FEBs[i]].hitCharge_pe->at(check);
                    }
                  }
                }
              }
            }
          }



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

          Double_t xSD = hitsXY->GetStdDev(1);
          Double_t ySD = hitsXY->GetStdDev(2);
          event->SetFEB12ch(FEB[12].hitsChannel->at(TOFtrigger));
          event->SetFEB12LeadTime(FEB[12].hitLeadTime->at(TOFtrigger));
          event->SetFEB12hitTfromSpill(FEB[12].hitTimefromSpill->at(TOFtrigger));
          event->SetMicropulse(micropulse);
          event->SetRange(range);
          event->SetMaxCharge(event->FindMaxCharge());
          event->SetEventID(eventNum);
          event->SetEmptyZ(empty);
          event->SetXStdDev(xSD);
          event->SetYStdDev(ySD);

          AllEvents.Fill();

          delete hitsXY;
          event->Clear();
          eventNum++;
          micropulse++;

        }
      }
      //cout <<"_Getting Spill Number " << Spilltag  << " of "<<TotalSpills-1<<"..."<<"100% done"<<endl;
    }
    //cout<<"subspill & micropulse number & eventNum : "<<subSpill<<" "<<micropulse<<" "<<eventNum<<endl;
  }

  cout<<"eventNum : "<<eventNum<<endl;

  wfile.cd();
  AllEvents.Write("",TObject::kOverwrite);
  //AllEvents.Delete();
  wfile.Close();
  FileInput->Close();


  return 0;
}
