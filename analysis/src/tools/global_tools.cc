#include <TStyle.h>
#include <TColor.h>
#include <TH1.h>

TStyle* SetT2KStyle(Int_t WhichStyle = 1, TString styleName = "T2K") {
  TStyle *t2kStyle= new TStyle(styleName, "T2K approved plots style");

  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper

  Int_t FontStyle = 22;
  Float_t FontSizeLabel = 0.035;
  Float_t FontSizeTitle = 0.05;
  Float_t YOffsetTitle = 1.3;

  switch(WhichStyle) {
  case 1:
    FontStyle = 42;
    FontSizeLabel = 0.05;
    FontSizeTitle = 0.065;
    YOffsetTitle = 1.19;
    break;
  case 2:
    FontStyle = 42;
    FontSizeLabel = 0.035;
    FontSizeTitle = 0.05;
    YOffsetTitle = 1.6;
    break;
  case 3:
    FontStyle = 132;
    FontSizeLabel = 0.035;
    FontSizeTitle = 0.05;
    YOffsetTitle = 1.6;
    break;
  }

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetCanvasBorderSize(0);
  t2kStyle->SetFrameBorderSize(0);
  t2kStyle->SetDrawBorder(0);
  t2kStyle->SetTitleBorderSize(0);

  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  t2kStyle->SetFillColor(0);

  t2kStyle->SetEndErrorSize(4);
  t2kStyle->SetStripDecimals(kFALSE);

  t2kStyle->SetLegendBorderSize(0);
  t2kStyle->SetLegendFont(FontStyle);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20, 26);
  t2kStyle->SetPadTopMargin(0.1);
  t2kStyle->SetPadBottomMargin(0.15);
  t2kStyle->SetPadRightMargin(0.13); // 0.075 -> 0.13 for colz option
  t2kStyle->SetPadLeftMargin(0.16);//to include both large/small font options

  // Fonts, sizes, offsets
  t2kStyle->SetTextFont(FontStyle);
  t2kStyle->SetTextSize(0.08);

  t2kStyle->SetLabelFont(FontStyle, "x");
  t2kStyle->SetLabelFont(FontStyle, "y");
  t2kStyle->SetLabelFont(FontStyle, "z");
  t2kStyle->SetLabelFont(FontStyle, "t");
  t2kStyle->SetLabelSize(FontSizeLabel, "x");
  t2kStyle->SetLabelSize(FontSizeLabel, "y");
  t2kStyle->SetLabelSize(FontSizeLabel, "z");
  t2kStyle->SetLabelOffset(0.015, "x");
  t2kStyle->SetLabelOffset(0.015, "y");
  t2kStyle->SetLabelOffset(0.015, "z");

  t2kStyle->SetTitleFont(FontStyle, "x");
  t2kStyle->SetTitleFont(FontStyle, "y");
  t2kStyle->SetTitleFont(FontStyle, "z");
  t2kStyle->SetTitleFont(FontStyle, "t");
  t2kStyle->SetTitleSize(FontSizeTitle, "y");
  t2kStyle->SetTitleSize(FontSizeTitle, "x");
  t2kStyle->SetTitleSize(FontSizeTitle, "z");
  t2kStyle->SetTitleOffset(1.14, "x");
  t2kStyle->SetTitleOffset(YOffsetTitle, "y");
  t2kStyle->SetTitleOffset(1.2, "z");

  t2kStyle->SetTitleStyle(0);
  t2kStyle->SetTitleFontSize(0.06);//0.08
  t2kStyle->SetTitleFont(FontStyle, "pad");
  t2kStyle->SetTitleBorderSize(0);
  t2kStyle->SetTitleX(0.1f);
  t2kStyle->SetTitleW(0.8f);

  // use bold lines and markers
  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth( Width_t(2.5) );
  t2kStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  // -- color --
  // functions blue
  //t2kStyle->SetFuncColor(600-4);
  t2kStyle->SetFuncColor(2);
  t2kStyle->SetFuncWidth(2);

  t2kStyle->SetFillColor(1); // make color fillings (not white)
  // - color setup for 2D -
  // - "cold"/ blue-ish -
  Double_t red[]   = { 0.00, 0.00, 0.00 };
  Double_t green[] = { 1.00, 0.00, 0.00 };
  Double_t blue[]  = { 1.00, 1.00, 0.25 };
  // - "warm" red-ish colors -
  //  Double_t red[]   = {1.00, 1.00, 0.25 };
  //  Double_t green[] = {1.00, 0.00, 0.00 };
  //  Double_t blue[]  = {0.00, 0.00, 0.00 };

  Double_t stops[] = { 0.25, 0.75, 1.00 };
  const Int_t NRGBs = 3;
  const Int_t NCont = 500;

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  t2kStyle->SetNumberContours(NCont);

  // - Rainbow -
  //  t2kStyle->SetPalette(1);  // use the rainbow color set

  // -- axis --
  t2kStyle->SetStripDecimals(kFALSE); // don't do 1.0 -> 1
  //  TGaxis::SetMaxDigits(3); // doesn't have an effect
  // no supressed zeroes!
  t2kStyle->SetHistMinimumZero(kTRUE);


 return(t2kStyle);
}

void CenterHistoTitles(TH1 *thisHisto){
  thisHisto->GetXaxis()->CenterTitle();
  thisHisto->GetYaxis()->CenterTitle();
  thisHisto->GetZaxis()->CenterTitle();
}


int AddGridLinesToPad(TPad *thisPad) {
  thisPad->SetGridx();
  thisPad->SetGridy();
  return(0);
}


//Copied from Unpacking

//Function to get the directory of the input file in order to create output files in same directory
std::string GetDir(std::string str)
{
  int i = str.rfind("/");
  std::string way = str.substr(0,i);
  return way;
}

//Function to check that the input file is indeed a _NewStructure.root file
bool CheckFile(string str)
{
  if (str.find("_events") != string::npos) return 0;
  return 1;
}

// Function to get the name and directory of the source root file
std::string GetName(std::string str)
{
  int i = str.rfind("_events");
  std::string way = str.substr(0,i);
  return way;
}

// Setting The Event Selection Criteria
bool Selection (Event* event, TH2F* event_XY, array<Double_t, 24> dEdxXZ, array<Double_t, 8> dEdyZY)
{
  Double_t stdx = event_XY->GetStdDev(1);
  Double_t stdy = event_XY->GetStdDev(2);

  if (event->GetMaxCharge()>250 && event->GetdEdz(47)==0 && stdx<1 && stdy<1 && dEdxXZ[0]==0 && dEdxXZ[23]==0 && dEdyZY[0]<20 && dEdyZY[7]<20) return 0;
  return 1;
}

Hit* AddHit(TClonesArray* HitsArray, int HitNum){
    TClonesArray &hits = *HitsArray;
    Hit *hit = new(hits[HitNum++]) Hit();
    return hit;
}

void parseArguments(){
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


void linkFilesToTTrees(){

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

std::vector <ND280SFGDHit*> getEventMPPChits(int iev, int mode=0){

    dataIn->GetEntry(iev);
    bool jumpNextEvent = false;

    vector <ND280SFGDHit*> listOfHits;


    if(IsMC) {
        listOfHits = inputEvent->GetHits();
        for (auto hit : listOfHits){
            if(hit->GetView() == 0){
                hit->SetX(hit->GetX()-0.5+1.*SFGD_X/2);
                hit->SetY(hit->GetY()-0.5+1.*SFGD_Y/2);
            }
            if(hit->GetView() == 1){
                hit->SetX(hit->GetX()-0.5+1.*SFGD_X/2);
                hit->SetZ(hit->GetZ()-0.5+1.*SFGD_Z/2);
            }
            if(hit->GetView() == 2){
                hit->SetY(hit->GetY()-0.5+1.*SFGD_Y/2);
                hit->SetZ(hit->GetZ()-0.5+1.*SFGD_Z/2);
            }
        }
        if (IsReversed)
        {
            for (auto hit : listOfHits){
                if (hit->GetView() == 0)
                {
                    hit->SetX(SFGD_X - hit->GetX());
                }
                if (hit->GetView() == 1)
                {
                    hit->SetX(SFGD_X - hit->GetX());
                    hit->SetZ(SFGD_Z - hit->GetZ());
                }
                if (hit->GetView() == 2)
                {
                    hit->SetZ(SFGD_Z - hit->GetZ());
                }
                hit->SetY(hit->GetY()+0.5+1.*SFGD_Y/2);
                hit->SetZ(hit->GetZ()+0.5+1.*SFGD_Z/2);
            }
        }
    }
    else{
        TClonesArray * unpackHits = unpackEvent->GetHits();
        for(Int_t ihit=0; ihit<unpackEvent->GetNHits(); ihit++)
        {
            Hit *hit = (Hit*) unpackHits->At(ihit);
            //if(hit->GetDt() < -100 || hit->GetDt() > -60) continue;
            //if(hit->GetPE() <= 0 ) {jumpNextEvent = true; break;}
            ND280SFGDHit* sfgdHit = new ND280SFGDHit();
            sfgdHit->SetX(hit->GetX());
            sfgdHit->SetY(hit->GetY());
            sfgdHit->SetZ(hit->GetZ());
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
        if(jumpNextEvent) {listOfHits.clear();}
        //delete unpackHits;
    }
    return listOfHits;
}

void writeOutput(){

    //FileOutput->cd();
    FileOutput->Write("",TObject::kOverwrite);
    FileOutput->Close();
    FileInput->Close();

}

void handleEndOfExecution(){
    float time = ( clock() - start_time ) / CLOCKS_PER_SEC;
    cout << endl << std::setprecision(8) << "Elapsed time: " << time << endl;
    if(!batch){
        cout << "This macro remains open until it is manually closed with Ctrl+C\nallowing to watch the output on the output Canvas.\nTo avoid this, use option '-b'." << endl;
        return;
    }
    else{
        cout << "Execution completed successfully." << endl << endl;
        exit(1);
    }
}

void convertCoordinates(std::vector<ND280SFGDVoxel*> listOfVoxels){
  for(auto v:listOfVoxels){
    v->SetXYZ(v->GetX()+0.5+SFGD_X/2,v->GetY()+0.5+SFGD_Y/2,v->GetZ()+0.5+SFGD_Z/2);
    v->GetHits()[0]->SetX(v->GetX()+0.5+1.*SFGD_X/2); v->GetHits()[0]->SetY(v->GetY()+0.5+1.*SFGD_Y/2);
    v->GetHits()[1]->SetX(v->GetX()+0.5+1.*SFGD_X/2); v->GetHits()[1]->SetZ(v->GetZ()+0.5+1.*SFGD_Z/2);
    v->GetHits()[2]->SetY(v->GetY()+0.5+1.*SFGD_Y/2); v->GetHits()[2]->SetZ(v->GetZ()+0.5+1.*SFGD_Z/2);
  }
}

void convertVtxCoordiantes(ND280SFGDVoxel* v){
  v->SetXYZ(v->GetX()+0.5+SFGD_X/2,v->GetY()+0.5+SFGD_Y/2,v->GetZ()+0.5+SFGD_Z/2);
}

double att(double L){
    // computes the attenuation factor for a given distance L to the MPPC.
    const double a             = 0.77;
    const double LongAtt       = 4634.;
    const double ShortAtt      = 332.;
    const double d             = 0;//41; ¿?
    double factor = a*exp((-L-d)/LongAtt) + (1-a)*exp((-L-d)/ShortAtt);
    return factor;
}

void printLightInformation(std::vector<ND280SFGDVoxel*> voxels){
    for(auto v:voxels){
        std::cout << "truePE-recoPE: "  << std::left << std::setw(6) << v->GetTruePE() << "," << std::left << std::setw(6) << v->GetRecoPE() << "\t";
        std::cout << "multiplicity: ["  << std::left << std::setw(2) << v->GetHits()[0]->GetMultiplicity()           << "," << std::left << std::setw(2) << v->GetHits()[1]->GetMultiplicity()           << "," << std::left << std::setw(2) << v->GetHits()[2]->GetMultiplicity() << "]\t";
        std::cout << "PE: ["            << std::left << std::setw(4) << v->GetHits()[0]->GetPE()                     << "," << std::left << std::setw(4) << v->GetHits()[1]->GetPE()                     << "," << std::left << std::setw(4) << v->GetHits()[2]->GetPE() << "]\t";
        std::cout << "attenuation: ["   << std::left << std::setw(4) << v->GetHits()[0]->GetPE()/att(v->GetZ()*10)   << "," << std::left << std::setw(4) << v->GetHits()[1]->GetPE()/att(v->GetY()*10)   << "," << std::left << std::setw(4) << v->GetHits()[2]->GetPE()/att(v->GetX()*10) << "]\n\n";
    }
}

// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z, double x0=0, double y0=0, double z0=0) {
    // p[0],p[1],p[2], p[3],p[4],p[5] are free parameters
    if (x0 && y0 && z0){
        x = x0 + p[3]*t;
        y = y0 + p[4]*t;
        z = z0 + p[5]*t;
    }
    else{
        x = p[0] + p[3]*t;
        y = p[1] + p[4]*t;
        z = p[2] + p[5]*t;
    }
}

std::vector<ND280SFGDVoxel*> ApplyFakeNN(std::vector<ND280SFGDVoxel*> voxels){
  double prob [3] = {0.95,0.08,0.1};
  TRandom3* rndm = new TRandom3(0);
  std::vector<ND280SFGDVoxel*> selected;
  for (auto v:voxels) if (prob[v->GetTrueType()]>rndm->Uniform()) selected.push_back(v);
  delete rndm;
  return selected; 
}

std::vector<ND280SFGDVoxel*> cleanVoxels(std::vector<ND280SFGDVoxel*> voxels){
  std::vector<ND280SFGDVoxel*> selected;
  for (auto v:voxels) if (!v->GetTrueType()) selected.push_back(v);
  return selected; 
}

double computeMaxDistanceOfSets(std::vector<ND280SFGDVoxel*> v1, std::vector<ND280SFGDVoxel*> v2){
    double maxDist = -1; 
    for (int i=0; i<(int) v1.size(); ++i)
        for (int j=0; j<(int) v2.size(); ++j){
            double dist = v1[i]->DistToVoxel(v2[j]);
            if ( dist > maxDist) maxDist = dist;
        }
    return maxDist;
}

double computeMinDistanceOfSets(std::vector<ND280SFGDVoxel*> v1, std::vector<ND280SFGDVoxel*> v2){
    double minDist = 1E6; 
    for (int i=0; i<(int) v1.size(); ++i)
        for (int j=0; j<(int) v2.size(); ++j){
            double dist = v1[i]->DistToVoxel(v2[j]);
            if ( dist < minDist) minDist = dist;
        }
    return minDist;
}


//Function to store voxel information in a csv file
void DumpToCSVfile(std::ofstream &outCSVfile, ND280SFGDVoxelSet* EventVoxels, int iev, int opt=0){

    if(!(outCSVfile)) {std::cerr << "No filename specified!" << std::endl; exit(1);}

    for (auto n:EventVoxels->GetVoxels()){
      double QxyCorr = n->GetHits()[0]->GetPE()/att(n->GetZ());
      double QxzCorr = n->GetHits()[1]->GetPE()/att(n->GetY());
      double QyzCorr = n->GetHits()[2]->GetPE()/att(n->GetX());
      double QavCorr = (QxyCorr+QxzCorr+QyzCorr)/3;
      double Chi2 = (pow(QxyCorr-QavCorr, 2)+pow(QxzCorr-QavCorr, 2)+pow(QyzCorr-QavCorr, 2))/QavCorr;

        if(opt ==1){
            std::cout << iev << ", "
                      << n->GetX() << ", " << n->GetY() << ", " << n->GetZ() << ", "
                      << n->GetHits()[0]->GetPE() << ", " << n->GetHits()[1]->GetPE() << ", " << n->GetHits()[2]->GetPE() << ", "
                      << n->GetHits()[0]->GetMultiplicity() << ", " << n->GetHits()[1]->GetMultiplicity() << ", " << n->GetHits()[2]->GetMultiplicity() << ", "
                      << QxyCorr << ", " << QxzCorr << ", " << QyzCorr << ", "
                      << Chi2 << ", " << (QxyCorr-QxzCorr)/(QxyCorr+QxzCorr) << ", " << (QxyCorr-QyzCorr)/(QxyCorr+QyzCorr) << ", " << (QxzCorr-QyzCorr)/(QxzCorr+QyzCorr) << ", "
                      << n->GetTrueType() << "\n";
        }
        outCSVfile    << iev << ", "
                      << n->GetX() << ", " << n->GetY() << ", " << n->GetZ() << ", "
                      << n->GetHits()[0]->GetPE() << ", " << n->GetHits()[1]->GetPE() << ", " << n->GetHits()[2]->GetPE() << ", "
                      << n->GetHits()[0]->GetMultiplicity() << ", " << n->GetHits()[1]->GetMultiplicity() << ", " << n->GetHits()[2]->GetMultiplicity() << ", "
                      << QxyCorr << ", " << QxzCorr << ", " << QyzCorr << ", "
                      << Chi2 << ", " << (QxyCorr-QxzCorr)/(QxyCorr+QxzCorr) << ", " << (QxyCorr-QyzCorr)/(QxyCorr+QyzCorr) << ", " << (QxzCorr-QyzCorr)/(QxzCorr+QyzCorr) << ", "
                      << n->GetTrueType() << "\n";


    }
}

bool compareDistance(ND280SFGDVoxel* v1, ND280SFGDVoxel* v2) 
{ 
    return (v1->GetDistance() < v2->GetDistance()); 
}
