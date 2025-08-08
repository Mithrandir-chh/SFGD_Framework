{
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetTitleOffset(0.9,"y");

  gStyle->SetStripDecimals(kFALSE);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(99);

  gStyle->SetHistLineWidth(2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFuncWidth(2);

  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0000);

  // no hole x shift -0.8 y shift +0.35
  // one hole x shift 0  y shift -0.3
  // x0 y0 = no hole
  double radiusCut0 = 0.45;
  double radiusCut1 = 0.55;
  double radiusCut2 = 0.65;
  double radiusCut3 = 0.75;
  double shiftX0 = 0.5;
  double shiftX = 0.5;
  double shiftY0 = 0.5;
  double shiftY = 0.5;
  double rot = 1.;
  double threshold0 = 0.5;
  double elwAmount = 0.0000000001;
  double elwAmount2 = 1.5;

  double range1low  = 0;
  double range1high = 50;
  double range2low  = 100;
  double range2high = 150;
  double range3low  = 300;
  double range3high = 350;

  double d = 90;
  double e = 1;
  double m = 939;
  double c = 0.3;

  double gt = d/c;
  double gamm  = e/m + 1.;
  double v = sqrt(-( (1/gamm)*(1/gamm) * c * c - c * c));
  double tof = d / v;
  cout<<tof<<endl;

  double binEdge[10000];
  binEdge[0] = 1.;
  int nbin =0;
  for(int i=1;i<10000;i++){
    if (e>=300) break;
    tof=tof - 2.5;
    v = d/tof ;
    gamm = 1/ (sqrt((c*c-v*v)/(c*c)));
    e = (gamm -1)*m;
    binEdge[i] = e;
    cout<<tof<<" "<<v<<" "<<gamm<<" "<<binEdge[i]<<endl;
    nbin ++;
  }

  //TH1F* hist = new TH1F("","",nbin,binEdge);

  int nhour = 8;
  int nday = 5;
  double dataAmount = 5e7 * (3/2.) * nhour * nday;
  double perMacro = 346/2.;

  double perHour = (346/2.)* 100 * 3600 * 30;

TGraph* bFlux = new TGraph(41);
bFlux->SetPoint(1, 0 ,  1.26E-001);
//bFlux->SetPoint(1, 1.13E-001 ,	1.26E-001);
bFlux->SetPoint(2, 1.42E-001 ,	1.35E-001);
bFlux->SetPoint(3, 1.79E-001 ,	1.35E-001);
bFlux->SetPoint(4, 2.25E-001 ,	1.40E-001);
bFlux->SetPoint(5, 2.84E-001 ,	1.44E-001);
bFlux->SetPoint(6, 3.57E-001 ,	1.42E-001);
bFlux->SetPoint(7, 4.50E-001 ,	1.41E-001);
bFlux->SetPoint(8, 5.66E-001 ,	1.42E-001);
bFlux->SetPoint(9, 7.13E-001 ,	1.30E-001);
bFlux->SetPoint(10,8.97E-001 ,	1.17E-001);
bFlux->SetPoint(11,1.13E+000 ,	9.66E-002);
bFlux->SetPoint(12,1.42E+000 ,	7.64E-002);
bFlux->SetPoint(13,1.79E+000 ,	6.59E-002);
bFlux->SetPoint(14,2.25E+000 ,	5.34E-002);
bFlux->SetPoint(15,2.84E+000 ,	4.27E-002);
bFlux->SetPoint(16,3.57E+000 ,	3.36E-002);
bFlux->SetPoint(17,4.50E+000 ,	2.73E-002);
bFlux->SetPoint(18,5.66E+000 ,	1.97E-002);
bFlux->SetPoint(19,7.13E+000 ,	1.39E-002);
bFlux->SetPoint(20,8.97E+000 ,	9.51E-003);
bFlux->SetPoint(21,1.13E+001 ,	6.26E-003);
bFlux->SetPoint(22,1.42E+001 ,	4.00E-003);
bFlux->SetPoint(23,1.79E+001 ,	2.61E-003);
bFlux->SetPoint(24,2.26E+001 ,	2.13E-003);
bFlux->SetPoint(25,2.84E+001 ,	1.59E-003);
bFlux->SetPoint(26,3.57E+001 ,	1.30E-003);
bFlux->SetPoint(27,4.50E+001 ,	1.05E-003);
bFlux->SetPoint(28,5.66E+001 ,	8.59E-004);
bFlux->SetPoint(29,7.13E+001 ,	7.47E-004);
bFlux->SetPoint(30,8.97E+001 ,  5.93E-004);
bFlux->SetPoint(31,1.13E+002 ,  5.53E-004);
bFlux->SetPoint(32,1.42E+002 ,  4.69E-004);
bFlux->SetPoint(33,1.79E+002 ,  3.93E-004);
bFlux->SetPoint(34,2.25E+002 ,  3.24E-004);
bFlux->SetPoint(35,2.84E+002 ,  2.67E-004);
bFlux->SetPoint(36,3.57E+002 ,  1.88E-004);
bFlux->SetPoint(37,4.50E+002 ,  9.76E-005);
bFlux->SetPoint(38,5.38E+002 ,  9.08E-005);
bFlux->SetPoint(39,6.13E+002 ,  1.38E-004);
bFlux->SetPoint(40,6.88E+002 ,  9.64E-005);
bFlux->SetPoint(41,7.63E+002 ,  1.57E-005);

TCanvas* c0 = new TCanvas();
bFlux->SetTitle("Flux at 90 meter location");
bFlux->GetYaxis()->SetTitle("per proton per sr per MeV");
bFlux->GetXaxis()->SetTitle("neutron energy (MeV)");
c0->SetLogy();
c0->SetLogx();
c0->SetGridx();
c0->SetGridy();
bFlux->SetLineWidth(4);
bFlux->Draw();

TChain t("tree");
t.Add("build/neutron/*3*chunk*.root");

//TFile f("out.root");
//TTree* t = (TTree*)f.Get("tree");

Float_t neutronHitX[1000]; Float_t neutronHitY[1000]; Float_t neutronHitZ[1000];
Float_t neutronHitE[1000]; Float_t neutronHitPDG[1000]; Float_t initialE[1000];
Float_t vtx[3];
Float_t neutronHitT[1000];
Int_t delay[1000];
Float_t neutron_time[1000]; Float_t gamma_time[1000];
Float_t t0[1000]; Float_t spill[1000];
Float_t neutronParentId[1000];
Float_t neutronParentPDG[1000];

TH1D* eff1 = new TH1D("","",nbin,binEdge);
TH1D* eff2 = new TH1D("","",nbin,binEdge);

TH1D* stat1 = new TH1D("","",nbin,binEdge);
TH1D* stat2 = new TH1D("","",nbin,binEdge);
TH1D* stat3 = new TH1D("","",nbin,binEdge);
TH1D* stat4 = new TH1D("","",nbin,binEdge);

TH2D* stat21 = new TH2D("","",20,0,800,100,-50,50);
TH2D* stat22 = new TH2D("","",20,0,800,100,-50,50);
TH2D* stat23 = new TH2D("","",20,0,800,100,-50,50);
TH2D* stat24 = new TH2D("","",20,0,800,100,-50,50);

TH2D* xy = new TH2D("","",24,-12,12,8,-4,4);
TH2D* xz = new TH2D("","",24,-12,12,100,-50,50);
TH2D* yz = new TH2D("","",100,-50,50,8,-4,4);

TH2D* xy2 = new TH2D("","",24,-12,12,8,-4,4);
TH2D* xz2 = new TH2D("","",24,-12,12,100,-50,50);
TH2D* yz2 = new TH2D("","",100,-50,50,8,-4,4);

TH1D* protonPenDenom = new TH1D("","",20,0,800);
TH1D* protonPen = new TH1D("","Interaction at upstream half of the detector; Neutron Energy (MeV); Fraction of penetrating protons",20,0,800);
TH1D* pionPenDenom = new TH1D("","",20,0,800);
TH1D* pionPen = new TH1D("","Interaction at upstream half of the detector; Neutron Energy (MeV); Fraction of penetrating pions",20,0,800);

TH2D* xy3[3000];
TH2D* xz3[3000];
TH2D* yz3[3000];
TH2D* xy4[3000];
TH2D* xz4[3000];
TH2D* yz4[3000];
float elWeight[100]={};

TH1D* expHandler[100];
TH1D* expHandler2[100];
TH1D* expHandler_fluxW[100];
TH1D* expHandler_fluxWsctW[100];
TH1D* expHandler2_fluxW[100];
TH1D* expHandler2_fluxWsctW[100];
TH1D* expHandler_fluxWsctW2[100];
TH1D* expHandler2_fluxWsctW2[100];
TH1D* expHandler_fluxWsctW3[100];
TH1D* expHandler2_fluxWsctW3[100];

TH1D* timeSpan[100];
TH1D* timeSpanS[100];
TH2D* timeSpan2D[100];
TH2D* pdg2D[100];
TH1D* lay[100];
TH1D* h_nelastic[100];
for (int i=0;i<100;i++){
  lay[i] = new TH1D("","",100,-48,0);
  expHandler[i] = new TH1D("","min loc.; Z (cm, along beam) ;counts",24,-12,12);
  expHandler2[i] = new TH1D("","0th loc.; Z (cm, along beam) ;counts",24,-12,12);
  expHandler_fluxW[i] = new TH1D("","min loc. flux weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler_fluxWsctW[i] = new TH1D("","min loc. flux and sct weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler2_fluxW[i] = new TH1D("","0th loc. flux weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler2_fluxWsctW[i] = new TH1D("","0th loc. flux and sct weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler_fluxWsctW2[i] = new TH1D("","min loc. flux weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler_fluxWsctW3[i] = new TH1D("","min loc. flux and sct weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler2_fluxWsctW2[i] = new TH1D("","0th loc. flux weighted; Z (cm, along beam) ;counts",24,-12,12);
  expHandler2_fluxWsctW3[i] = new TH1D("","0th loc. flux and sct weighted; Z (cm, along beam) ;counts",24,-12,12);

  timeSpan[i] = new TH1D("","; Time span (ns) ;counts",100,-100,5000);
  timeSpanS[i] = new TH1D("","; Time span (ns) ;counts",100,-10,100);
  timeSpan2D[i] = new TH2D("","; Time span (ns) ;Neutron energy",100,-100,5000,100,0,800);
  pdg2D[i] = new TH2D("","",100,0,100,100,0,100);
  h_nelastic[i] = new TH1D("","",40,-20,20);
}
for (int i=0;i<3000;i++){
  xy3[i] = new TH2D("","",24,-12,12,8,-4,4);
  xz3[i] = new TH2D("","",32,-12,20,100,-50,50);
  yz3[i] = new TH2D("","",100,-50,50,8,-4,4);  
  xy4[i] = new TH2D("","",24,-12,12,8,-4,4);
  xz4[i] = new TH2D("","",32,-12,20,100,-50,50);
  yz4[i] = new TH2D("","",100,-50,50,8,-4,4);
}

TH2D* h2_test1 = new TH2D("","testing with x dist.; 0th; time",24,-12,12,24,-12,12);
TH2D* h_delay = new TH2D("","",100,0,100,20,0,20);
TH1D* h1_delay = new TH1D("","",20,0,20);
TH1D* h_initialE = new TH1D("","",nbin,binEdge);
TH1D* h_initialE2 = new TH1D("","",100,0,800);
TH1D* h_ts = new TH1D("","",625000/2.5,0,625000);

TH1D* h_layer0[3];
TH1D* h_layer1[3];
TH1D* h_layer2[3];
TH1D* h_layer3[3];
for(int ii=0;ii<3;ii++){
  h_layer0[ii] = new TH1D("","layer first hit",48,-24,24);
  h_layer1[ii] = new TH1D("","",48,-24,24);
  h_layer2[ii] = new TH1D("","",48,-24,24);
  h_layer3[ii] = new TH1D("","",48,-24,24);
}
TH1D* h_layerr0[3];
TH1D* h_layerr1[3];
TH1D* h_layerr2[3];
TH1D* h_layerr3[3];
for(int ii=0;ii<3;ii++){
  h_layerr0[ii] = new TH1D("","layer all hits",48,-24,24);
  h_layerr1[ii] = new TH1D("","",48,-24,24);
  h_layerr2[ii] = new TH1D("","",48,-24,24);
  h_layerr3[ii] = new TH1D("","",48,-24,24);
}

t.SetBranchAddress("neutronHitX",&neutronHitX);
t.SetBranchAddress("neutronHitY",&neutronHitY);
t.SetBranchAddress("neutronHitZ",&neutronHitZ);
t.SetBranchAddress("neutronHitE",&neutronHitE);
t.SetBranchAddress("neutronHitT",&neutronHitT);
t.SetBranchAddress("neutronHitPDG",&neutronHitPDG);
t.SetBranchAddress("initialE",&initialE);
t.SetBranchAddress("vtx",&vtx);
t.SetBranchAddress("neutronParentId",&neutronParentId);
t.SetBranchAddress("neutronParentPDG",&neutronParentPDG);
t.SetBranchAddress("delay",&delay);
t.SetBranchAddress("t0",&t0);
t.SetBranchAddress("gamma_time",&gamma_time);
t.SetBranchAddress("neutron_time",&neutron_time);
t.SetBranchAddress("spill",&spill);
t.SetBranchAddress("nElastic",&elWeight);

bool nTag=false;

cout<<t.GetEntries()<<endl;

int Nentry = t.GetEntries();
double tota = 0;
double coll1[10]={};	
double coll2[100]={};
int counter1 = 0;
int counter2 = 0;
int counter3 = 0;

for(Int_t i=0;i< Nentry; i++){
  t.GetEntry(i);

  nTag = false;

  if(i%1000==0)
    cout<<i<<endl;

  if( sqrt((vtx[0]-shiftX0)* (vtx[0]-shiftX0) + (vtx[1]-shiftY0)* (vtx[1]-shiftY0))>radiusCut0) continue;
  //if(sqrt((neutronHitX[0]*neutronHitX[0])+(neutronHitY[0]*neutronHitY[0]))>4 || neutronHitZ[0]> -24 ) continue;

  double wei = abs(bFlux->Eval(initialE[0]));
  //double wei = 1;

  double totalE=0;

  bool useIt = true;

  for(int j=0;j<100;j++){
    totalE += neutronHitE[j];
  }

  bool goingThru = false;
  bool goingThru2 = false;
  bool filled1 = false;
  bool filled2 = false;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//
  for(int iii=0;iii<99;iii++){
    h_nelastic[iii]->Fill(elWeight[iii]-1);
  }
  int minLoc10 = 0;
  double minHitT10 = 1e7;
  int minLoc11 = 0;
  double minHitT11 = 1e7;  
  bool filled10 = false;
  bool filled11 = false;  

  int maxLoc10 = 0;
  double maxHitT10 = -1e7;
  int maxLoc11 = 0;
  double maxHitT11 = -1e7;
  bool filled12 = false;
  bool filled13 = false;

  int minLoc12 = 0;
  double minHitT12 = 1e7;
  int maxLoc12 = 0;
  double maxHitT12 = -1e7;
  bool filled14 = false;
  bool filled15 = false;  

  for(int j=0;j<1000;j++){
    if( neutronHitE[j] > 0. && neutronHitT[j]< minHitT10){
      minLoc10 = j;
      minHitT10 = neutronHitT[j];
      filled10 = true;
    }
    if( neutronHitE[j] > 0.5 && neutronHitT[j]< minHitT11 ){
      minLoc11 = j;
      minHitT11 = neutronHitT[j];
      filled11 = true;
    }
    if( neutronHitE[j] > 5 && neutronHitT[j]< minHitT12 ){
      minLoc12 = j;
      minHitT12 = neutronHitT[j];
      filled14 = true;
    }

  if(neutronHitE[j] > 0.){
  if ( sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2))<radiusCut0 && neutronHitZ[j] !=0  ){
    h_layerr0[0] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut0 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut1 && neutronHitZ[j] !=0 ){
    h_layerr1[0] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut1 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut2 && neutronHitZ[j] !=0 ) {
    h_layerr2[0] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut2 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut3 && neutronHitZ[j] !=0 ){
    h_layerr3[0] -> Fill(neutronHitZ[j]);
  }
  }
  if(neutronHitE[j] > 0.5){
  if ( sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2))<radiusCut0 && neutronHitZ[j] !=0 ){
    h_layerr0[1] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut0 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut1 && neutronHitZ[j] !=0 ){
    h_layerr1[1] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut1 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut2 && neutronHitZ[j] !=0 ) {
    h_layerr2[1] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut2 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut3 && neutronHitZ[j] !=0 ){
    h_layerr3[1] -> Fill(neutronHitZ[j]);
  }
  }
  if(neutronHitE[j] > 5){
  if ( sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2))<radiusCut0 && neutronHitZ[j] !=0 ){
    h_layerr0[2] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut0 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut1 && neutronHitZ[j] !=0 ){
    h_layerr1[2] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut1 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut2 && neutronHitZ[j] !=0 ) {
    h_layerr2[2] -> Fill(neutronHitZ[j]);
  }
  else if (sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) > radiusCut2 && sqrt( pow(neutronHitX[j] - shiftX0,2)+pow(neutronHitY[j] - shiftY0,2)) < radiusCut3 && neutronHitZ[j] !=0 ){
    h_layerr3[2] -> Fill(neutronHitZ[j]);
  }
  }

  }

  h2_test1 -> Fill(neutronHitX[minLoc10],neutronHitX[0]);

  int minLoc = 0;
  double minHitT = 1e7;
  int minLoc2 = 0;
  double minHitT2 = 1e7;
  bool filled3 = false;
  bool filled4 = false;

  if(sqrt((vtx[0]-shiftX0)* (vtx[0]-shiftX0) + (vtx[1]-shiftY0)* (vtx[1]-shiftY0))<radiusCut0 && elWeight[9]>0){
    expHandler[int(neutronHitZ[minLoc11]+24)] -> Fill(neutronHitX[minLoc11]);
    expHandler2[int(neutronHitZ[0]+24)] -> Fill(neutronHitX[0]);

    expHandler_fluxW[int(neutronHitZ[minLoc11]+24)] -> Fill(neutronHitX[minLoc11], wei);
    expHandler_fluxWsctW[int(neutronHitZ[minLoc11]+24)] -> Fill(neutronHitX[minLoc11],  wei*TMath::Power(elwAmount, elWeight[9]-1));

    expHandler2_fluxW[int(neutronHitZ[0]+24)] -> Fill(neutronHitX[0], wei);
    expHandler2_fluxWsctW[int(neutronHitZ[0]+24)] -> Fill(neutronHitX[0],  wei*TMath::Power(elwAmount, elWeight[9]-1));

    if(elWeight[9] == 1){
      expHandler_fluxWsctW2[int(neutronHitZ[minLoc11]+24)] -> Fill(neutronHitX[minLoc11], 1); //wei*TMath::Power(elwAmount, elWeight[9]-1));	    
      expHandler2_fluxWsctW2[int(neutronHitZ[0]+24)] -> Fill(neutronHitX[0],  1);//wei*TMath::Power(elwAmount, elWeight[9]-1));
    }
    expHandler_fluxWsctW3[int(neutronHitZ[minLoc11]+24)] -> Fill(neutronHitX[minLoc11],  wei*TMath::Power(elwAmount2, elWeight[9]-1)); 
    expHandler2_fluxWsctW3[int(neutronHitZ[0]+24)] -> Fill(neutronHitX[0],  wei*TMath::Power(elwAmount2, elWeight[9]-1));
  }

  int locc1=0;
  int locc2=0;
  int locc3=0;

  if ( sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2))<radiusCut0 && neutronHitZ[locc1] !=0  ){
    h_layer0[0] -> Fill(neutronHitZ[locc1]);
  }
  else if (sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2)) > radiusCut0 && sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2)) < radiusCut1 && neutronHitZ[locc1] !=0 ){
    h_layer1[0] -> Fill(neutronHitZ[locc1]);
  }
  else if (sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2)) > radiusCut1 && sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2)) < radiusCut2 && neutronHitZ[locc1] !=0 ) {
    h_layer2[0] -> Fill(neutronHitZ[locc1]);
  } 
  else if (sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2)) > radiusCut2 && sqrt( pow(neutronHitX[locc1] - shiftX0,2)+pow(neutronHitY[locc1] - shiftY0,2)) < radiusCut3 && neutronHitZ[locc1] !=0 ){
    h_layer3[0] -> Fill(neutronHitZ[locc1]);
  }

  if ( sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2))<radiusCut0 && neutronHitZ[locc2] !=0 ){
    h_layer0[1] -> Fill(neutronHitZ[locc2]);
  }
  else if (sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2)) > radiusCut0 && sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2)) < radiusCut1 && neutronHitZ[locc2] !=0 ){
    h_layer1[1] -> Fill(neutronHitZ[locc2]);
  }
  else if (sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2)) > radiusCut1 && sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2)) < radiusCut2 && neutronHitZ[locc2] !=0 ) {
    h_layer2[1] -> Fill(neutronHitZ[locc2]);
  } 
  else if (sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2)) > radiusCut2 && sqrt( pow(neutronHitX[locc2] - shiftX0,2)+pow(neutronHitY[locc2] - shiftY0,2)) < radiusCut3 && neutronHitZ[locc2] !=0 ){
    h_layer3[1] -> Fill(neutronHitZ[locc2]);
  }

  if ( sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2))<radiusCut0 && neutronHitZ[locc3] !=0 ){
    h_layer0[2] -> Fill(neutronHitZ[locc3]);
  }
  else if (sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2)) > radiusCut0 && sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2)) < radiusCut1 && neutronHitZ[locc3] !=0 ){
    h_layer1[2] -> Fill(neutronHitZ[locc3]);
  }
  else if (sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2)) > radiusCut1 && sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2)) < radiusCut2 && neutronHitZ[locc3] !=0 ) {
    h_layer2[2] -> Fill(neutronHitZ[locc3]);
  } 
  else if (sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2)) > radiusCut2 && sqrt( pow(neutronHitX[locc3] - shiftX0,2)+pow(neutronHitY[locc3] - shiftY0,2)) < radiusCut3 && neutronHitZ[locc3] !=0 ){
    h_layer3[2] -> Fill(neutronHitZ[locc3]);
  }
  
}
for(int iii=0;iii<3;iii++){
  for(int ii=0;ii<h_layer0[iii]->GetNbinsX();ii++){
    h_layer1[iii]->SetBinContent(ii+1,h_layer1[iii]->GetBinContent(ii+1)/h_layer0[iii]->GetBinContent(ii+1));
    h_layer2[iii]->SetBinContent(ii+1,h_layer2[iii]->GetBinContent(ii+1)/h_layer0[iii]->GetBinContent(ii+1));
    h_layer3[iii]->SetBinContent(ii+1,h_layer3[iii]->GetBinContent(ii+1)/h_layer0[iii]->GetBinContent(ii+1));
    h_layerr1[iii]->SetBinContent(ii+1,h_layerr1[iii]->GetBinContent(ii+1)/h_layerr0[iii]->GetBinContent(ii+1));
    h_layerr2[iii]->SetBinContent(ii+1,h_layerr2[iii]->GetBinContent(ii+1)/h_layerr0[iii]->GetBinContent(ii+1));
    h_layerr3[iii]->SetBinContent(ii+1,h_layerr3[iii]->GetBinContent(ii+1)/h_layerr0[iii]->GetBinContent(ii+1));
  }
}

h_nelastic[0]->Scale(1./h_nelastic[0]->Integral(21,40));
h_nelastic[9]->Scale(1./h_nelastic[9]->Integral(21,40));
h_nelastic[98]->Scale(1./h_nelastic[98]->Integral(21,40));

new TCanvas();
h_layer1[0]->SetLineColor(1);
h_layer1[0]->SetLineWidth(3);
h_layer2[0]->SetLineColor(2);
h_layer2[0]->SetLineWidth(3);
h_layer3[0]->SetLineColor(4);
h_layer3[0]->SetLineWidth(3);
h_layer1[1]->SetLineColor(1);
h_layer1[1]->SetLineWidth(3);
h_layer1[1]->SetLineStyle(2);
h_layer2[1]->SetLineColor(2);
h_layer2[1]->SetLineWidth(3);
h_layer2[1]->SetLineStyle(2);
h_layer3[1]->SetLineColor(4);
h_layer3[1]->SetLineWidth(3);
h_layer3[1]->SetLineStyle(2);
h_layer1[2]->SetLineColor(1);
h_layer1[2]->SetLineWidth(3);
h_layer1[2]->SetLineStyle(3);
h_layer2[2]->SetLineColor(2);
h_layer2[2]->SetLineWidth(3);
h_layer2[2]->SetLineStyle(3);
h_layer3[2]->SetLineColor(4);
h_layer3[2]->SetLineWidth(3);
h_layer3[2]->SetLineStyle(3);
h_layer1[0]->Draw();
h_layer2[0]->Draw("same");
h_layer3[0]->Draw("same");
h_layer1[1]->Draw("same");
h_layer2[1]->Draw("same");
h_layer3[1]->Draw("same");
h_layer1[2]->Draw("same");
h_layer2[2]->Draw("same");
h_layer3[2]->Draw("same");
TLegend *lgg_layer = new TLegend(0.7, 0.7, 0.9, 0.9);
lgg_layer->AddEntry(h_layer1[0], "0 MeV threshold inner ratio", "l");
lgg_layer->AddEntry(h_layer2[0], "0 MeV threshold middle ratio", "l");
lgg_layer->AddEntry(h_layer3[0], "0 MeV threshold outer ratio", "l");
lgg_layer->AddEntry(h_layer1[1], "0.5 MeV threshold inner ratio", "l");
lgg_layer->AddEntry(h_layer2[1], "0.5 MeV threshold middle ratio", "l");
lgg_layer->AddEntry(h_layer3[1], "0.5 MeV threshold outer ratio", "l");
lgg_layer->AddEntry(h_layer1[2], "5 MeV threshold inner ratio", "l");
lgg_layer->AddEntry(h_layer2[2], "5 MeV threshold middle ratio", "l");
lgg_layer->AddEntry(h_layer3[2], "5 MeV threshold outer ratio", "l");
lgg_layer->Draw();

new TCanvas();
h_layerr1[0]->SetLineColor(1);
h_layerr1[0]->SetLineWidth(3);
h_layerr2[0]->SetLineColor(2);
h_layerr2[0]->SetLineWidth(3);
h_layerr3[0]->SetLineColor(4);
h_layerr3[0]->SetLineWidth(3);
h_layerr1[1]->SetLineColor(1);
h_layerr1[1]->SetLineWidth(3);
h_layerr1[1]->SetLineStyle(2);
h_layerr2[1]->SetLineColor(2);
h_layerr2[1]->SetLineWidth(3);
h_layerr2[1]->SetLineStyle(2);
h_layerr3[1]->SetLineColor(4);
h_layerr3[1]->SetLineWidth(3);
h_layerr3[1]->SetLineStyle(2);
h_layerr1[2]->SetLineColor(1);
h_layerr1[2]->SetLineWidth(3);
h_layerr1[2]->SetLineStyle(3);
h_layerr2[2]->SetLineColor(2);
h_layerr2[2]->SetLineWidth(3);
h_layerr2[2]->SetLineStyle(3);
h_layerr3[2]->SetLineColor(4);
h_layerr3[2]->SetLineWidth(3);
h_layerr3[2]->SetLineStyle(3);
h_layerr1[0]->Draw();
h_layerr2[0]->Draw("same");
h_layerr3[0]->Draw("same");
h_layerr1[1]->Draw("same");
h_layerr2[1]->Draw("same");
h_layerr3[1]->Draw("same");
h_layerr1[2]->Draw("same");
h_layerr2[2]->Draw("same");
h_layerr3[2]->Draw("same");
TLegend *lgg_layerr = new TLegend(0.7, 0.7, 0.9, 0.9);
lgg_layerr->AddEntry(h_layerr1[0], "0 MeV threshold inner ratio", "l");
lgg_layerr->AddEntry(h_layerr2[0], "0 MeV threshold middle ratio", "l");
lgg_layerr->AddEntry(h_layerr3[0], "0 MeV threshold outer ratio", "l");
lgg_layerr->AddEntry(h_layerr1[1], "0.5 MeV threshold inner ratio", "l");
lgg_layerr->AddEntry(h_layerr2[1], "0.5 MeV threshold middle ratio", "l");
lgg_layerr->AddEntry(h_layerr3[1], "0.5 MeV threshold outer ratio", "l");
lgg_layerr->AddEntry(h_layerr1[2], "5 MeV threshold inner ratio", "l");
lgg_layerr->AddEntry(h_layerr2[2], "5 MeV threshold middle ratio", "l");
lgg_layerr->AddEntry(h_layerr3[2], "5 MeV threshold outer ratio", "l");
lgg_layerr->Draw();

new TCanvas();
h_nelastic[0]->SetLineColor(1);
h_nelastic[9]->SetLineColor(2);
h_nelastic[98]->SetLineColor(4);
h_nelastic[0]->SetLineWidth(3);
h_nelastic[9]->SetLineWidth(3);
h_nelastic[98]->SetLineWidth(3);
h_nelastic[0]->Draw("hist");
h_nelastic[9]->Draw("same");
h_nelastic[98]->Draw("same");
TLegend *lgg_nelastic = new TLegend(0.7, 0.7, 0.9, 0.9);
lgg_nelastic->AddEntry(h_nelastic[0], "0 MeV threshold", "l");
lgg_nelastic->AddEntry(h_nelastic[9], "0.5 MeV threshold", "l");
lgg_nelastic->AddEntry(h_nelastic[98], "5 MeV threshold", "l");
lgg_nelastic->Draw();

new TCanvas();
h2_test1->Draw("colz");

TCanvas* c1 = new TCanvas();
c1->Divide(6,6);
for(int i=0;i<36;i++){
  c1->cd(i+1);
  expHandler[i]->Draw();
}

TCanvas* c2 = new TCanvas();
c2->Divide(6,6);
for(int i=0;i<36;i++){
  c2->cd(i+1);
  expHandler2[i]->Draw();
}

TGraph* gLeft = new TGraph(48);
TGraph* gRight = new TGraph(48);
TGraph* gLeft2 = new TGraph(48);
TGraph* gRight2 = new TGraph(48);
TGraph* gLR   = new TGraph(48);
TGraph* gTB   = new TGraph(48);

TGraph* gLeft_fluxW = new TGraph(48);
TGraph* gRight_fluxW = new TGraph(48);
TGraph* gLeft2_fluxW = new TGraph(48);
TGraph* gRight2_fluxW = new TGraph(48);
TGraph* gLR_fluxW   = new TGraph(48);
TGraph* gTB_fluxW   = new TGraph(48);

TGraph* gLeft_fluxWsctW = new TGraph(48);
TGraph* gRight_fluxWsctW = new TGraph(48);
TGraph* gLeft2_fluxWsctW = new TGraph(48);
TGraph* gRight2_fluxWsctW = new TGraph(48);
TGraph* gLR_fluxWsctW   = new TGraph(48);
TGraph* gTB_fluxWsctW   = new TGraph(48);

TGraph* gLeft_fluxWsctW2 = new TGraph(48);
TGraph* gRight_fluxWsctW2 = new TGraph(48);
TGraph* gLeft2_fluxWsctW2 = new TGraph(48);
TGraph* gRight2_fluxWsctW2 = new TGraph(48);
TGraph* gLR_fluxWsctW2   = new TGraph(48);
TGraph* gTB_fluxWsctW2   = new TGraph(48);

TGraph* gLeft_fluxWsctW3 = new TGraph(48);
TGraph* gRight_fluxWsctW3 = new TGraph(48);
TGraph* gLeft2_fluxWsctW3 = new TGraph(48);
TGraph* gRight2_fluxWsctW3 = new TGraph(48);
TGraph* gLR_fluxWsctW3   = new TGraph(48);
TGraph* gTB_fluxWsctW3   = new TGraph(48);

for(int i=0;i<48;i++){
  if(expHandler[i]->GetBinContent(13)>0){
    gLeft->SetPoint(i,i, expHandler[i]->GetBinContent(12)/expHandler[i]->GetBinContent(13));
    gRight->SetPoint(i,i, expHandler[i]->GetBinContent(14)/expHandler[i]->GetBinContent(13));
    TF1* fb = new TF1("fb","gaus(0)",-3,3);
    expHandler[i] -> Fit("fb","R");
    //gLR  ->SetPoint(i, i, fb->GetParameter(2));
    gLR ->SetPoint(i,i, expHandler[i]->GetRMS());
  }
  if(expHandler2[i]->GetBinContent(13)>0){
    gLeft2->SetPoint(i,i, expHandler2[i]->GetBinContent(12)/expHandler2[i]->GetBinContent(13));
    gRight2->SetPoint(i,i, expHandler2[i]->GetBinContent(14)/expHandler2[i]->GetBinContent(13));
    TF1* fb2 = new TF1("fb2","gaus(0)",-3,3);
    expHandler2[i] -> Fit("fb2","R");
    //gTB  ->SetPoint(i, i, fb2->GetParameter(2));
    gTB  ->SetPoint(i, i, expHandler2[i]->GetRMS());
  }

  //cout<<"testing : "<<expHandler_fluxW[i]->GetBinContent(12)<<" "<<expHandler_fluxW[i]->GetBinContent(13)<<endl;
  //cout<<expHandler_fluxWsctW[i]->GetBinContent(12)<<" "<<expHandler_fluxWsctW[i]->GetBinContent(13)<<endl;
  if(expHandler_fluxW[i]->GetBinContent(13)>0){
    gLeft_fluxW->SetPoint(i,i, expHandler_fluxW[i]->GetBinContent(12)/expHandler_fluxW[i]->GetBinContent(13));
    gRight_fluxW->SetPoint(i,i, expHandler_fluxW[i]->GetBinContent(14)/expHandler_fluxW[i]->GetBinContent(13));
    TF1* fb = new TF1("fb","gaus(0)",-3,3);
    expHandler_fluxW[i] -> Fit("fb","R");
    gLR_fluxW  ->SetPoint(i, i, fb->GetParameter(2));
  }
  if(expHandler2_fluxW[i]->GetBinContent(13)>0){
    gLeft2_fluxW->SetPoint(i,i, expHandler2_fluxW[i]->GetBinContent(12)/expHandler2_fluxW[i]->GetBinContent(13));
    gRight2_fluxW->SetPoint(i,i, expHandler2_fluxW[i]->GetBinContent(14)/expHandler2_fluxW[i]->GetBinContent(13));
    TF1* fb2 = new TF1("fb2","gaus(0)",-3,3);
    expHandler2_fluxW[i] -> Fit("fb2","R");
    gTB_fluxW  ->SetPoint(i, i, fb2->GetParameter(2));
  }

  if(expHandler_fluxWsctW[i]->GetBinContent(13)>0){
    gLeft_fluxWsctW->SetPoint(i,i, expHandler_fluxWsctW[i]->GetBinContent(12)/expHandler_fluxWsctW[i]->GetBinContent(13));
    gRight_fluxWsctW->SetPoint(i,i, expHandler_fluxWsctW[i]->GetBinContent(14)/expHandler_fluxWsctW[i]->GetBinContent(13));
    TF1* fb = new TF1("fb","gaus(0)",-3,3);
    expHandler_fluxWsctW[i] -> Fit("fb","R");
    gLR_fluxWsctW  ->SetPoint(i, i, fb->GetParameter(2));
  }
 if(expHandler2_fluxWsctW[i]->GetBinContent(13)>0){
    gLeft2_fluxWsctW->SetPoint(i,i, expHandler2_fluxWsctW[i]->GetBinContent(12)/expHandler2_fluxWsctW[i]->GetBinContent(13));
    gRight2_fluxWsctW->SetPoint(i,i, expHandler2_fluxWsctW[i]->GetBinContent(14)/expHandler2_fluxWsctW[i]->GetBinContent(13));
    TF1* fb2 = new TF1("fb2","gaus(0)",-3,3);
    expHandler2_fluxWsctW[i] -> Fit("fb2","R");
    gTB_fluxWsctW  ->SetPoint(i, i, fb2->GetParameter(2));
  }
  

  if(expHandler_fluxWsctW2[i]->GetBinContent(13)>0){
    gLeft_fluxWsctW2->SetPoint(i,i, expHandler_fluxWsctW2[i]->GetBinContent(12)/expHandler_fluxWsctW2[i]->GetBinContent(13));
    gRight_fluxWsctW2->SetPoint(i,i, expHandler_fluxWsctW2[i]->GetBinContent(14)/expHandler_fluxWsctW2[i]->GetBinContent(13));
    TF1* fb = new TF1("fb","gaus(0)",-3,3);
    expHandler_fluxWsctW2[i] -> Fit("fb","R");
    gLR_fluxWsctW2  ->SetPoint(i, i, fb->GetParameter(2));
  }
  if(expHandler2_fluxWsctW2[i]->GetBinContent(13)>0){
    gLeft2_fluxWsctW2->SetPoint(i,i, expHandler2_fluxWsctW2[i]->GetBinContent(12)/expHandler2_fluxWsctW2[i]->GetBinContent(13));
    gRight2_fluxWsctW2->SetPoint(i,i, expHandler2_fluxWsctW2[i]->GetBinContent(14)/expHandler2_fluxWsctW2[i]->GetBinContent(13));
    TF1* fb2 = new TF1("fb2","gaus(0)",-3,3);
    expHandler2_fluxWsctW2[i] -> Fit("fb2","R");
    gTB_fluxWsctW2  ->SetPoint(i, i, fb2->GetParameter(2));
  }

  if(expHandler_fluxWsctW3[i]->GetBinContent(13)>0){
    gLeft_fluxWsctW3->SetPoint(i,i, expHandler_fluxWsctW3[i]->GetBinContent(12)/expHandler_fluxWsctW3[i]->GetBinContent(13));
    gRight_fluxWsctW3->SetPoint(i,i, expHandler_fluxWsctW3[i]->GetBinContent(14)/expHandler_fluxWsctW3[i]->GetBinContent(13));
    TF1* fb = new TF1("fb","gaus(0)",-3,3);
    expHandler_fluxWsctW3[i] -> Fit("fb","R");
    gLR_fluxWsctW3  ->SetPoint(i, i, fb->GetParameter(2));
  }
  if(expHandler2_fluxWsctW3[i]->GetBinContent(13)>0){
    gLeft2_fluxWsctW3->SetPoint(i,i, expHandler2_fluxWsctW3[i]->GetBinContent(12)/expHandler2_fluxWsctW3[i]->GetBinContent(13));
    gRight2_fluxWsctW3->SetPoint(i,i, expHandler2_fluxWsctW3[i]->GetBinContent(14)/expHandler2_fluxWsctW3[i]->GetBinContent(13));
    TF1* fb2 = new TF1("fb2","gaus(0)",-3,3);
    expHandler2_fluxWsctW3[i] -> Fit("fb2","R");
    gTB_fluxWsctW3  ->SetPoint(i, i, fb2->GetParameter(2));
  }

}

  new TCanvas();
  gRight->SetMarkerStyle(20);
  gRight->SetMarkerSize(2);
  gLeft->SetMarkerStyle(21);
  gLeft->SetMarkerSize(2);
  gLeft->SetMarkerColor(4);
  TMultiGraph* mgg = new TMultiGraph();
  mgg->SetTitle("min loc. no weight.");
  mgg->Add(gLeft,"AP");
  mgg->Add(gRight,"AP");
  mgg->Draw("AP");
  mgg->GetXaxis()->SetTitle("layer");
  mgg->GetYaxis()->SetTitle("fraction of neighbour bins");
  TLegend *lgg = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg->AddEntry(gLeft, "cube on the left - min", "p");
  lgg->AddEntry(gRight, "cube on the right - min", "p");
  lgg->Draw();

  new TCanvas();
  gRight2->SetMarkerStyle(20);
  gRight2->SetMarkerSize(2);
  gLeft2->SetMarkerStyle(21);
  gLeft2->SetMarkerSize(2);
  gLeft2->SetMarkerColor(4);
  TMultiGraph* mgg2 = new TMultiGraph();
  mgg2->SetTitle("0th loc. no weight.");
  mgg2->Add(gLeft2,"AP");
  mgg2->Add(gRight2,"AP");
  mgg2->Draw("AP");
  mgg2->GetXaxis()->SetTitle("layer");
  mgg2->GetYaxis()->SetTitle("fraction of neighbour bins");
  lgg = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg->AddEntry(gLeft2, "cube on the left - 0th", "p");
  lgg->AddEntry(gRight2, "cube on the right - 0th", "p");
  lgg->Draw();

  new TCanvas();
  gTB->SetMarkerStyle(20);
  gTB->SetMarkerSize(2);
  gLR->SetMarkerStyle(21);
  gLR->SetMarkerSize(2);
  gLR->SetMarkerColor(4);
  TMultiGraph* mgg3 = new TMultiGraph();
  mgg3->SetTitle("min loc. no weight.");
  mgg3->Add(gTB,"AP");
  mgg3->Add(gLR,"AP");
  mgg3->Draw("AP");
  mgg3->GetXaxis()->SetTitle("layer");
  mgg3->GetYaxis()->SetTitle("Gaussian width");
  TLegend *lgg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg3->AddEntry(gTB, "0th ", "p");
  lgg3->AddEntry(gLR, "min", "p");
  lgg3->Draw();

  new TCanvas();
  gRight_fluxW->SetMarkerStyle(20);
  gRight_fluxW->SetMarkerSize(2);
  gLeft_fluxW->SetMarkerStyle(21);
  gLeft_fluxW->SetMarkerSize(2);
  gLeft_fluxW->SetMarkerColor(4);
  TMultiGraph* mgg_fluxW = new TMultiGraph();
  mgg_fluxW->SetTitle("min loc. flux weight");
  mgg_fluxW->Add(gLeft_fluxW,"AP");
  mgg_fluxW->Add(gRight_fluxW,"AP");
  mgg_fluxW->Draw("AP");
  mgg_fluxW->GetXaxis()->SetTitle("layer");
  mgg_fluxW->GetYaxis()->SetTitle("fraction of neighbour bins");
  TLegend *lgg_fluxW = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxW->AddEntry(gLeft_fluxW, "cube on the left - min", "p");
  lgg_fluxW->AddEntry(gRight_fluxW, "cube on the right - min", "p");
  lgg_fluxW->Draw();

  new TCanvas();
  gRight2_fluxW->SetMarkerStyle(20);
  gRight2_fluxW->SetMarkerSize(2);
  gLeft2_fluxW->SetMarkerStyle(21);
  gLeft2_fluxW->SetMarkerSize(2);
  gLeft2_fluxW->SetMarkerColor(4);
  TMultiGraph* mgg2_fluxW = new TMultiGraph();
  mgg2_fluxW->SetTitle("0th loc. flux weight");
  mgg2_fluxW->Add(gLeft2_fluxW,"AP");
  mgg2_fluxW->Add(gRight2_fluxW,"AP");
  mgg2_fluxW->Draw("AP");
  mgg2_fluxW->GetXaxis()->SetTitle("layer");
  mgg2_fluxW->GetYaxis()->SetTitle("fraction of neighbour bins");
  lgg_fluxW = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxW->AddEntry(gLeft2_fluxW, "cube on the left - 0th", "p");
  lgg_fluxW->AddEntry(gRight2_fluxW, "cube on the right - 0th", "p");
  lgg_fluxW->Draw();

  new TCanvas();
  gTB_fluxW->SetMarkerStyle(20);
  gTB_fluxW->SetMarkerSize(2);
  gLR_fluxW->SetMarkerStyle(21);
  gLR_fluxW->SetMarkerSize(2);
  gLR_fluxW->SetMarkerColor(4);
  TMultiGraph* mgg3_fluxW = new TMultiGraph();
  mgg3->SetTitle("min loc. flux weight");
  mgg3_fluxW->Add(gTB_fluxW,"AP");
  mgg3_fluxW->Add(gLR_fluxW,"AP");
  mgg3_fluxW->Draw("AP");
  mgg3_fluxW->GetXaxis()->SetTitle("layer");
  mgg3_fluxW->GetYaxis()->SetTitle("Gaussian width");
  TLegend *lgg3_fluxW = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg3_fluxW->AddEntry(gTB_fluxW, "0th ", "p");
  lgg3_fluxW->AddEntry(gLR_fluxW, "min", "p");
  lgg3_fluxW->Draw();  
  
  new TCanvas();
  gRight_fluxWsctW->SetMarkerStyle(20);
  gRight_fluxWsctW->SetMarkerSize(2);
  gLeft_fluxWsctW->SetMarkerStyle(21);
  gLeft_fluxWsctW->SetMarkerSize(2);
  gLeft_fluxWsctW->SetMarkerColor(4);
  TMultiGraph* mgg_fluxWsctW = new TMultiGraph();
  mgg_fluxWsctW->SetTitle("min loc. flux+sct weight");
  mgg_fluxWsctW->Add(gLeft_fluxWsctW,"AP");
  mgg_fluxWsctW->Add(gRight_fluxWsctW,"AP");
  mgg_fluxWsctW->Draw("AP");
  mgg_fluxWsctW->GetXaxis()->SetTitle("layer");
  mgg_fluxWsctW->GetYaxis()->SetTitle("fraction of neighbour bins");
  TLegend *lgg_fluxWsctW = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxWsctW->AddEntry(gLeft_fluxWsctW, "cube on the left - min", "p");
  lgg_fluxWsctW->AddEntry(gRight_fluxWsctW, "cube on the right - min", "p");
  lgg_fluxWsctW->Draw();

  new TCanvas();
  gRight2_fluxWsctW->SetMarkerStyle(20);
  gRight2_fluxWsctW->SetMarkerSize(2);
  gLeft2_fluxWsctW->SetMarkerStyle(21);
  gLeft2_fluxWsctW->SetMarkerSize(2);
  gLeft2_fluxWsctW->SetMarkerColor(4);
  TMultiGraph* mgg2_fluxWsctW = new TMultiGraph();
  mgg2_fluxWsctW->SetTitle("0th loc. flux+sct weight");
  mgg2_fluxWsctW->Add(gLeft2_fluxWsctW,"AP");
  mgg2_fluxWsctW->Add(gRight2_fluxWsctW,"AP");
  mgg2_fluxWsctW->Draw("AP");
  mgg2_fluxWsctW->GetXaxis()->SetTitle("layer");
  mgg2_fluxWsctW->GetYaxis()->SetTitle("fraction of neighbour bins");
  lgg_fluxWsctW = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxWsctW->AddEntry(gLeft2_fluxWsctW, "cube on the left - 0th", "p");
  lgg_fluxWsctW->AddEntry(gRight2_fluxWsctW, "cube on the right - 0th", "p");
  lgg_fluxWsctW->Draw();

  new TCanvas();
  gTB_fluxWsctW->SetMarkerStyle(20);
  gTB_fluxWsctW->SetMarkerSize(2);
  gLR_fluxWsctW->SetMarkerStyle(21);
  gLR_fluxWsctW->SetMarkerSize(2);
  gLR_fluxWsctW->SetMarkerColor(4);
  TMultiGraph* mgg3_fluxWsctW = new TMultiGraph();
  mgg3_fluxWsctW->SetTitle("min loc. flux+sct weight");
  mgg3_fluxWsctW->Add(gTB_fluxWsctW,"AP");
  mgg3_fluxWsctW->Add(gLR_fluxWsctW,"AP");
  mgg3_fluxWsctW->Draw("AP");
  mgg3_fluxWsctW->GetXaxis()->SetTitle("layer");
  mgg3_fluxWsctW->GetYaxis()->SetTitle("Gaussian width");
  TLegend *lgg3_fluxWsctW = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg3_fluxWsctW->AddEntry(gTB_fluxWsctW, "0th ", "p");
  lgg3_fluxWsctW->AddEntry(gLR_fluxWsctW, "min", "p");
  lgg3_fluxWsctW->Draw();


  new TCanvas();
  gRight_fluxWsctW2->SetMarkerStyle(20);
  gRight_fluxWsctW2->SetMarkerSize(2);
  gLeft_fluxWsctW2->SetMarkerStyle(21);
  gLeft_fluxWsctW2->SetMarkerSize(2);
  gLeft_fluxWsctW2->SetMarkerColor(4);
  TMultiGraph* mgg_fluxWsctW2 = new TMultiGraph();
  mgg_fluxWsctW2->SetTitle("min loc. flux weight + no sct");
  mgg_fluxWsctW2->Add(gLeft_fluxWsctW2,"AP");
  mgg_fluxWsctW2->Add(gRight_fluxWsctW2,"AP");
  mgg_fluxWsctW2->Draw("AP");
  mgg_fluxWsctW2->GetXaxis()->SetTitle("layer");
  mgg_fluxWsctW2->GetYaxis()->SetTitle("fraction of neighbour bins");
  TLegend *lgg_fluxWsctW2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxWsctW2->AddEntry(gLeft_fluxWsctW2, "cube on the left - min", "p");
  lgg_fluxWsctW2->AddEntry(gRight_fluxWsctW2, "cube on the right - min", "p");
  lgg_fluxWsctW2->Draw();

  new TCanvas();
  gRight2_fluxWsctW2->SetMarkerStyle(20);
  gRight2_fluxWsctW2->SetMarkerSize(2);
  gLeft2_fluxWsctW2->SetMarkerStyle(21);
  gLeft2_fluxWsctW2->SetMarkerSize(2);
  gLeft2_fluxWsctW2->SetMarkerColor(4);
  TMultiGraph* mgg2_fluxWsctW2 = new TMultiGraph();
  mgg2_fluxWsctW2->SetTitle("0th loc. flux weight + no sct");
  mgg2_fluxWsctW2->Add(gLeft2_fluxWsctW2,"AP");
  mgg2_fluxWsctW2->Add(gRight2_fluxWsctW2,"AP");
  mgg2_fluxWsctW2->Draw("AP");
  mgg2_fluxWsctW2->GetXaxis()->SetTitle("layer");
  mgg2_fluxWsctW2->GetYaxis()->SetTitle("fraction of neighbour bins");
  lgg_fluxWsctW2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxWsctW2->AddEntry(gLeft2_fluxWsctW2, "cube on the left - 0th", "p");
  lgg_fluxWsctW2->AddEntry(gRight2_fluxWsctW2, "cube on the right - 0th", "p");
  lgg_fluxWsctW2->Draw();

  new TCanvas();
  gTB_fluxWsctW2->SetMarkerStyle(20);
  gTB_fluxWsctW2->SetMarkerSize(2);
  gLR_fluxWsctW2->SetMarkerStyle(21);
  gLR_fluxWsctW2->SetMarkerSize(2);
  gLR_fluxWsctW2->SetMarkerColor(4);
  TMultiGraph* mgg3_fluxWsctW2 = new TMultiGraph();
  mgg3_fluxWsctW2->SetTitle("min loc. flux weight + no sct");
  mgg3_fluxWsctW2->Add(gTB_fluxWsctW2,"AP");
  mgg3_fluxWsctW2->Add(gLR_fluxWsctW2,"AP");
  mgg3_fluxWsctW2->Draw("AP");
  mgg3_fluxWsctW2->GetXaxis()->SetTitle("layer");
  mgg3_fluxWsctW2->GetYaxis()->SetTitle("Gaussian width");
  TLegend *lgg3_fluxWsctW2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg3_fluxWsctW2->AddEntry(gTB_fluxWsctW2, "0th ", "p");
  lgg3_fluxWsctW2->AddEntry(gLR_fluxWsctW2, "min", "p");
  lgg3_fluxWsctW2->Draw();

  new TCanvas();
  gRight_fluxWsctW3->SetMarkerStyle(20);
  gRight_fluxWsctW3->SetMarkerSize(2);
  gLeft_fluxWsctW3->SetMarkerStyle(21);
  gLeft_fluxWsctW3->SetMarkerSize(2);
  gLeft_fluxWsctW3->SetMarkerColor(4);
  TMultiGraph* mgg_fluxWsctW3 = new TMultiGraph();
  mgg_fluxWsctW3->SetTitle("min loc. flux+sct (larger) weight");
  mgg_fluxWsctW3->Add(gLeft_fluxWsctW3,"AP");
  mgg_fluxWsctW3->Add(gRight_fluxWsctW3,"AP");
  mgg_fluxWsctW3->Draw("AP");
  mgg_fluxWsctW3->GetXaxis()->SetTitle("layer");
  mgg_fluxWsctW3->GetYaxis()->SetTitle("fraction of neighbour bins");
  TLegend *lgg_fluxWsctW3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxWsctW3->AddEntry(gLeft_fluxWsctW3, "cube on the left - min", "p");
  lgg_fluxWsctW3->AddEntry(gRight_fluxWsctW3, "cube on the right - min", "p");
  lgg_fluxWsctW3->Draw();

  new TCanvas();
  gRight2_fluxWsctW3->SetMarkerStyle(20);
  gRight2_fluxWsctW3->SetMarkerSize(2);
  gLeft2_fluxWsctW3->SetMarkerStyle(21);
  gLeft2_fluxWsctW3->SetMarkerSize(2);
  gLeft2_fluxWsctW3->SetMarkerColor(4);
  TMultiGraph* mgg2_fluxWsctW3 = new TMultiGraph();
  mgg2_fluxWsctW3->SetTitle("0th loc. flux+sct (larger) weight");
  mgg2_fluxWsctW3->Add(gLeft2_fluxWsctW3,"AP");
  mgg2_fluxWsctW3->Add(gRight2_fluxWsctW3,"AP");
  mgg2_fluxWsctW3->Draw("AP");
  mgg2_fluxWsctW3->GetXaxis()->SetTitle("layer");
  mgg2_fluxWsctW3->GetYaxis()->SetTitle("fraction of neighbour bins");
  lgg_fluxWsctW3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg_fluxWsctW3->AddEntry(gLeft2_fluxWsctW3, "cube on the left - 0th", "p");
  lgg_fluxWsctW3->AddEntry(gRight2_fluxWsctW3, "cube on the right - 0th", "p");
  lgg_fluxWsctW3->Draw();

  new TCanvas();
  gTB_fluxWsctW3->SetMarkerStyle(20);
  gTB_fluxWsctW3->SetMarkerSize(2);
  gLR_fluxWsctW3->SetMarkerStyle(21);
  gLR_fluxWsctW3->SetMarkerSize(2);
  gLR_fluxWsctW3->SetMarkerColor(4);
  TMultiGraph* mgg3_fluxWsctW3 = new TMultiGraph();
  mgg3_fluxWsctW3->SetTitle("min loc. flux+sct (larger) weight");
  mgg3_fluxWsctW3->Add(gTB_fluxWsctW3,"AP");
  mgg3_fluxWsctW3->Add(gLR_fluxWsctW3,"AP");
  mgg3_fluxWsctW3->Draw("AP");
  mgg3_fluxWsctW3->GetXaxis()->SetTitle("layer");
  mgg3_fluxWsctW3->GetYaxis()->SetTitle("Gaussian width");
  TLegend *lgg3_fluxWsctW3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  lgg3_fluxWsctW3->AddEntry(gTB_fluxWsctW3, "0th ", "p");
  lgg3_fluxWsctW3->AddEntry(gLR_fluxWsctW3, "min", "p");
  lgg3_fluxWsctW3->Draw();
  

}

