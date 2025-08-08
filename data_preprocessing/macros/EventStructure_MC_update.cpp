
#include "../src/global_header.hh"
#include <cstring>
#include "TMath.h"
#include "TGraph.h"
#include "TRandom3.h"


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

  TChain t("tree");
  t.Add("input.root");
  //t.Add("/home/guang/work/sfgd_framework/data_preprocessing/build/FHC_200_24x8x48_wFiber_wDead_realSize_wTyvek_noGap_v2_dump_Jun16_bert2.root");

  const double cubeSize  = 1.01;
  const double cubeSizeY = 1.025;
  const double tickSize = 2.5;
  const double ifNeutronSample = false;
  const bool IsRotated = false;

  //constants for energy calibration
  const double CBIRKS = 0.00208; // mm/MeV
  //const double EdepToPhotConv_FGD = 70.8; // CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 

  const double EdepToPhotConv_FGD = 59.7; // from the cern beam test paper
  const double EdepToPhotConv_FGD_1_vert = 58.62; // from the cern beam test paper
  const double EdepToPhotConv_FGD_1_hori = 57.67; // from the cern beam test paper
  const double EdepToPhotConv_FGD_2 = EdepToPhotConv_FGD_1_vert * (51.56/52.53); // from the cern beam test paper
  const double EdepToPhotConv_FGD_3 = EdepToPhotConv_FGD_1_vert * (42.14/52.53); // from the cern beam test paper
  const double EdepToPhotConv_FGD_std_1_vert = 7.563; // from the cern beam test paper
  const double EdepToPhotConv_FGD_std_1_hori = 6.82; // from the cern beam test paper
  const double EdepToPhotConv_FGD_std_2 = EdepToPhotConv_FGD_std_1_vert * (51.56/52.53); // from the cern beam test paper
  const double EdepToPhotConv_FGD_std_3 = EdepToPhotConv_FGD_std_1_vert * (42.14/52.53); // from the cern beam test paper

  const double DistMPPCscint_FGD = 41; //*CLHEP::mm;
  //const double LongCompFrac_FGD = 0.816;
  const double LongCompFrac_FGD = 0.86; // from the cern beam test paper
  //const double LongAtt_FGD = 11926.; //*CLHEP::mm;
  //const double ShortAtt_FGD = 312.; //*CLHEP::mm;
  const double LongAtt_FGD = 400.0; //cm from cern beam test paper
  const double ShortAtt_FGD = 63.1; //cm from cern beam test paper
  const double DecayLength_FGD = 0.0858; // CLHEP::mm;
  const double Lbar_FGD = 1864.3; //* CLHEP::mm;
  const double TransTimeInFiber = 1./28. *10.; //mm/ns
  const double speedInFiber = 30;
  // SuperFGD constants
  const double MPPCEff_SuperFGD = 0.38;
  // Approximate collection factors from PDG2016 (Detectors and accelerators section)
  const double CollFactor_SingleClad = 0.06;
  //const double CollFactor_DoubleClad = 0.054; // from Licciardi's thesis  
  const double CollFactor_DoubleClad = 0.10;
  // cross talk effect: Cesar: for Tyvek is not so clear, measurements report 56% Tyvek, so you could migrate 56% of 2.8% which is: 1.6%
  const double crossFraction = 0.028;
  const double crossFractionTyvek = 0.016;

  // fiber half length
  const int xx = 12;
  const int yy = 4;
  const int zz = 24;

  const double Pedestal = 0;//145;  // pedeltal of ADC counts
  const double Gain = 10;  // Gain ADC counts of high gain channel
  const double LowGain  = 1;  // Gain ADC counts of low gain channel
  const double ElecNoise = 1.7;  // sigma of high gain electronics noise
  const double LowElecNoise = 1.2;  // sigma of low gain electronics noise
  const double PixelGainVari = 0.031;  // gain variation among pixels

  double a=0.;        // long attenuation component fraction
  double d=0.;        // distance MPPC-scint outside the bar
  double LongAtt=0.;  // long attenuation length
  double ShortAtt=0.; // short attenuation length
  double Ldecay=0.;   // decay length
  double Lbar=0.;     // bar length

  double hitLocation[3]={},hitPE[3]={},hitT[3]={};
  //for 180 should be (0.279, 0.736) and center for 0 should be (0.990,0.908);
  //double beamCenterX = 0.90; double beamCenterY = 0.90;
  //double beamCenterX = 0.5; double beamCenterY = 0.5;
  ///gps/pos/centre -3.0895007 -0.17907120 -25.0 cm
  double beamCenterX = -3.0895007;
  double beamCenterY = -0.17907120;
  // configure 0 
  double rotation_x = 0.01 ;
  double rotation_y = 0.0083;
  double offsetX = 0;
  double offsetY = 0.1;
  double detHalfSizeZ = 24.; 

  Float_t neutronHitX[1000]; Float_t neutronHitY[1000]; Float_t neutronHitZ[1000];
  Float_t neutronHitE[1000]; Float_t neutronHitPDG[1000]; Float_t initialE[1000];
  Float_t vtx[3];
  Float_t neutronHitT[1000];
  Int_t delay[1000]; Int_t id[1000];
  Float_t neutron_time[1000]; Float_t gamma_time[1000];
  Float_t t0[1000]; Float_t spill[1000];
  Float_t neutronParentPDG[1000]; Float_t neutronParentId[1000];
  Float_t fsPx[100], fsPy[100], fsPz[100];
  Float_t neutronPx[1000], neutronPy[1000], neutronPz[1000], nElastic[100];
  Float_t pointPositionX[100], pointPositionY[100], pointPositionZ[100];
  Int_t pointProcess[100];

  TGraph* bFlux = new TGraph(41);
  bFlux->SetPoint(1, 0 ,  1.26E-001);
  bFlux->SetPoint(2, 1.42E-001 ,  1.35E-001);
  bFlux->SetPoint(3, 1.79E-001 ,  1.35E-001);
  bFlux->SetPoint(4, 2.25E-001 ,  1.40E-001);
  bFlux->SetPoint(5, 2.84E-001 ,  1.44E-001);
  bFlux->SetPoint(6, 3.57E-001 ,  1.42E-001);
  bFlux->SetPoint(7, 4.50E-001 ,  1.41E-001);
  bFlux->SetPoint(8, 5.66E-001 ,  1.42E-001);
  bFlux->SetPoint(9, 7.13E-001 ,  1.30E-001);
  bFlux->SetPoint(10,8.97E-001 ,  1.17E-001);
  bFlux->SetPoint(11,1.13E+000 ,  9.66E-002);
  bFlux->SetPoint(12,1.42E+000 ,  7.64E-002);
  bFlux->SetPoint(13,1.79E+000 ,  6.59E-002);
  bFlux->SetPoint(14,2.25E+000 ,  5.34E-002);
  bFlux->SetPoint(15,2.84E+000 ,  4.27E-002);
  bFlux->SetPoint(16,3.57E+000 ,  3.36E-002);
  bFlux->SetPoint(17,4.50E+000 ,  2.73E-002);
  bFlux->SetPoint(18,5.66E+000 ,  1.97E-002);
  bFlux->SetPoint(19,7.13E+000 ,  1.39E-002);
  bFlux->SetPoint(20,8.97E+000 ,  9.51E-003);
  bFlux->SetPoint(21,1.13E+001 ,  6.26E-003);
  bFlux->SetPoint(22,1.42E+001 ,  4.00E-003);
  bFlux->SetPoint(23,1.79E+001 ,  2.61E-003);
  bFlux->SetPoint(24,2.26E+001 ,  2.13E-003);
  bFlux->SetPoint(25,2.84E+001 ,  1.59E-003);
  bFlux->SetPoint(26,3.57E+001 ,  1.30E-003);
  bFlux->SetPoint(27,4.50E+001 ,  1.05E-003);
  bFlux->SetPoint(28,5.66E+001 ,  8.59E-004);
  bFlux->SetPoint(29,7.13E+001 ,  7.47E-004);
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

  t.SetBranchAddress("fsPx",&fsPx);
  t.SetBranchAddress("fsPy",&fsPy);
  t.SetBranchAddress("fsPz",&fsPz);  
  t.SetBranchAddress("neutronHitX",&neutronHitX);
  t.SetBranchAddress("neutronHitY",&neutronHitY);
  t.SetBranchAddress("neutronHitZ",&neutronHitZ);
  t.SetBranchAddress("neutronHitE",&neutronHitE);
  t.SetBranchAddress("neutronHitT",&neutronHitT);
  t.SetBranchAddress("neutronHitPDG",&neutronHitPDG);
  t.SetBranchAddress("initialE",&initialE);
  t.SetBranchAddress("vtx",&vtx);
  t.SetBranchAddress("delay",&delay);
  t.SetBranchAddress("t0",&t0);
  t.SetBranchAddress("id",&id);
  t.SetBranchAddress("gamma_time",&gamma_time);
  t.SetBranchAddress("neutron_time",&neutron_time);
  t.SetBranchAddress("spill",&spill);
  t.SetBranchAddress("neutronParentPDG",&neutronParentPDG);
  t.SetBranchAddress("neutronParentId",&neutronParentId);
  t.SetBranchAddress("neutronPx",&neutronPx);
  t.SetBranchAddress("neutronPy",&neutronPy);
  t.SetBranchAddress("neutronPz",&neutronPz);
  t.SetBranchAddress("nElastic",&nElastic);
  t.SetBranchAddress("pointPositionX",&pointPositionX);
  t.SetBranchAddress("pointPositionY",&pointPositionY);
  t.SetBranchAddress("pointPositionZ",&pointPositionZ);
  t.SetBranchAddress("pointProcess",&pointProcess);

  int TotalSpills=0;
  int MinSpills=1e9;
  int febx=0;

  int cubeFlag[101]={};
  int cubeFlag2[101]={}; 
  int cubeFlag3[101]={};
  int cubeFlag4[10000]={};
  int cubeFlagXY[111][111]={};
  int cubeFlagXZ[111][111]={};
  int cubeFlagYZ[111][111]={};
  double cubeFlagXY_c[111][111]={};
  double cubeFlagXZ_c[111][111]={};
  double cubeFlagYZ_c[111][111]={};
  double cubeFlagXY_t[111][111]={};
  double cubeFlagXZ_t[111][111]={};
  double cubeFlagYZ_t[111][111]={};  
  int timeFlagXY[30000]={};
  int timeFlagXZ[30000]={};
  int timeFlagYZ[30000]={};
  int timeFlag[30000]={};

  TRandom3* grandom = new TRandom3();

  ostringstream sFEBnum;
  string sFEB;

  vector<int> FEBnumbers;
  FEBnumbers.clear();

  //AllEvents tree
  TTree AllEvents("AllEvents", "The ROOT tree of events");
  Event* event = new Event();
  AllEvents.Branch("Event", "Event", event);

  Int_t eventNum=0;
  bool SpillMissed = false;
  Int_t XZocc[48], ZYocc[48], EDepXZ[48],EDepZY[48],EDep[48];

  Double_t vtxSave[3]={};
  Double_t initialESave=0;

  double adc_tmp[3]={},loadc_tmp[3]={},Q[3]={},loQ[3]={},adc[3]={};
  double loadc[3]={};
  int aveXSave=-1000; int aveYSave=-1000; int aveZSave=-1000;
  int aveX=-1000; int aveY=-1000; int aveZ=-1000;

  int counter1 = 0;
  int counter2 = 0;
  int counter3 = 0;
  int counter4 = 0;
  int saveCount[1000]={};
  bool somethingIn = false;
  bool notOnlyN = false;  
  bool onlyP = true;
  int flagST = 0;

  TFile wfile("output.root","recreate");
  //TFile wfile(Form("test_r%f.root",atof(argv[1])), "recreate");

  //loop over spills
  //for (Int_t subSpill = atof(argv[2]) * 0.01 * t.GetEntries() ; subSpill< (atof(argv[2])+1)* 0.01 * t.GetEntries(); subSpill++) {
  for(Int_t subSpill = 0; subSpill< t.GetEntries() ; subSpill++){

    t.GetEntry(subSpill); 
    //cout<<"spill "<<subSpill<<endl;
    if(subSpill%1000 == 0 ) std::cout<<"event "<<subSpill<<" vertex x and y : "<<vtx[0]<<" "<<vtx[1]<<" "<<sqrt((vtx[0] - beamCenterX)*(vtx[0] - beamCenterX)+(vtx[1]- beamCenterY+1)*(vtx[1]-beamCenterY+1) )<<std::endl;

    // for the large detector
    //if(sqrt((vtx[0] - beamCenterX)*(vtx[0] - beamCenterX)+(vtx[1]- beamCenterY+1)*(vtx[1]-beamCenterY+1) )<0.4)
    // for the normal size detector
    if(sqrt((vtx[0] - beamCenterX)*(vtx[0] - beamCenterX)+(vtx[1]- beamCenterY)*(vtx[1]-beamCenterY) )<atof(argv[1]))
    {	  
    Int_t micropulse=0;
    double timeSave = 1e10;
    double edep = 0;
    aveXSave = -1000; aveYSave =-1000; aveZSave=-1000;
    aveX = -1100; aveY = -1100; aveZ = -1100;
    int locSave = 0;

    somethingIn = false;
    notOnlyN = false;
    onlyP = true;
    int occ[100]={};
    int occ2[100]={};
    int occ3[100]={};
    counter1 =0;
    counter3 =0;
    counter4 =0;
    flagST = 0;
    int locCounter = 0;
    for (int veryTemp =0; veryTemp<1000; veryTemp ++){
      neutronHitX[veryTemp] = neutronHitX[veryTemp];// + offsetX + TMath::Sin(rotation_x)*(neutronHitZ[veryTemp]+detHalfSizeZ);
      neutronHitY[veryTemp] = neutronHitY[veryTemp];// + offsetY + TMath::Sin(rotation_y)*(neutronHitZ[veryTemp]+detHalfSizeZ);
      neutronHitZ[veryTemp] = neutronHitZ[veryTemp];//TMath::Cos(rotation_x) * (neutronHitZ[veryTemp] + detHalfSizeZ) - detHalfSizeZ;
    }

    for(int lll=0;lll<101;lll++){
      for(int lll2=0;lll2<101;lll2++){
        cubeFlagXY[lll][lll2] = 0;	      
        cubeFlagXZ[lll][lll2] = 0;
	cubeFlagYZ[lll][lll2] = 0;
        cubeFlagXY_c[lll][lll2] = 0;
        cubeFlagXZ_c[lll][lll2] = 0;
        cubeFlagYZ_c[lll][lll2] = 0;
        cubeFlagXY_t[lll][lll2] = 1e10;
        cubeFlagXZ_t[lll][lll2] = 1e10;
        cubeFlagYZ_t[lll][lll2] = 1e10;	
      }      
    }

    for(int lll=0;lll<101;lll++){
      cubeFlag[lll] = 0;
      cubeFlag2[lll] = 0;
      cubeFlag3[lll] = 0;
    }
    for(int lll=0;lll<10000;lll++){
      cubeFlag4[lll] = 0;
      timeFlagXY[lll] = 0;
      timeFlagXZ[lll] = 0;
      timeFlagYZ[lll] = 0;
      timeFlag[lll] = 0;
    }

    double mark[101]={};
    for(int j=0;j<200;j++){
      if(id[j]>=0){
        if ( neutronParentPDG[j] == 2112 && neutronHitPDG[j] == 2212){
	//if ( neutronHitPDG[j] == 211){
          occ[id[j]] = 1;
        }
        else{
	  if(neutronHitE[j]>0.1)
	    occ3[id[j]] = 1;
	  if(neutronHitE[j]>2){	  
            int tempz = (int)neutronHitZ[j]+50;
	    if(tempz< 100){
              if (mark[tempz] == 0) 
	        occ2[id[j]] = 1;
	      mark[tempz] = 1;
	    }  
	  }  
        }
      }
    }

    for(int j=0;j<100;j++){
      if(occ[j] == 1) counter1 ++;
      if(occ2[j] == 1 ) counter3 ++;
      if(occ3[j] == 1 ) counter4 ++;
    }
    for(int j=0;j<100;j++){
      if( neutronHitPDG[j] != 2112 && neutronHitPDG[j]> 0) notOnlyN = true;
      if( neutronHitPDG[j] != 2212 ) onlyP = false;
    }

    for(int j=0;j<100;j++){
      if(neutronHitX[j]<12 && neutronHitX[j]>-12 && neutronHitY[j]<4 && neutronHitY[j]>-4 && neutronHitZ[j]<-2 && neutronHitZ[j]>-50)
      {
        somethingIn = true;
      }
    }

    if ( counter1 == 1 && counter4 < 1 )
      flagST = 1 ;
    else if (counter3 > 3)
      flagST = 2 ;    
    else
      flagST = 0 ;
    
    double current_time[10];
    std::vector<double> true_vtx;
    for (int ivtx = 0;ivtx<10;ivtx++){
      true_vtx.push_back(-999);
      current_time[ivtx] = 1e9;
    }  

    for(Int_t hitLoop=0;hitLoop<200;hitLoop++){

    for(int ivtx=0;ivtx<10;ivtx++){
      if (neutronHitE[hitLoop]>0.1*ivtx && neutronHitT[hitLoop]< current_time[ivtx]){
        current_time[ivtx] = neutronHitT[hitLoop];
        true_vtx[ivtx] = neutronHitZ[hitLoop];
      }
    }

    if (neutronHitE[hitLoop] < 0.001) continue;
    double hh1 = ((neutronHitX[hitLoop]+ cubeSize*12. ) + 0.0001)/cubeSize;
    double hh2 = ((neutronHitY[hitLoop]+ cubeSizeY*4. )  + 0.0001)/cubeSize;
    double hh3 = ((neutronHitZ[hitLoop]+ cubeSize*24. ) + 0.0001)/cubeSize;
    if (hh3<0) {cout<<"ah!! "<<neutronHitZ[hitLoop]<<" "<<hh3<<endl;  continue;} //exit(0);}
    aveX = hh1; aveY = hh2; aveZ = hh3;

    edep = neutronHitE[hitLoop];
    double Nphot_y; 

    if (!IsRotated){
      if (hh3 <= 15 )  Nphot_y = edep * (EdepToPhotConv_FGD_2 + gRandom->Gaus(0,EdepToPhotConv_FGD_std_2));
      else if (hh3 >15 && hh3 <=23)  Nphot_y = edep * (EdepToPhotConv_FGD_3 + gRandom->Gaus(0,EdepToPhotConv_FGD_std_3));
      else  Nphot_y = edep * (EdepToPhotConv_FGD_1_vert + gRandom->Gaus(0,EdepToPhotConv_FGD_std_1_vert));
    }
    else{
      if (hh3 <= 23 )  Nphot_y = edep *( EdepToPhotConv_FGD_1_vert + gRandom->Gaus(0,EdepToPhotConv_FGD_std_1_vert));
      else if (hh3 >23 && hh3 <=31)  Nphot_y = edep *( EdepToPhotConv_FGD_3 + gRandom->Gaus(0,EdepToPhotConv_FGD_std_3));
      else  Nphot_y = edep * (EdepToPhotConv_FGD_2 + gRandom->Gaus(1,EdepToPhotConv_FGD_std_2));
    }

    double Nphot_x = edep * (EdepToPhotConv_FGD_1_hori + gRandom->Gaus(0,EdepToPhotConv_FGD_std_1_hori));
    double Nphot_z = edep * (EdepToPhotConv_FGD_1_hori + gRandom->Gaus(0,EdepToPhotConv_FGD_std_1_hori));

    a = LongCompFrac_FGD;
    LongAtt = LongAtt_FGD;
    ShortAtt = ShortAtt_FGD;

    double xNphot; double yNphot; double zNphot;
    double xTimeDelay; double yTimeDelay; double zTimeDelay;
    if(timeSave == 1e10) timeSave = neutronHitT[hitLoop];

    // was 0, set to 1 on Jun 28 2020
    if (aveX%2 == 0) {
      zNphot = Nphot_z * ( a*exp((-(TMath::Abs(aveZ)))/LongAtt) + (1-a)*exp((-(TMath::Abs(aveZ)))/ShortAtt) );
      zTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveZ);
    }
    else {
      zNphot = Nphot_z * ( a*exp((-(TMath::Abs(cubeSize*2*zz- aveZ)))/LongAtt) + (1-a)*exp((-(TMath::Abs(cubeSize*2*zz- aveZ)))/ShortAtt) );
      zTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(cubeSize*2*zz- aveZ);
    }
    if (aveZ%2 == 0) {
      yNphot = Nphot_y * ( a*exp((-(TMath::Abs(aveY)))/LongAtt) + (1-a)*exp((-(TMath::Abs(aveY)))/ShortAtt) );
      xNphot = Nphot_x * ( a*exp((-(TMath::Abs(aveX)))/LongAtt) + (1-a)*exp((-(TMath::Abs(aveX)))/ShortAtt) );
      yTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveY);
      xTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveX);
    }
    else{
      yNphot = Nphot_y * ( a*exp((-(TMath::Abs(cubeSize*2*yy - aveY)))/LongAtt) + (1-a)*exp((-(TMath::Abs(cubeSize*2*yy -aveY)))/ShortAtt) );
      xNphot = Nphot_x * ( a*exp((-(TMath::Abs(cubeSize*2*xx - aveX)))/LongAtt) + (1-a)*exp((-(TMath::Abs(cubeSize*2*xx -aveX)))/ShortAtt) );
      yTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(cubeSize*2*yy -aveY);
      xTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(cubeSize*2*xx -aveX);
    }
   
    //cout<<"::::energy check point:::: ====== >> energy deposit : "<<edep<<" x fiber readout  "<<xNphot<<" y fiber readout "<<yNphot<<" z fiber readout "<<zNphot<<endl; 
    //cout<<"checking more... aveX aveY aveZ : "<<aveX<<" "<<aveY<<" "<<aveZ<<"  Nphot_x Nphot_y Nphot_z : "<<Nphot_x<<" "<<Nphot_y<<" "<<Nphot_z<<endl;
    // no need for MPPC effect, already included in the LY
    double xpe = xNphot;
    double ype = yNphot;
    double zpe = zNphot;

    // neutron_time : spill time + tof
    // subtract spill starting time, remains tof + neutron scattering time
    hitT[0]=xTimeDelay + neutron_time[0] -t0[0];
    hitT[1]=yTimeDelay + neutron_time[0] -t0[0];
    hitT[2]=zTimeDelay + neutron_time[0] -t0[0];

    cubeFlagXY_c[aveX][aveY] += zpe * (1- 4*crossFraction - 2*crossFractionTyvek);
    if(hitT[2]<cubeFlagXY_t[aveX][aveY])
      cubeFlagXY_t[aveX][aveY] = hitT[2];

    cubeFlagXZ_c[aveX][aveZ] += ype * (1- 4*crossFraction - 2*crossFractionTyvek);
    if(hitT[1]<cubeFlagXZ_t[aveX][aveZ])
      cubeFlagXZ_t[aveX][aveZ] = hitT[1];

    cubeFlagYZ_c[aveY][aveZ] += xpe * (1- 4*crossFraction - 2*crossFractionTyvek);
    if(hitT[0]<cubeFlagYZ_t[aveY][aveZ])
      cubeFlagYZ_t[aveY][aveZ] = hitT[0];

    cubeFlagXY_c[aveX+1][aveY] += zpe * crossFraction;
    cubeFlagXZ_c[aveX+1][aveZ] += zpe * crossFraction;
    cubeFlagYZ_c[aveY][aveZ] += zpe * crossFraction;

    cubeFlagXY_c[aveX-1][aveY] += zpe * crossFraction;
    cubeFlagXZ_c[aveX-1][aveZ] += zpe * crossFraction;
    cubeFlagYZ_c[aveY][aveZ] += zpe * crossFraction;

    cubeFlagXY_c[aveX][aveY] += zpe * crossFraction;
    cubeFlagXZ_c[aveX][aveZ+1] += zpe * crossFraction;
    cubeFlagYZ_c[aveY][aveZ+1] += zpe * crossFraction;

    cubeFlagXY_c[aveX][aveY] += zpe * crossFraction;
    cubeFlagXZ_c[aveX][aveZ-1] += zpe * crossFraction;
    cubeFlagYZ_c[aveY][aveZ-1] += zpe * crossFraction;

    cubeFlagXY_c[aveX][aveY+1] += zpe * crossFractionTyvek;
    cubeFlagXZ_c[aveX][aveZ] += zpe * crossFractionTyvek;
    cubeFlagYZ_c[aveY+1][aveZ] += zpe * crossFractionTyvek;

    cubeFlagXY_c[aveX][aveY-1] += zpe * crossFractionTyvek;
    cubeFlagXZ_c[aveX][aveZ] += zpe * crossFractionTyvek;
    cubeFlagYZ_c[aveY-1][aveZ] += zpe * crossFractionTyvek;

    cubeFlagXY_t[aveX+1][aveY] = cubeFlagXY_t[aveX][aveY];
    cubeFlagXY_t[aveX-1][aveY] = cubeFlagXY_t[aveX][aveY];
    cubeFlagXY_t[aveX][aveY+1] = cubeFlagXY_t[aveX][aveY];
    cubeFlagXY_t[aveX][aveY-1] = cubeFlagXY_t[aveX][aveY];

    cubeFlagXZ_t[aveX+1][aveZ] = cubeFlagXZ_t[aveX][aveZ];
    cubeFlagXZ_t[aveX-1][aveZ] = cubeFlagXZ_t[aveX][aveZ];
    cubeFlagXZ_t[aveX][aveZ+1] = cubeFlagXZ_t[aveX][aveZ];
    cubeFlagXZ_t[aveX][aveZ-1] = cubeFlagXZ_t[aveX][aveZ];

    cubeFlagYZ_t[aveY+1][aveZ] = cubeFlagYZ_t[aveY][aveZ];
    cubeFlagYZ_t[aveY-1][aveZ] = cubeFlagYZ_t[aveY][aveZ];
    cubeFlagYZ_t[aveY][aveZ+1] = cubeFlagYZ_t[aveY][aveZ];
    cubeFlagYZ_t[aveY][aveZ-1] = cubeFlagYZ_t[aveY][aveZ];

    }

    timeSave = 1e10;
    for(int loop1=0;loop1<24;loop1++){
      for(int loop2=0;loop2<8;loop2++){
	
	if(cubeFlagXY_c[loop1][loop2] < 0.01) continue;
        Hit* hit = event->AddHit();
	hit->SetView(0);
	hit->SetX(loop1);
	hit->SetY(loop2);
	hit->SetZ(-1);
        hit->SetPE(cubeFlagXY_c[loop1][loop2]);
	hit->SetDt(cubeFlagXY_t[loop1][loop2]/tickSize);
	hit->SetMC(true);
      }
    }
    for(int loop1=0;loop1<24;loop1++){
      for(int loop2=0;loop2<48;loop2++){
        if(cubeFlagXZ_c[loop1][loop2] < 0.01) continue;
        Hit* hit = event->AddHit();
        hit->SetView(1);
        hit->SetX(loop1);
        hit->SetY(-1);
        hit->SetZ(loop2);	
        hit->SetPE(cubeFlagXZ_c[loop1][loop2]);
        hit->SetDt(cubeFlagXZ_t[loop1][loop2]/tickSize);
        hit->SetMC(true);
      }
    }
    for(int loop1=0;loop1<8;loop1++){
      for(int loop2=0;loop2<48;loop2++){
        if(cubeFlagYZ_c[loop1][loop2] < 0.01) continue;
        Hit* hit = event->AddHit();
        hit->SetView(2);
        hit->SetX(-1);
        hit->SetY(loop1);
        hit->SetZ(loop2);	
        hit->SetPE(cubeFlagYZ_c[loop1][loop2]);
        hit->SetDt(cubeFlagYZ_t[loop1][loop2]/tickSize);
        hit->SetMC(true);
      }
    }

    // five digits: first two: baseline  last three: rotation
    event->SetEventID(90000);
    event->SetFEB12ch(0);
    event->SetFEB12LeadTime(0);
    event->SetFEB12hitTfromSpill(0);
    event->SetMicropulse(micropulse);
    //cout<<flagST<<endl;
    event->SetRange(-999);
    if(flagST == 1)
      event->SetRange(0);
    if(flagST == 2)
      event->SetRange(1);
    event->SetMaxCharge(event->FindMaxCharge());
    event->SetEventID(eventNum);
    event->SetEmptyZ(0);
    event->SetdEdz(0,fsPx[0]);
    event->SetdEdz(1,fsPy[0]);
    event->SetdEdz(2,fsPz[0]);
    event->SetdEdz(3,neutronPx[0]);
    event->SetdEdz(4,neutronPy[0]);
    event->SetdEdz(5,neutronPz[0]);
    for(int inloop=0 ;inloop<100 ;inloop++){
      event->SetnElastic(inloop, nElastic[inloop]);
      event->SetPointPosition(inloop, 0, pointPositionX[inloop]);
      event->SetPointPosition(inloop, 1, pointPositionY[inloop]);
      event->SetPointPosition(inloop, 2, pointPositionZ[inloop]);
      event->SetPointPosition(inloop, 3, pointProcess[inloop]);
    }

    for(int ivt=0;ivt<1000;ivt++){
    //std::cout<<":::::::::::::::::  "<<trueLocationX[ivt]<<" "<<trueLocationY[ivt]<<" "<<trueLocationZ[ivt]<<std::endl;
      event->SetTrueLocation(ivt, 0, (double)neutronHitX[ivt]);
      event->SetTrueLocation(ivt, 1, (double)neutronHitY[ivt]);
      event->SetTrueLocation(ivt, 2, (double)neutronHitZ[ivt]);
      event->SetTrueLocation(ivt, 3, (double)neutronHitT[ivt]);
    }

    event->SetColmt(0, vtx[0]);
    event->SetColmt(1, vtx[1]);
    for(int ivt=0;ivt<10;ivt++)
      event->SetTrueVtx(ivt,true_vtx[ivt]);
    //cout<<"initial momentum along x y z : "<<event->GetdEdz(0)<<" "<<event->GetdEdz(1)<<" "<<event->GetdEdz(2)<<endl;

    // SetXStdDev is the neutron true energy
    event->SetXStdDev(initialE[0]);

    // SetYStdDev is the flux weight
    //cout<<"flux weight "<<bFlux->Eval(initialE[0])<<endl;
    event->SetYStdDev(bFlux->Eval(initialE[0]));

    AllEvents.Fill();

    event->Clear();
    eventNum++;
    micropulse++;
  }
  }
  wfile.cd();
  AllEvents.Write("",TObject::kOverwrite);
  //AllEvents.Delete();
  wfile.Close();
  //t.Delete();
  //return 0;
}
