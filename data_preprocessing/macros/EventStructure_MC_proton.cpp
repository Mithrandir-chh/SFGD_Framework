

#include "../src/global_header.hh"
#include <cstring>
#include "TMath.h"

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
  //t.Add("~/work/sfgd_framework/data_preprocessing/build/proton750/FHC_*_wFiber_wDead_realSize_wTyvek_noGap_realRot_0.root");
  //t.Add("~/work/sfgd_framework/data_preprocessing/build/neutron/FHC_*_24x8x48_wFiber_wDead_realSize_wTyvek_noGap_realRot_0_*.root");
  //t.Add("~/work/sfgd_framework/data_preprocessing/build/neutron/FHC_*_24x8x48_wFiber_wDead_realSize_wTyvek_noGap_traj3_centerX0.5_Y0.5.root");
  //t.Add("~/work/sfgd_framework/data_preprocessing/build/neutron/FHC_499_24x8x48_wFiber_wDead_realSize_wTyvek_noGap_realRot_0_traj3_centerX0.5_Y0.5.root");
  //t.Add("~/work/sfgd_framework/data_preprocessing/build/neutron/FHC_*_noDead_chunk_traj3_centerX0.5_Y0.5.root");
  t.Add("~/work/sfgd_framework/data_preprocessing/build/proton/FHC_3*_wFiber_wDead_realSize_wTyvek_noGap_realRot_0.root");

  //t.Add("~/work/sfgd_framework/data_preprocessing/build/neutron/FHC*100MeV.root");
  //t.Add("~/work/edep-sim/FHC_wFiber_wDead_realSize_wTyvek_noGap_realRot_0_sub*.root");
  //t.Add("~/work/sfgd_framework/data_preprocessing/build/pion2/FHC_*_wFiber_wDead_realSize_wTyvek_noGap_realRot_0.root");
  //const double cubeSize = 1.026;
  const double cubeSize = 1.;
  const double tickSize = 2.5;
  const double ifNeutronSample = false;

  //constants for energy calibration
  const double CBIRKS = 0.00208; // mm/MeV
  const double EdepToPhotConv_FGD = 70.8; // CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
  const double DistMPPCscint_FGD = 41; //*CLHEP::mm;
  const double LongCompFrac_FGD = 0.816;
  const double LongAtt_FGD = 11926.; //*CLHEP::mm;
  const double ShortAtt_FGD = 312.; //*CLHEP::mm;
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
  double beamCenterX = 0.5; double beamCenterY = 0.5;
  // configure 0 
  double rotation_x = 0.01 ;
  double rotation_y = 0.0083;
  double offsetX = 0;
  double offsetY = 0.1;
  double detHalfSizeZ = 24.; 
/*
  // configure 180
  double rotation_x = -0.0062;
  double rotation_y = -0.0062;
  double offsetX = 0.2;
  double offsetY = -0.2;
  double detHalfSizeZ = 24.
*/
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

  TFile wfile(Form("LANL_MC_proton_r%f.root",atof(argv[1])), "recreate");
  //TFile wfile("LANL_MC_normalSize_pencilBeam_100MeV.root", "recreate");
  //TFile wfile("LANL_MC_pion_flag.root","recreate");

  //loop over spills
  for (Int_t subSpill = 0; subSpill<t.GetEntries(); subSpill++) {

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
    int trueLocationX[1000] ;
    int trueLocationY[1000] ;
    int trueLocationZ[1000] ;
    double trueLocationT[1000] ;
    for (int veryTemp =0; veryTemp<1000; veryTemp ++){
      trueLocationX[veryTemp] = -100;
      trueLocationY[veryTemp] = -100;
      trueLocationZ[veryTemp] = -100;
      trueLocationT[veryTemp] = -100;

      neutronHitX[veryTemp] = neutronHitX[veryTemp];// + offsetX + TMath::Sin(rotation_x)*(neutronHitZ[veryTemp]+detHalfSizeZ);
      neutronHitY[veryTemp] = neutronHitY[veryTemp];// + offsetY + TMath::Sin(rotation_y)*(neutronHitZ[veryTemp]+detHalfSizeZ);
      neutronHitZ[veryTemp] = neutronHitZ[veryTemp];//TMath::Cos(rotation_x) * (neutronHitZ[veryTemp] + detHalfSizeZ) - detHalfSizeZ;
    }

    for(int lll=0;lll<101;lll++){
      for(int lll2=0;lll2<101;lll2++){
        cubeFlagXY[lll][lll2] = 0;	      
        cubeFlagXZ[lll][lll2] = 0;
	cubeFlagYZ[lll][lll2] = 0;
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
        if ( neutronHitPDG[j] == 2212){
	//if ( neutronHitPDG[j] == 211){
          occ[id[j]] = 1;
        }
        else{
	  if(neutronHitE[j]>0.1)
	    occ3[id[j]] = 1;
	  if(neutronHitE[j]>2 && neutronHitPDG[j] != 2112 ){	  
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
    /*
    if (counter3 > 3)
      flagST = 2 ; 
    else 
      flagST = 1;
    */
     
    //if( somethingIn == false || notOnlyN == false)
    //  continue;

    if(notOnlyN == false)
      continue;
    else if ( counter1 == 1 && counter4 < 1 )
      flagST = 1 ;
    else if (counter3 > 1)
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

    if (neutronHitE[hitLoop] < 0.0001) continue;
    int sameCube = 0;
    aveX = neutronHitX[hitLoop]/cubeSize+xx;
    aveY = neutronHitY[hitLoop]/cubeSize+yy;
    aveZ = neutronHitZ[hitLoop]/cubeSize+zz;
    aveX -= xx;
    aveY -= yy;
    aveZ -= zz;
    //cout<<"reference x y z "<<aveX<<" "<<aveY<<" "<<aveZ<<" "<<neutronHitE[hitLoop]<<endl;
    //if(neutronHitZ[hitLoop]<-49.99)
    //std::cout<<neutronHitX[hitLoop]<<" "<<neutronHitY[hitLoop]<<" "<<neutronHitZ[hitLoop]<<" "<<aveZ<<std::endl;
    //std::cout<<hitLoop<<" "<<neutronHitX[hitLoop]<<" "<<neutronHitY[hitLoop]<<" "<<neutronHitZ[hitLoop]<<" "<<neutronHitPDG[hitLoop]<<std::endl;
/*
    std::cout<<"--------------- compare -------------"<<endl;
    std::cout<<aveX<<" "<<aveY<<" "<<aveZ<<std::endl;
    std::cout<<aveXSave<<" "<<aveYSave<<" "<<aveZSave<<std::endl;
    std::cout<<"--------------- endl -------------"<<endl;
*/
    if ( (aveX == aveXSave && aveY == aveYSave && aveZ == aveZSave) ) sameCube = 1;
    if ( !(aveX == aveXSave && aveY == aveYSave && aveZ == aveZSave) && aveXSave!=-1000 ) sameCube =2;
    if (aveXSave == -1000 && aveYSave == -1000 && aveZSave == -1000) sameCube = 1;

    //std::cout<<neutronHitX[hitLoop]<<std::endl;
    //std::cout<<"aveX and aveXSave : "<<aveX<<" "<<aveXSave<<" aveY and aveYSave : "<<aveY<<" "<<aveYSave<< " aveZ and aveZSave : "<<aveZ<<" "<<aveZSave<<"  "<<neutronHitT[hitLoop]<<std::endl;

    //sameCube = 1;
    //cout<<"neutron hit time "<<neutronHitT[hitLoop]<<endl;
    if (sameCube ==1 ){ 	    
      edep += neutronHitE[hitLoop];
      //cout<<"in sameCube == 1 : " <<neutronHitT[hitLoop]<<" "<<timeSave<<endl;
      if (neutronHitT[hitLoop]< timeSave){
        timeSave = neutronHitT[hitLoop];
      }
    }	    
    std::cout<<"if same cube "<<sameCube<<std::endl;
    //loop over hits

    //sameCube = 2;
    if( sameCube ==2 )
    {      

      bool viewXYFlag = false;
      bool viewXZFlag = false;
      bool viewYZFlag = false;

      int veryTempT=0;
      if(ifNeutronSample)
        veryTempT = (neutron_time[0]+neutronHitT[locSave]-t0[0]) /tickSize;
      else
        veryTempT = (neutronHitT[locSave]-1) /tickSize;

      //cout<<aveXSave+24<<" "<<aveYSave+24<<" "<<aveZSave+50<<" "<<neutron_time[0]<<" "<<neutronHitT[locSave]<<" "<<t0[0]<<" "<<veryTempT<<endl;
      if (veryTempT < 0) veryTempT = 0; 
      //cout<<"testing.. "<<aveXSave+24<<" "<<aveYSave+24<<" "<<aveZSave+50<<" "<<cubeFlag[aveXSave+24]<<" "<<cubeFlag2[aveYSave+24]<<" "<<cubeFlag3[aveZSave+50]<<endl;
/*
      if (cubeFlag[aveXSave+24]==1 && cubeFlag2[aveYSave+24] == 1 && cubeFlag3[aveZSave+50] == 1 && cubeFlag4[veryTempT] == 1){
        aveXSave = aveX;
        aveYSave = aveY;
        aveZSave = aveZ;
        locSave = hitLoop;
        edep = neutronHitE[hitLoop];
        timeSave = 1e10;
        continue;
      }
      //cout<<locSave<<endl;
*/
      //cout<<"-------------------------- "<<aveXSave+24<<" "<<aveYSave+24<<" "<<veryTempT<<endl;
/*
      if (cubeFlagXY[aveXSave+24][aveYSave+24]==1 && cubeFlag4[veryTempT] == 1){
        viewXYFlag = true;
	//cout<<"!!!!!!!!!!!!!!!!!!!!! caught one overlapping hit "<<endl;
      }

      if (cubeFlagXZ[aveXSave+24][aveZSave+50]==1 && cubeFlag4[veryTempT] == 1){
        viewXZFlag = true;
      }

      if (cubeFlagYZ[aveYSave+24][aveZSave+50]==1 && cubeFlag4[veryTempT] == 1){
        viewYZFlag = true;
      }

      if(aveXSave+24<100 && aveYSave+24<100 && aveZSave+50< 101 && veryTempT< 10000){
        cubeFlag[aveXSave+24] = 1;
        cubeFlag2[aveYSave+24] = 1;
        cubeFlag3[aveZSave+50] = 1;
        cubeFlag4[veryTempT]  = 1;  
        cubeFlagXY[aveXSave+24][aveYSave+24] =1 ;
	cubeFlagXZ[aveXSave+24][aveZSave+50] =1 ;
	cubeFlagYZ[aveYSave+24][aveZSave+50] =1 ;
      }
*/
      //else{
        //cout<<aveXSave+24<<" "<<aveYSave+24<<" "<<aveZSave+50<<" "<<neutron_time[0]<<" "<<neutronHitT[locSave]<<" "<<t0[0]<<" "<<veryTempT<<endl;
      //}

      double collfact = CollFactor_DoubleClad;
      double fact_fib1 = collfact;
      double fact_fib2 = (1-fact_fib1)*collfact;
      double fact_fib3 = (1-fact_fib2)*collfact;
      double CollFactAve = (fact_fib1+fact_fib2+fact_fib3)/3.;
      double NormShadowLight = CollFactAve / collfact; // fraction 
      double Nphot = edep * EdepToPhotConv_FGD * NormShadowLight;
      //cout<<edep<<" "<<EdepToPhotConv_FGD <<" "<< NormShadowLight<<" "<<Nphot<<endl;

      if(Nphot < 0.01) continue;

      a = LongCompFrac_FGD;
      d = DistMPPCscint_FGD;
      LongAtt = LongAtt_FGD;
      ShortAtt = ShortAtt_FGD;
      Ldecay= DecayLength_FGD;
      Lbar = Lbar_FGD;

      double xNphot; double yNphot; double zNphot;
      double xTimeDelay; double yTimeDelay; double zTimeDelay;

      if(timeSave == 1e10) timeSave = neutronHitT[hitLoop];

      // was 0, set to 1 on Jun 28 2020
      if ((int)neutronHitX[locSave]%2 == 0) {
        xNphot = Nphot * ( a*exp((-(TMath::Abs(aveXSave-xx))-d)/LongAtt) + (1-a)*exp((-(TMath::Abs(aveXSave-xx))-d)/ShortAtt) );
        yNphot = Nphot * ( a*exp((-(TMath::Abs(aveYSave-yy))-d)/LongAtt) + (1-a)*exp((-(TMath::Abs(aveYSave-yy))-d)/ShortAtt) );
        zNphot = Nphot * ( a*exp((-(TMath::Abs(aveZSave-zz))-d)/LongAtt) + (1-a)*exp((-(TMath::Abs(aveZSave-zz))-d)/ShortAtt) );
        xTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveXSave-xx);
        yTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveYSave-yy);
        zTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveZSave-zz);
      }
      else {
        xNphot = Nphot * ( a*exp((-(TMath::Abs(aveXSave-xx))-d)/LongAtt) + (1-a)*exp((-(TMath::Abs(aveXSave-xx))-d)/ShortAtt) );
        yNphot = Nphot * ( a*exp((-(TMath::Abs(aveYSave-yy))-d)/LongAtt) + (1-a)*exp((-(TMath::Abs(aveYSave-yy))-d)/ShortAtt) );
        zNphot = Nphot * ( a*exp((-(TMath::Abs(aveZSave+zz))-d)/LongAtt) + (1-a)*exp((-(TMath::Abs(aveZSave+zz))-d)/ShortAtt) );
        xTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveXSave-xx);
        yTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveYSave-yy);
        zTimeDelay =  timeSave+(1./speedInFiber) * TMath::Abs(aveZSave+zz);
      }

      double xpe = xNphot * MPPCEff_SuperFGD;
      double ype = yNphot * MPPCEff_SuperFGD;
      double zpe = zNphot * MPPCEff_SuperFGD;

      //cout<<"1 "<<zTimeDelay<<" "<<timeSave<<" "<<t0[0]<<endl;

      hitT[0]=xTimeDelay + neutron_time[0] ;
      hitT[1]=yTimeDelay + neutron_time[0] ;
      hitT[2]=zTimeDelay + neutron_time[0] ;
      //hitT[0] = timeSave + neutron_time[0];
      //hitT[1] = timeSave + neutron_time[0];
      //hitT[2] = timeSave + neutron_time[0];

      hitPE[0]=xpe;
      hitPE[1]=ype;
      hitPE[2]=zpe;
      //cout<<"test 1 "<<endl;
      /////////////////////////////////////////////////////////////////////////////////
      // checking the duplicate of single fiber
      /////////////////////////////////////////////////////////////////////////////////
     
      //if(timeFlag[veryTempT] == 1 && cubeFlagXY[aveXSave+24][aveYSave+24]==1  && cubeFlagXZ[aveXSave+24][aveZSave+50]==1 && cubeFlagYZ[aveYSave+24][aveZSave+50]==1)
      if( cubeFlagXY[aveXSave+24][aveYSave+24]==1 && cubeFlagXZ[aveXSave+24][aveZSave+50]==1 && cubeFlagYZ[aveYSave+24][aveZSave+50]==1)
      {
        aveXSave = aveX;
        aveYSave = aveY;
        aveZSave = aveZ;
        locSave = hitLoop;
        edep = neutronHitE[hitLoop];
        timeSave = 1e10;

      	continue;
      }
      //cout<<"test 2"<<endl;
      //cout<<" ... "<<aveXSave+24<<" "<<aveYSave+24<<" "<<veryTempT<<endl;
      if(veryTempT >= 20000) veryTempT = 0;
      if (cubeFlagXY[aveXSave+24][aveYSave+24]==1 && timeFlagXY[veryTempT] == 1){
        viewXYFlag = true;
        //cout<<"!!!!!!!!!!!!!!!!!!!!! caught one overlapping hit "<<endl;
      }

      if (cubeFlagXZ[aveXSave+24][aveZSave+50]==1 && timeFlagXZ[veryTempT] == 1){
        viewXZFlag = true;
      }

      if (cubeFlagYZ[aveYSave+24][aveZSave+50]==1 && timeFlagYZ[veryTempT] == 1){
        viewYZFlag = true;
      }
      
      // y = a + b * x ; x = PE , y = tot tick, a = 8.8, b = 0.2; above 80 pe, all 25 ticks and below 6 pe, all 10 ticks
      int durationX = 8.8 + 0.2 * hitPE[0];
      int durationY = 8.8 + 0.2 * hitPE[1];
      int durationZ = 8.8 + 0.2 * hitPE[2];

      cout<<"duartion of X Y Z fibers : "<<durationX<<" "<<durationY<<" "<<durationZ<<endl;

      if(aveXSave+24<100 && aveYSave+24<100 && aveZSave+50< 101 && veryTempT< 10000){
        //cubeFlag[aveXSave+24] = 1;
        //cubeFlag2[aveYSave+24] = 1;
        //cubeFlag3[aveZSave+50] = 1;
        //cubeFlag4[veryTempT]  = 1;  
        cubeFlagXY[aveXSave+24][aveYSave+24] =1 ;
        cubeFlagXZ[aveXSave+24][aveZSave+50] =1 ;
        cubeFlagYZ[aveYSave+24][aveZSave+50] =1 ;
        timeFlag[veryTempT] = 1;

        for(int lll=0;lll<durationX; lll++)
	  timeFlagYZ[veryTempT+lll] = 1;
        for(int lll=0;lll<durationY; lll++)
          timeFlagXZ[veryTempT+lll] = 1;
        for(int lll=0;lll<durationZ; lll++)
          timeFlagXY[veryTempT+lll] = 1;      
      }      

      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////

      std::cout<<"In event "<<subSpill<<" with X Y Z : "<<aveXSave<<" "<<aveYSave<<" "<<aveZSave<<" "<<edep<<std::endl;
      trueLocationX[locCounter] = aveXSave;
      trueLocationY[locCounter] = aveYSave;
      trueLocationZ[locCounter] = aveZSave;
      trueLocationT[locCounter] = neutron_time[0]+neutronHitT[locSave]-t0[0];

      locCounter ++;

      for(Int_t dim =0;dim<3;dim++){
        //PE to ADC
        adc_tmp[dim] = Pedestal + (hitPE[dim])*Gain;
        loadc_tmp[dim] = Pedestal + (hitPE[dim])*LowGain*14.29/13.55;

        //Electronics noise
        adc_tmp[dim] = grandom->Gaus(adc_tmp[dim],ElecNoise);
        loadc_tmp[dim] = grandom->Gaus(loadc_tmp[dim],LowElecNoise);

        //ADC to Charge
        //Q=(adc_tmp+53)/217;
        //loQ=(loadc_tmp+82)/26;
        Q[dim]=(adc_tmp[dim])/135.5;
        loQ[dim]=(loadc_tmp[dim])/14.29;

        //Non linearlity of high gain ADC
        if(Q[dim]<0.65) adc[dim]=135.5*Q[dim];
        else if(Q[dim]<3.2)  adc[dim]=217*Q[dim]-53;
        else if(Q[dim]<4.2)  adc[dim]=158.6*Q[dim]+133.9;
        else if(Q[dim]<14)  adc[dim]=5.1*Q[dim]+778.6;
        else  adc[dim]=850;

        //Non linearlity of low gain ADC
        if(loQ[dim]<7)  loadc[dim]=14.29*loQ[dim];
        else if(loQ[dim]<27)  loadc[dim]=26*loQ[dim]-82;
        else if(loQ[dim]<35.5)  loadc[dim]=21.12*loQ[dim]+48.24;
        else if(loQ[dim]<178.4)  loadc[dim]=0.7*loQ[dim]+775.1;
        else  loadc[dim]=900;

        Hit* hit = event->AddHit();

        hit->SetPE(hitPE[dim]*(1-crossFraction));
        hit->SetHG_pe(hitPE[dim]*(1-crossFraction));
        hit->SetLG_pe(hitPE[dim]*(1-crossFraction));
        hit->SetToT_pe(hitPE[dim]*(1-crossFraction));
        hit->SetHG_ADC(adc[dim]);
        hit->SetLG_ADC(loadc[dim]);
        hit->SetRE(0);
        hit->SetFE(0);
        hit->SetToT(0);
	double tempT = (hitT[dim]-t0[0])/tickSize;
        hit->SetDt(tempT);

	//cout<<hitT[dim]<<" "<<t0[0]<<endl;

	hit->SetSpillTag(spill[0]);
        // spilltime set to true time
        hit->SetSpillTime(hitT[dim]-t0[0]);

        hit->SetSpillTrailTime(0);
        hit->SetTfromSpill(neutron_time[0]);
        hit->SetFEB(0);
        hit->SetCh(0);
        hit->SetGTrigTag((int)hitT[dim]/tickSize);
        hit->SetGTrigTime((int)hitT[dim]/tickSize);

	if(dim == 2){
          hit->SetView(0);		
          hit->SetX(aveXSave);
          hit->SetY(aveYSave);
          hit->SetZ(-1);
	  if (viewXYFlag == true) 
	    hit->SetRecoHit(false);
	  else hit->SetRecoHit(true);
	  //std::cout<<"Z fiber in event "<<subSpill<<" with X and Y "<<aveXSave<<" "<<aveYSave<<std::endl;
	}
        if(dim == 1){
	  hit->SetView(1);
          hit->SetX(aveXSave);
          hit->SetZ(aveZSave);
          hit->SetY(-1);
          if (viewXZFlag == true)
            hit->SetRecoHit(false);
          else hit->SetRecoHit(true);	  
        }
        if(dim == 0){
	  hit->SetView(2);
          hit->SetZ(aveZSave);
          hit->SetY(aveYSave);
          hit->SetX(-1);
          if (viewYZFlag == true)
            hit->SetRecoHit(false);
          else hit->SetRecoHit(true);	  
        }

        hit->SetTrueXTalk(false);
	hit->SetTrueX(aveXSave);
	hit->SetTrueY(aveYSave);
	hit->SetTrueZ(aveZSave);
	hit->SetMC(true);
	hit->SetPDG(neutronHitPDG[locSave]);
	hit->SetParentID(neutronParentId[locSave]);
	hit->SetTrackID(neutronParentPDG[locSave]);
	hit->SetEdep(edep);


      //cross talk cubes
      for (Int_t iXs = 0; iXs<2; iXs++){

        Hit* hit = event->AddHit();

        hit->SetPE(hitPE[dim]*(crossFraction));
        hit->SetHG_pe(hitPE[dim]*(crossFraction));
        hit->SetLG_pe(hitPE[dim]*(crossFraction));
        hit->SetToT_pe(hitPE[dim]*(crossFraction));
        hit->SetHG_ADC(adc[dim]);
        hit->SetLG_ADC(loadc[dim]);
        hit->SetRE(0);
        hit->SetFE(0);
        hit->SetToT(0);

        tempT = (hitT[dim]-t0[0])/tickSize;
        hit->SetSpillTag(spill[0]);

	// set spilltime as true time 
        hit->SetSpillTime(hitT[dim]-t0[0]);

        hit->SetSpillTrailTime(0);
        hit->SetTfromSpill(neutron_time[0]);
        hit->SetFEB(0);
        hit->SetCh(0);
        hit->SetGTrigTag((int)hitT[dim]/tickSize);
        hit->SetGTrigTime((int)hitT[dim]/tickSize);

        if(dim == 2){
	  hit->SetView(0);
          hit->SetX(aveXSave-1+iXs*2);
          hit->SetY(aveYSave);
          // these cross talk numbers are extracted by eye from : https://indico.cern.ch/event/930603/contributions/3928119/attachments/2068887/3472703/LANL_discussion.pdf	  
          double xtalkDelay = gRandom->Gaus((2.5*7)/tickSize, (3*2.5)/tickSize);
          hit->SetDt(tempT+xtalkDelay);	  
          hit->SetZ(-1);

          if (viewXYFlag == true)
            hit->SetRecoHit(false);
          else hit->SetRecoHit(true);

    	  // y = a + b * x ; x = PE , y = tot tick, a = 8.8, b = 0.2; above 80 pe, all 25 ticks and below 6 pe, all 10 ticks
          durationZ = 8.8 + 0.2 * hitPE[2];

          if(aveXSave+24<100 && aveYSave+24<100 && aveZSave+50< 101 && veryTempT< 10000){
            cubeFlagXY[aveXSave+24][aveYSave+24] =1 ;
            timeFlag[(int)(veryTempT+xtalkDelay)] = 1;

            for(int lll=0;lll<durationZ; lll++)
              timeFlagXY[(int)(veryTempT+xtalkDelay+lll)] = 1;
          }
	  
	}
        if(dim == 1){
	  hit->SetView(1);	
          hit->SetX(aveXSave);
          hit->SetZ(aveZSave-1+iXs*2);
          double xtalkDelay = gRandom->Gaus(10./tickSize, 7.5/tickSize);
          hit->SetDt(tempT+xtalkDelay);	  
          hit->SetY(-1);

          if (viewXZFlag == true)
            hit->SetRecoHit(false);
          else hit->SetRecoHit(true);

    	  // y = a + b * x ; x = PE , y = tot tick, a = 8.8, b = 0.2; above 80 pe, all 25 ticks and below 6 pe, all 10 ticks
          durationY = 8.8 + 0.2 * hitPE[1];

          if(aveXSave+24<100 && aveYSave+24<100 && aveZSave+50< 101 && veryTempT< 10000){
            cubeFlagXZ[aveXSave+24][aveZSave+50] =1 ;
            timeFlag[(int)(veryTempT+xtalkDelay)] = 1;

            for(int lll=0;lll<durationY; lll++)
              timeFlagXZ[(int)(veryTempT+xtalkDelay+lll)] = 1;
          }
	  
        }
        if(dim == 0){
	  hit->SetView(2);	
          hit->SetZ(aveZSave);
          hit->SetY(aveYSave-1+iXs*2);
          double xtalkDelay = gRandom->Gaus((2.5*7)/tickSize, (3*2.5)/tickSize);
          hit->SetDt(tempT+xtalkDelay);
	  hit->SetX(-1);

          if (viewYZFlag == true)
            hit->SetRecoHit(false);
          else hit->SetRecoHit(true);

          // y = a + b * x ; x = PE , y = tot tick, a = 8.8, b = 0.2; above 80 pe, all 25 ticks and below 6 pe, all 10 ticks
          durationX = 8.8 + 0.2 * hitPE[0];

          if(aveXSave+24<100 && aveYSave+24<100 && aveZSave+50< 101 && veryTempT< 10000){
            cubeFlagYZ[aveYSave+24][aveZSave+50] =1 ;
            timeFlag[(int)(veryTempT+xtalkDelay)] = 1;

            for(int lll=0;lll<durationX; lll++)
              timeFlagYZ[(int)(veryTempT+xtalkDelay+lll)] = 1;
          }
	  
        }

        hit->SetTrueXTalk(true);
        hit->SetTrueX(aveXSave-1+iXs*2);
        hit->SetTrueY(aveYSave-1+iXs*2);
        hit->SetTrueZ(aveZSave-1+iXs*2);
        hit->SetMC(true);
        hit->SetPDG(-1);
        hit->SetParentID(-1);
        hit->SetTrackID(-1);
        hit->SetEdep(0);

      }
      
      }
      edep = neutronHitE[hitLoop];
      timeSave = 1e10;
      //std::cout<<"added a hit "<<std::endl;
      
    }
    //aveXSave = neutronHitX[hitLoop]/cubeSize+xx;
    //aveYSave = neutronHitY[hitLoop]/cubeSize+yy;
    //aveZSave = neutronHitZ[hitLoop]/cubeSize+zz;    
    //aveXSave -= xx;
    //aveYSave -= yy;
    //aveZSave -= zz;
    aveXSave = aveX;
    aveYSave = aveY;
    aveZSave = aveZ;
    locSave = hitLoop; 
    //cout<<"check again x y z "<<aveXSave<<" "<<aveYSave<<" "<<aveZSave<<endl;
    }

    if(1){

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
    }

    event->SetTrueVoxelNumber(locCounter);
    for(int ivt=0;ivt<1000;ivt++){
	//std::cout<<":::::::::::::::::  "<<trueLocationX[ivt]<<" "<<trueLocationY[ivt]<<" "<<trueLocationZ[ivt]<<std::endl;
	event->SetTrueLocation(ivt, 0, (double)trueLocationX[ivt]);
	event->SetTrueLocation(ivt, 1, (double)trueLocationY[ivt]);
	event->SetTrueLocation(ivt, 2, (double)trueLocationZ[ivt]);
	event->SetTrueLocation(ivt, 3, (double)trueLocationT[ivt]);
    }

    event->SetColmt(0, vtx[0]);
    event->SetColmt(1, vtx[1]);
    for(int ivt=0;ivt<10;ivt++)
      event->SetTrueVtx(ivt,true_vtx[ivt]);
    cout<<"initial momentum along x y z : "<<event->GetdEdz(0)<<" "<<event->GetdEdz(1)<<" "<<event->GetdEdz(2)<<endl;

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
  }

  wfile.cd();
  AllEvents.Write("",TObject::kOverwrite);
  //AllEvents.Delete();
  wfile.Close();
  //t.Delete();

  //return 0;
}
