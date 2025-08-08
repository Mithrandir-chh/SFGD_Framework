#define THIS_NAME cubeLightYield
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS
#define _CRTDBG_MAP_ALLOC

#include "../src/tools/global_header.hh"
#include "TSystem.h"
#include "TMacro.h"
#include <sstream>
#include <iterator>
#include <algorithm>
#include "TImage.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "Rtypes.h"
#include <TPaveLabel.h>

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool    batch          = false;
bool    IsMC           = true;
bool    IsReversed     = false;
bool    IsTwenty       = false;
int     maxEvents = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "~/Desktop/templateAna_Out.root";
TString fileIn         = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
int IsPrototype        = true;
int SFGD_X             = 204;
int SFGD_Y             = 56;
int SFGD_Z             = 188;
bool useClean          = false;
bool useNN             = false;
int dTimeIni           = -100;
int dTimeFin           = -100;
size_t evtIni             = 0;
size_t evtFin             = 0;

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

void cubeLightYield() {


  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent,
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  cout<<"total events = "<<nEvents<<endl;
  cout<<"displaying "<<evtFin-evtIni<<" events"<<endl;

  // Create Histograms________________________________________________________________________________

  TH1D *CubeLY_X = new TH1D("CubeLY_X","Uncorrected Cube Light Yields along X;Light Yield [p.e.];Number of entries", 120,0,120);
  TH1D *CubeLY_Xlayer[24];
  TH1D *CubeLY_Xlayercorr[24];
  TH1D *CubeLY_Xdlayer[24];
  TH1D *CubeLY_Xdlayercorr[24];
  TH1D *CubeLY_Xfiber[48][8];
  TH1D *CubeLY_Ylayer[24];
  TH1D *CubeLY_Ylayercorr[24];
  TH1D *CubeLY_Ydlayer[24];
  TH1D *CubeLY_Ydlayercorr[24];
  TH1D *CubeLY_Yfiber[48][24];
  // TH1D *CubeLY_Xdfiber[24][48][8];
  for (Int_t i=0; i<24; i++){
    CubeLY_Xlayer[i] = new TH1D(("CubeLY_X"+to_string(i)+"layer").c_str(),("Uncorrected Cube Light Yields along X, X="+to_string(i)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Xlayercorr[i] = new TH1D(("CubeLY_X"+to_string(i)+"layercorr").c_str(),("Corrected Cube Light Yields along X, X="+to_string(i)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Xdlayer[i] = new TH1D(("CubeLY_Xd"+to_string(i+2)+"layer").c_str(),("Uncorrected Cube Light Yields along X, d="+to_string(i+2)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Xdlayercorr[i] = new TH1D(("CubeLY_Xd"+to_string(i+2)+"layercorr").c_str(),("Corrected Cube Light Yields along X, d="+to_string(i+2)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Ylayer[i] = new TH1D(("CubeLY_Y"+to_string(i)+"layer").c_str(),("Uncorrected Cube Light Yields along X, X="+to_string(i)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Ylayercorr[i] = new TH1D(("CubeLY_Y"+to_string(i)+"layercorr").c_str(),("Corrected Cube Light Yields along Y, Y="+to_string(i)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Ydlayer[i] = new TH1D(("CubeLY_Yd"+to_string(i+2)+"layer").c_str(),("Uncorrected Cube Light Yields along Y, d="+to_string(i+2)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    CubeLY_Ydlayercorr[i] = new TH1D(("CubeLY_Yd"+to_string(i+2)+"layercorr").c_str(),("Corrected Cube Light Yields along Y, d="+to_string(i+2)+";Light Yield [p.e.];Number of entries").c_str(), 120,0,120);
    // for (Int_t iz=0; iz<48; iz++){
      // for (Int_t iy=0; iy<8; iy++){
        // CubeLY_Xdfiber[i][iz][iy] = new TH1D(("CubeLY_Xdfiber_"+to_string(i+2)+"_"+to_string(iz)+"_"+to_string(iy)).c_str(),("Average cube LY values along x fiber "+to_string(iz)+"_"+to_string(iy)+" at d="+to_string(i+2)+"; Light Yield [p.e.]; Number of Entries").c_str(), 120,0,120); 
      // }
    // }
  }
  for (Int_t iz=0; iz<48; iz++){
    for (Int_t iy=0; iy<8; iy++){
      CubeLY_Xfiber[iz][iy] = new TH1D(("CubeLY_Xfiber_"+to_string(iz)+"_"+to_string(iy)).c_str(),("x fiber "+to_string(iz)+"_"+to_string(iy)+"; d [cm]; Light Yield [p.e.]").c_str(), 24,2,26); 
    }
    for (Int_t ix=0; ix<24; ix++){
      CubeLY_Yfiber[iz][ix] = new TH1D(("CubeLY_Yfiber_"+to_string(iz)+"_"+to_string(ix)).c_str(),("y fiber "+to_string(iz)+"_"+to_string(ix)+"; d [cm]; Light Yield [p.e.]").c_str(), 8,2,10); 
    }
  }
  TH1D *CubeLY_Y = new TH1D("CubeLY_Y","Uncorrected Cube Light Yields along Y;Light Yield [p.e.];Number of entries", 120,0,120);
  TH1D *CubeLY_Xcorr = new TH1D("CubeLY_Xcorr","Corrected Cube Light Yields along X;Light Yield [p.e.];Number of entries", 120,0,120);
  TH1D *CubeLY_Ycorr = new TH1D("CubeLY_Ycorr","Corrected Cube Light Yields along Y;Light Yield [p.e.];Number of entries", 120,0,120);
  TH1D *Z_dist = new TH1D("Z_dist","Average LY per Z layer; Z [cm]; Average LY [p.e.]", 48,0,48);
  TH1D *Z_distCorr = new TH1D("Z_distCorr","Average corrected LY per Z layer; Z [cm]; Average LY [p.e.]", 48,0,48);
  // TH2D *Q_XZ = new TH2D("Q_XZ","Average charge per MPPC in XZ plane; X [cm]; Z [cm]; Average Charge [p.e.]", 24,0,24, 48,0,48);
  // TH2D *Qcorr_XZ = new TH2D("Qcorr_XZ","Average corrected charge per MPPC in XZ plane; X [cm]; Z [cm]; Average Charge [p.e.]", 24,0,24, 48,0,48);
  // TH2D *Q_ZY = new TH2D("Q_XZ","Average charge per MPPC in ZY plane; Z [cm]; Y [cm]; Average Charge [p.e.]", 48,0,48, 8,0,8);
  // TH2D *Qcorr_ZY = new TH2D("Qcorr_XZ","Average corrected charge per MPPC in ZY plane; Z [cm]; Y [cm]; Average Charge [p.e.]", 48,0,48, 8,0,8);
  Double_t cubeX_Q[24][8][48] = {0};
  Double_t cubeX_Qcorr[24][8][48] = {0};
  Int_t cubeX_N[24][8][48] = {0};
  Double_t cubeXd_Q[24][8][48] = {0};
  Double_t cubeXd_Qcorr[24][8][48] = {0};
  Double_t cubeYd_Q[24][8][48] = {0};
  Double_t cubeYd_Qcorr[24][8][48] = {0};
  Int_t cubeXd_N[24][8][48] = {0};
  Int_t cubeYd_N[24][8][48] = {0};
  Double_t cubeY_Q[24][8][48] = {0};
  Double_t cubeY_Qcorr[24][8][48] = {0};
  Int_t cubeY_N[24][8][48] = {0};
  Double_t zQ[48] = {0};
  Double_t zQcorr[48] = {0};
  Int_t zN[48] = {0};

  // Double_t Ls = 0.0854; //m
  Double_t Ls = 0.05883; //m
  // Double_t Ll = 4.13; //m
  Double_t Ll = 4.; //m
  // Double_t alpha = 0.1716;
  Double_t alpha = 0.1453;
  // Double_t mppcOffset = 0.02; //m
  Double_t mppcOffset = 0.023; //m

  // Filling Histograms_______________________________________________________________________________

  Int_t prevEventTime = 0;
  Int_t passed = 0;
  //loop over events
  for (size_t iev=evtIni; iev<evtFin; iev++){
    if (iev == maxEvents-1) break;
    cout<<"_Processing events..."<<(int)(iev*100/(Double_t)(evtFin-evtIni))<<"% done         \r"<<flush;
    
    std::vector<ND280SFGDHit*> mppc_hits;
    mppc_hits = getEventMPPChits(iev);

    //skip FEB12 triggers that are too close in time to the previous trigger
    if (abs(unpackEvent->GetFEB12hitTfromSpill() - prevEventTime) < 250){
      for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
      continue;
    }

    //Skip event if it doesn't look like a muon going straight through the detector along the z axis
    if (mppc_hits.size() < 20 || unpackEvent->GetMaxCharge() < 50 || unpackEvent->GetEmptyZ() > 2 || unpackEvent->GetXStdDev() > 0.5 || unpackEvent->GetYStdDev() > 0.5){ 
      for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
      continue;
    }
    passed++;
    
    //check all main hits have the same x, y coords
    // Double_t hitX;
    // Double_t hitY;
    // Double_t hitX_p;
    // Double_t hitY_p;
    // bool first = true;
    // bool skipEvent = false;
    // for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      // ND280SFGDHit *hit = mppc_hits[ihit];
      // if (hit->GetPE() < 20 || hit->GetView() == 0 || hit->GetDt() < -130 || hit->GetDt() > -107) continue;
      // hitX = hit->GetX();
      // hitY = hit->GetY();
      // if (! first){
        // if (hitX != -1){
          // if (hitX != hitX_p){
            // skipEvent = true;
            // break;
          // }
        // }
        // if (hitY != -1){
          // if (hitY != hitY_p){
            // skipEvent = true;
            // break;
          // }
        // }
      // }
      // first = false;
      // hitX_p = hitX;
      // hitY_p = hitY;
    // }
    // if (skipEvent) continue;

    //loop over hits
    for (size_t ihit=0; ihit<mppc_hits.size(); ihit++){
      Double_t travelDist;
      ND280SFGDHit *hit = mppc_hits[ihit];
      //only keep type 1 mppc hits
      if (hit->GetFEB() == 8 || hit->GetFEB() == 11 || hit->GetFEB() == 3 || hit->GetFEB() == 4 || hit->GetFEB() == 9 || hit->GetFEB() == 10)
        continue;

      // if (hit->GetDt() < -130 || hit->GetDt() > -107) //time window where we expect muons
        // continue;

      if (hit->GetPE() < 20) continue; //cut out crosstalk with charge cut
      
      if (hit->GetView() == 1){ //XZ view
        //find y coordinate
        Double_t y;
        Double_t pairCharge = 0;
        for (size_t pHit=0; pHit<mppc_hits.size(); pHit++){
          if (mppc_hits[pHit]->GetView() == 2 && hit->GetZ() == mppc_hits[pHit]->GetZ() && mppc_hits[pHit]->GetPE() > pairCharge && mppc_hits[pHit]->GetPE() > 20){
            pairCharge = mppc_hits[pHit]->GetPE();
            y = mppc_hits[pHit]->GetY();
            Int_t hitFEB = (int)hit->GetFEB();
            if (hitFEB == 3 || hitFEB == 4 || hitFEB == 8 || hitFEB == 18 || hitFEB == 19 || hitFEB == 20 ) // top channels
              travelDist = mppcOffset+(8.-y-0.5)*0.01; // meters
            else // bottom channels
              travelDist = mppcOffset+(y+0.5)*0.01; // meters
          }
        }

        cubeY_Q[(int)hit->GetX()][(int)y][(int)hit->GetZ()] += hit->GetPE();
        cubeY_Qcorr[(int)hit->GetX()][(int)y][(int)hit->GetZ()] += hit->GetPE()/(alpha*exp(-travelDist/Ls)+(1-alpha)*exp(-travelDist/Ll));
        cubeY_N[(int)hit->GetX()][(int)y][(int)hit->GetZ()]++;
        cubeYd_Q[(int)hit->GetX()][(int)((travelDist-mppcOffset)*100)][(int)hit->GetZ()] += hit->GetPE();
        cubeYd_Qcorr[(int)hit->GetX()][(int)((travelDist-mppcOffset)*100)][(int)hit->GetZ()] += hit->GetPE()/(alpha*exp(-travelDist/Ls)+(1-alpha)*exp(-travelDist/Ll));
        cubeYd_N[(int)hit->GetX()][(int)((travelDist-mppcOffset)*100)][(int)hit->GetZ()]++;
        zQ[(int)hit->GetZ()] += hit->GetPE();
        zQcorr[(int)hit->GetZ()] += hit->GetPE()/(alpha*exp(-travelDist/Ls)+(1-alpha)*exp(-travelDist/Ll));
        zN[(int)hit->GetZ()]++;
        // CubeLY_Y->Fill(hit->GetPE());
        // CubeLY_Ycorr->Fill(hit->GetPE()/(alpha*exp();
      }
      else if (hit->GetView() == 2){ //ZY view
        //find x coordinate
        Double_t x;
        Double_t pairCharge = 0;
        for (size_t pHit=0; pHit<mppc_hits.size(); pHit++){
          if (mppc_hits[pHit]->GetView() == 1 && hit->GetZ() == mppc_hits[pHit]->GetZ() && mppc_hits[pHit]->GetPE() > pairCharge && mppc_hits[pHit]->GetPE() > 20){
            pairCharge = mppc_hits[pHit]->GetPE();
            x = mppc_hits[pHit]->GetX();
            Int_t hitFEB = (int)hit->GetFEB();
            if (hitFEB == 1 || hitFEB == 24){ // right channels
              travelDist = mppcOffset+(24.-x-0.5)*0.01; // meters
            }
            else{ // left channels
              travelDist = mppcOffset+(x+0.5)*0.01; // meters
            }
          }
        }
        // cout << x << " " << travelDist*100 << endl;

        cubeX_Q[(int)x][(int)hit->GetY()][(int)hit->GetZ()] += hit->GetPE();
        cubeX_Qcorr[(int)x][(int)hit->GetY()][(int)hit->GetZ()] += hit->GetPE()/(alpha*exp(-travelDist/Ls)+(1-alpha)*exp(-travelDist/Ll));
        cubeX_N[(int)x][(int)hit->GetY()][(int)hit->GetZ()]++;
        cubeXd_Q[(int)((travelDist-mppcOffset)*100)][(int)hit->GetY()][(int)hit->GetZ()] += hit->GetPE();
        cubeXd_Qcorr[(int)((travelDist-mppcOffset)*100)][(int)hit->GetY()][(int)hit->GetZ()] += hit->GetPE()/(alpha*exp(-travelDist/Ls)+(1-alpha)*exp(-travelDist/Ll));
        cubeXd_N[(int)((travelDist-mppcOffset)*100)][(int)hit->GetY()][(int)hit->GetZ()]++;
        zQ[(int)hit->GetZ()] += hit->GetPE();
        zQcorr[(int)hit->GetZ()] += hit->GetPE()/(alpha*exp(-travelDist/Ls)+(1-alpha)*exp(-travelDist/Ll));
        zN[(int)hit->GetZ()]++;
        // CubeLY_X->Fill(hit->GetPE());
        // CubeLY_Xcorr->Fill(hit->GetPE());
      }
      prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
      //cout << hit->GetSpillTime()*2.5 << endl;
    }  

    for (size_t i=0; i<mppc_hits.size(); i++) delete mppc_hits[i];
  }
  cout<<"_Processing events...100% done"<<endl;
  cout << passed << " events passed cuts" << endl;

  //Draw histograms___________________________________________________________________________________

  Int_t minEntries = 450;
  for (int iz=0; iz<48; iz++){
    zQ[iz] = zQ[iz]/zN[iz];
    zQcorr[iz] = zQcorr[iz]/zN[iz];
    Z_dist->Fill(iz,zQ[iz]);
    Z_distCorr->Fill(iz,zQcorr[iz]);
    for (int ix=0; ix<24; ix++){
      for (int iy=0; iy<8; iy++){
        cubeY_Q[ix][iy][iz] = cubeY_Q[ix][iy][iz]/cubeY_N[ix][iy][iz];
        cubeY_Qcorr[ix][iy][iz] = cubeY_Qcorr[ix][iy][iz]/cubeY_N[ix][iy][iz];
        cubeYd_Q[ix][iy][iz] = cubeYd_Q[ix][iy][iz]/cubeYd_N[ix][iy][iz];
        cubeYd_Qcorr[ix][iy][iz] = cubeYd_Qcorr[ix][iy][iz]/cubeYd_N[ix][iy][iz];
        if (cubeY_N[ix][iy][iz] > minEntries){
          CubeLY_Ylayer[iy]->Fill(cubeY_Q[ix][iy][iz]);
          CubeLY_Ylayercorr[iy]->Fill(cubeY_Qcorr[ix][iy][iz]);
          CubeLY_Y->Fill(cubeY_Q[ix][iy][iz]);
          CubeLY_Ycorr->Fill(cubeY_Qcorr[ix][iy][iz]);
        }
        if (cubeYd_N[ix][iy][iz] > minEntries){
          CubeLY_Ydlayer[iy]->Fill(cubeYd_Q[ix][iy][iz]);
          CubeLY_Ydlayercorr[iy]->Fill(cubeYd_Qcorr[ix][iy][iz]);
          CubeLY_Yfiber[ix][iz]->Fill(ix+2.5, cubeYd_Q[ix][iy][iz]);
        }
        cubeX_Q[ix][iy][iz] = cubeX_Q[ix][iy][iz]/cubeX_N[ix][iy][iz];
        cubeX_Qcorr[ix][iy][iz] = cubeX_Qcorr[ix][iy][iz]/cubeX_N[ix][iy][iz];
        cubeXd_Q[ix][iy][iz] = cubeXd_Q[ix][iy][iz]/cubeXd_N[ix][iy][iz];
        cubeXd_Qcorr[ix][iy][iz] = cubeXd_Qcorr[ix][iy][iz]/cubeXd_N[ix][iy][iz];
        if (cubeX_N[ix][iy][iz] > minEntries){
          CubeLY_Xlayer[ix]->Fill(cubeX_Q[ix][iy][iz]);
          CubeLY_Xlayercorr[ix]->Fill(cubeX_Qcorr[ix][iy][iz]);
          CubeLY_X->Fill(cubeX_Q[ix][iy][iz]);
          CubeLY_Xcorr->Fill(cubeX_Qcorr[ix][iy][iz]);
        }
        if (cubeXd_N[ix][iy][iz] > minEntries){
          CubeLY_Xdlayer[ix]->Fill(cubeXd_Q[ix][iy][iz]);
          CubeLY_Xdlayercorr[ix]->Fill(cubeXd_Qcorr[ix][iy][iz]);
          CubeLY_Xfiber[iz][iy]->Fill(ix+2.5, cubeXd_Q[ix][iy][iz]);
        }
      }
    }
  }

  TString localStyleName = "T2K";
  // -- WhichStyle --
  // 1 = presentation large fonts
  // 2 = presentation small fonts
  // 3 = publication/paper
  Int_t localWhichStyle = 3;
  TStyle* t2kstyle = SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->ForceStyle(t2kstyle);

  gStyle->SetOptFit(111);
  gStyle->SetOptStat(1110);
  TCanvas *cc = new TCanvas("cc","cc", 0, 0, 800, 630);
  TF1 *f1 = new TF1("f1", "gaus", 0, 120);
  cc->cd();

  TF1 *fit;
  if (CubeLY_X->GetEntries() > 10){
    CubeLY_X->Scale(1/CubeLY_X->Integral("width"));
    CubeLY_X->Fit("f1","QSR");
    fit = (TF1*)CubeLY_X->GetListOfFunctions()->FindObject("f1");
    // fit->SetLineWidth(3);
    // CubeLY_X->SetFillColorAlpha(kGreen-8,0.559);
    // CubeLY_X->SetLineWidth(2);
    // CubeLY_Xdlayer[7]->Draw("same HIST");
    // CubeLY_Xdlayer[7]->SetFillColorAlpha(kRed,0.5);
    // CubeLY_Xdlayer[10]->Draw("same HIST");
    // CubeLY_Xdlayer[10]->SetFillColorAlpha(kGreen,0.5);
    // CubeLY_Xdlayer[14]->Draw("same HIST");
    // CubeLY_Xdlayer[14]->SetFillColorAlpha(kYellow,0.5);
    TLegend *leg = new TLegend(0.1,0.7,0.3,0.9);
    leg->AddEntry(CubeLY_X, "All cubes", "f");
    leg->AddEntry(CubeLY_Xdlayer[10], "Cubes where d = 9.5 cm", "f");
    // leg->AddEntry(CubeLY_Xdlayer[10], "Cubes where d = 9 cm", "f");
    leg->AddEntry(CubeLY_Xdlayer[14], "Cubes where d = 16.5 cm", "f");
    leg->SetBorderSize(0);
    // leg->Draw();
    cc->Modified();
    cc->Update();
    cc->SaveAs("CubeLY_X.C");
  }
  else cout << CubeLY_X->GetEntries() << endl;
  // cc->Clear();
  TH1D *meanxLYd = new TH1D("meanxLYd","Mean x LY for each d value; d [cm]; Mean LY [p.e.]", 24, 0, 24);
  TH1D *meanxLYdcorr = new TH1D("meanxLYdcorr","Mean corrected x LY for each d value; d [cm]; Mean LY [p.e.]", 24, 0, 24);
  TH1D *meanxLY = new TH1D("meanLxY","Mean LY for each x value; x [cm]; Mean LY [p.e.]", 24, 0, 24);
  TH1D *meanyLYd = new TH1D("meanyLYd","Mean y LY for each d value; d [cm]; Mean LY [p.e.]", 8, 2, 10);
  TH1D *meanyLYdcorr = new TH1D("meanyLYdcorr","Mean corrected y LY for each d value; d [cm]; Mean LY [p.e.]", 8, 2, 10);
  TH1D *meanyLY = new TH1D("meanyLY","Mean LY for each y value; y [cm]; Mean LY [p.e.]", 8, 2, 10);
  for (int i=0; i<24; i++){
    if (CubeLY_Xlayer[i]->GetEntries() > 10){
      TFitResultPtr r = CubeLY_Xlayer[i]->Fit("f1","QSR");
      meanxLY->Fill(i+2.5,r->Parameter(1));
      meanxLY->SetBinError(i+3,r->ParError(1));
      fit = (TF1*)CubeLY_Xlayer[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Xlayer[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Xlayer[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Xlayer"+to_string(i)+".C").c_str());
      CubeLY_Xlayercorr[i]->Fit("f1","QSR");
      fit = (TF1*)CubeLY_Xlayercorr[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Xlayercorr[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Xlayercorr[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Xlayercorr"+to_string(i)+".C").c_str());
    }
    if (CubeLY_Xdlayer[i]->GetEntries() > 10){
      TFitResultPtr r = CubeLY_Xdlayer[i]->Fit("f1","QSR");
      meanxLYd->Fill(i+2.5,r->Parameter(1));
      meanxLYd->SetBinError(i+3,r->ParError(1));
      fit = (TF1*)CubeLY_Xdlayer[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Xdlayer[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Xdlayer[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Xdlayer"+to_string(i)+".C").c_str());
      r = CubeLY_Xdlayercorr[i]->Fit("f1", "QSR");
      meanxLYdcorr->Fill(i+2.5,r->Parameter(1));
      meanxLYdcorr->SetBinError(i+3,r->ParError(1));
      fit = (TF1*)CubeLY_Xdlayercorr[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Xdlayercorr[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Xdlayercorr[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Xdlayercorr"+to_string(i)+".C").c_str());
    }
  }
  for (int i=0; i<8; i++){
    if (CubeLY_Ylayer[i]->GetEntries() > 10){
      TFitResultPtr r = CubeLY_Ylayer[i]->Fit("f1","QSR");
      meanyLY->Fill(i+2.5,r->Parameter(1));
      meanyLY->SetBinError(i+3,r->ParError(1));
      fit = (TF1*)CubeLY_Ylayer[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Ylayer[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Ylayer[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Ylayer"+to_string(i)+".C").c_str());
      CubeLY_Ylayercorr[i]->Fit("f1","QSR");
      fit = (TF1*)CubeLY_Ylayercorr[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Ylayercorr[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Ylayercorr[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Ylayercorr"+to_string(i)+".C").c_str());
    }
    if (CubeLY_Ydlayer[i]->GetEntries() > 10){
      TFitResultPtr r = CubeLY_Ydlayer[i]->Fit("f1","QSR");
      meanyLYd->Fill(i+2.5,r->Parameter(1));
      meanyLYd->SetBinError(i+3,r->ParError(1));
      fit = (TF1*)CubeLY_Ydlayer[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Ydlayer[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Ydlayer[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Ydlayer"+to_string(i)+".C").c_str());
      r = CubeLY_Ydlayercorr[i]->Fit("f1", "QSR");
      meanyLYdcorr->Fill(i+2.5,r->Parameter(1));
      meanyLYdcorr->SetBinError(i+3,r->ParError(1));
      fit = (TF1*)CubeLY_Ydlayercorr[i]->GetListOfFunctions()->FindObject("f1");
      fit->SetLineWidth(3);
      CubeLY_Ydlayercorr[i]->SetFillColorAlpha(kGreen-8,0.559);
      CubeLY_Ydlayercorr[i]->SetLineWidth(2);
      cc->Update();
      cc->SaveAs(("CubeLY_Ydlayercorr"+to_string(i)+".C").c_str());
    }
    // for (int iy=0; iy<8; iy++){
      // for (int iz=0; iz<48; iz++){
        // if (CubeLY_Xdfiber[i][iz][iy]->GetEntries() < 10) continue;
        // cout << "Fitting" << endl;
        // TFitResultPtr r = CubeLY_Xdfiber[i][iz][iy]->Fit("f1","QSR");
        // CubeLY_Xfiber[iz][iy]->Fill(i+2.5, r->Parameter(1));
      // }
    // }
  }
  meanxLY->Draw();
  cc->Update();
  cc->SaveAs("meanxLY.C");
  meanxLYd->Draw();
  cc->Update();
  cc->SaveAs("meanxLYd.C");
  meanyLY->Draw();
  cc->Update();
  cc->SaveAs("meanyLY.C");
  meanyLYd->Draw();
  cc->Update();
  cc->SaveAs("meanyLYd.C");
  TCanvas *cf = new TCanvas("cf","cf",0,0,1900,1030);
  cf->cd();
  cf->Divide(8,5);
  int canv = 40;
  gStyle->SetOptStat(0);
  for (int iy=0; iy<8; iy++){
    for (int iz=47; iz>=0; iz--){
      if (CubeLY_Xfiber[iz][iy]->GetEntries() < 3) continue;
      if (iz%6 != 0) continue;
      if (iy == 0 || iy == 1 || iy == 2 || iy == 6 || iy == 7){
        cf->cd(canv);
        canv--;
        CubeLY_Xfiber[iz][iy]->GetYaxis()->SetRangeUser(30,80);
        CubeLY_Xfiber[iz][iy]->Draw("hist p");
      }
      // cc->Update();
      // cc->SaveAs(("xfiber_"+to_string(iz)+"_"+to_string(iy)+".C").c_str());
    }
  }
  cf->SaveAs("meanxLYfibers.C");
  canv = 40;
  for (int ix=0; ix<24; ix++){
    for (int iz=47; iz>=0; iz--){
      if (CubeLY_Yfiber[iz][ix]->GetEntries() < 3) continue;
      if (iz%6 != 0) continue;
      if (ix == 8 || ix == 10 || ix == 12 || ix == 14 || ix == 16){
        cf->cd(canv);
        canv--;
        CubeLY_Yfiber[iz][ix]->GetYaxis()->SetRangeUser(30,80);
        CubeLY_Yfiber[iz][ix]->Draw("hist p");
      }
      // cc->Update();
      // cc->SaveAs(("xfiber_"+to_string(iz)+"_"+to_string(iy)+".C").c_str());
    }
  }
  cf->SaveAs("meanyLYfibers.C");
  gStyle->SetOptStat(1110);
  cc->cd();
  if (CubeLY_Y->GetEntries() > 10){
    CubeLY_Y->Scale(1/CubeLY_Y->Integral("width"));
    CubeLY_Y->Fit("f1","QSR");
    fit = (TF1*)CubeLY_Y->GetListOfFunctions()->FindObject("f1");
    // fit->SetLineWidth(3);
    // CubeLY_Y->SetFillColorAlpha(kGreen-8,0.559);
    // CubeLY_Y->SetLineWidth(2);
    cc->Update();
    cc->SaveAs("CubeLY_Y.C");
    CubeLY_Xcorr->Scale(1/CubeLY_Xcorr->Integral("width"));
    CubeLY_Xcorr->Fit("f1","QSR");
    fit = (TF1*)CubeLY_Xcorr->GetListOfFunctions()->FindObject("f1");
    fit->SetLineWidth(3);
    // CubeLY_Xdlayercorr[7]->Draw("same HIST");
    // CubeLY_Xdlayercorr[7]->SetFillColorAlpha(kRed,0.5);
    // CubeLY_Xdlayercorr[13]->Draw("same HIST");
    // CubeLY_Xdlayercorr[13]->SetFillColorAlpha(kYellow,0.5);
    // CubeLY_Xcorr->SetFillColorAlpha(kGreen-8,0.559);
    CubeLY_Xcorr->SetLineWidth(2);
    cc->Update();
    cc->SaveAs("CubeLY_Xcorr.C");
    CubeLY_Ycorr->Scale(1/CubeLY_Ycorr->Integral("width"));
    CubeLY_Ycorr->Fit("f1","QSR");
    fit = (TF1*)CubeLY_Ycorr->GetListOfFunctions()->FindObject("f1");
    fit->SetLineWidth(3);
    // CubeLY_Ycorr->SetFillColorAlpha(kGreen-8,0.559);
    CubeLY_Ycorr->SetLineWidth(2);
    cc->Update();
    cc->SaveAs("CubeLY_Ycorr.C");
  }
  Z_dist->Draw("hist");
  cc->Update();
  cc->SaveAs("Zdist.C");
  Z_distCorr->Draw("hist");
  cc->Update();
  cc->SaveAs("ZdistCorr.C");

  FileOutput->Close();
  FileInput->Close();
  handleEndOfExecution();
}
