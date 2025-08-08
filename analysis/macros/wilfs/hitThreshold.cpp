#define THIS_NAME hitThreshold
#define NOINTERACTIVE_OUTPUT
#define OVERRIDE_OPTIONS

#include "../src/tools/global_header.hh"
#include "TSystem.h"
#include "TMacro.h"
#include <sstream>
#include <iterator>
#include <algorithm>
#include "TImage.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "Rtypes.h"
#include <TPaveLabel.h>
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLegend.h"
#include "TImage.h"
#include "TLeaf.h"
#include "TColor.h"
#include "TPaletteAxis.h"
#include <vector>

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool batch = false;
bool IsMC = true;
bool IsReversed = false;
bool IsTwenty = false;
int maxEvents = std::numeric_limits<int>::max();
;
int maxSelEvents = std::numeric_limits<int>::max();
;
int selEvents = 0;
float start_time = clock();
bool RM_CROSSTALK = false;
bool SHOW_TRUE = true;
bool SHOW_RECO = true;
TString fileOut = "~/Desktop/templateAna_Out.root";
TString fileIn = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
bool IsPrototype = true;
int SFGD_X = 204;
int SFGD_Y = 56;
int SFGD_Z = 188;
bool useClean = false;
bool useNN = false;
int dTimeIni = -100;
int dTimeFin = -100;
int evtIni = 0;
int evtFin = 0;

TFile *FileInput;
TFile *FileOutput;
TTree *dataOut;
TTree *data;
ND280SFGDEvent *inputEvent;
ND280SFGDEvent *recoEvent;
TBranch *recoBranch;
TBranch *inputBranch;
int nEvents;
Event *unpackEvent;

#include "../src/tools/global_tools.cc"
#include "../src/tools/reconstruction_methods.hh"

void hitThreshold()
{

  // read in file names etc. defined on command line
  parseArguments();
  // open TFiles for input and output data and setup TTrees and TBranches
  // variables affected: FileOutput, FileInput, dataOut, data, recoEvent,
  // recoBranch, inputBranch, nEvents, unpackEvent
  linkFilesToTTrees();

  std::vector<ND280SFGDHit *> mppc_hits;

  // Specify directory for output Event Display PNG
  //string foutDir = GetDir(fileIn.Data()) + "/EventDisplay_PNG/";
  //string foutDir2 = GetDir(fileIn.Data()) + "/hitTimes_PNG/";
  //string createFolder = "mkdir -p ";
  //string createFolder2 = "mkdir -p ";
  //createFolder += foutDir.c_str();
  //createFolder2 += foutDir2.c_str();
  //system(createFolder.c_str());
  //system(createFolder2.c_str());

  //ostringstream foutPNGnum;
  //string foutPNG// set t2k style for plots
  
  TString localStyleName = "T2K";
  Int_t localWhichStyle = 3;
  TStyle *t2kstyle = SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->ForceStyle(t2kstyle);;

  cout << "total events = " << nEvents << endl;
  cout << "displaying " << evtFin - evtIni << " events" << endl;

  // Initialise Histograms_______________________________________________________________________________

  Int_t prevEventTime = -999;
  Int_t FEBs[18] = {0,1,2,3,4,8,9,10,11,16,17,18,19,20,24,25,26,27};
  Int_t FEBmap[28] = {0,1,2,3,4,0,0,0,5,6,7,8,0,0,0,0,9,10,11,12,13,0,0,0,14,15,16,17};
  Int_t FEBmap1[28] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
  Int_t FEBmap2[28] = {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,3,0,0,0};
  Int_t FEBmap3[28] = {0,0,0,0,1,0,0,0,2,3,4,5,0,0,0,0,0,0,6,7,8,0,0,0,0,9,10,11};
  Double_t binSize = 0.01;
  vector<Double_t> lowestHit_array;
  vector<Double_t> lowestHit_array1;
  vector<Double_t> lowestHit_peak_array1;
  vector<Double_t> lowestHit_peak_errors1;
  vector<Double_t> channels1;
  vector<Double_t> lowestHit_array2;
  vector<Double_t> lowestHit_peak_array2;
  vector<Double_t> lowestHit_peak_errors2;
  vector<Double_t> channels2;
  vector<Double_t> lowestHit_array3;
  vector<Double_t> lowestHit_peak_array3;
  vector<Double_t> lowestHit_peak_errors3;
  vector<Double_t> channels3;
  Double_t lowestHit_spill_array[1728*349];

  TH1D *lowestHit[1728];
  for (int ch=0; ch<1728; ch++){
    lowestHit_array.push_back(1e5);
    if (FEBs[(int)ch/96] == 0 || FEBs[(int)ch/96] == 16)
      channels1.push_back(ch);
    else if (FEBs[(int)ch/96] == 1 || FEBs[(int)ch/96] == 2 || FEBs[(int)ch/96] == 17 || FEBs[(int)ch/96] == 24)
      channels2.push_back(ch);
    else
      channels3.push_back(ch);

    if (ch < 24*48){
      lowestHit_array3.push_back(1e5);
      lowestHit_peak_array3.push_back(0);
      lowestHit_peak_errors3.push_back(0);
      if (ch < 48*8){
        lowestHit_array2.push_back(1e5);
        lowestHit_peak_array2.push_back(0);
        lowestHit_peak_errors2.push_back(0);
        if (ch < 24*8){
          lowestHit_array1.push_back(1e5);
          lowestHit_peak_array1.push_back(0);
          lowestHit_peak_errors1.push_back(0);
        }
      }
    }
    int FEB = FEBs[(int)(ch/96)];
    int Fch = ch - 96*(int)(ch/96);
    lowestHit[ch] = new TH1D(("FEB"+to_string(FEB)+"_channel"+to_string(Fch)).c_str(), ("FEB"+to_string(FEB)+" channel "+to_string(Fch)+" HG hit distribution;HG charge value [p.e.];Number of Entries").c_str(), 1000/binSize, 0, 1000);
    for (int sp=0; sp<349; sp++)
      lowestHit_spill_array[ch*349+sp-2] = 1e5;
  } 
  TH1D *lowestHit_spill[1728];
  for (int ch=0; ch<1728; ch++){
    int FEB = FEBs[(int)(ch/96)];
    int Fch = ch - 96*(int)(ch/96);
    lowestHit_spill[ch] = new TH1D(("FEB"+to_string(FEB)+"_channel"+to_string(Fch)+"_spill").c_str(), ("FEB"+to_string(FEB)+" channel "+to_string(Fch)+" HG hch distribution;HG charge value [p.e.];Number of Entries").c_str(), 50/0.1, 0, 50);
  }
  TH1D *hitDist = new TH1D("hitDist","Distribution of lowest HG hit value for each channel;HG charge value [p.e.];Number of Entries", 100, 0,10);
  TH1D *lowestHit_all = new TH1D("lowestHit_all","Distribution of lowest HG hit value for each channel;HG charge value [p.e.];Number of Entries", 100, 0,10);
  TH1D *hitDist1 = new TH1D("hitDist1","Distribution of lowest HG hit value for each channel with different MPPCs;HG charge value [p.e.];Number of Entries", 100, 0,10);
  TH1D *hitDist2 = new TH1D("hitDist2",";HG charge value [p.e.];Number of Entries", 100, 0,10);
  TH1D *hitDist3 = new TH1D("hitDist3",";HG charge value [p.e.];Number of Entries", 100, 0,10);  
  TH1D *hitDist_fine = new TH1D("hitDist_fine","Distribution of lowest HG hit value for each channel;HG charge value [p.e.];Number of Entries", 10/binSize, 0,10);
  TH1D *hitDist_fine1 = new TH1D("hitDist_fine1","Distribution of lowest HG hit value for each channel with different MPPCs;HG charge value [p.e.];Number of Entries", 10/binSize, 0,10);
  TH1D *hitDist_fine2 = new TH1D("hitDist_fine2",";HG charge value [p.e.];Number of Entries", 10/binSize, 0,10);
  TH1D *hitDist_fine3 = new TH1D("hitDist_fine3",";HG charge value [p.e.];Number of Entries", 10/binSize, 0,10);
  // TH1D *lowVchannel1 = new TH1D("lowVchannell","The lowest HG hit value per channel;Channel;Lowest HG charge value [p.e.]",1728,0,1728);
  // TH1D *lowVchannel2 = new TH1D("lowVchannel2",";Channel;Lowest HG charge value [p.e.]",1728,0,1728);
  // TH1D *lowVchannel3 = new TH1D("lowVchannel3",";Channel;Lowest HG charge value [p.e.]",1728,0,1728);

  // int window_size = 0.2e8;
  // int window[(int)(1.8e8/window_size-2)];
  // for (int i=1; i<1.8e8/window_size-1; i++)
    // window[i] = window_size*i;
  
  gStyle->SetOptStat(0);
  //loop over events
  for (int iev = evtIni; iev < evtFin; iev++)
  {
    if (iev == maxEvents - 1)
      break;
    cout << "_Processing events..." << (int)(iev * 100 / (Double_t)(evtFin - evtIni)) << "% done         \r" << flush;

    mppc_hits = getEventMPPChits(iev);

    //loop over hits
    for (size_t ihit = 0; ihit < mppc_hits.size(); ihit++)
    {
      ND280SFGDHit *hit = mppc_hits[ihit];
      //ignore hits with zero or negative charge
      if (hit->GetHG_pe() <= 0.) continue; 

      //build histogram of HG values for each channel
      lowestHit[FEBmap[(int)hit->GetFEB()]*96+(int)hit->GetCh()]->Fill(hit->GetHG_pe());

      if (lowestHit_array[(int)(FEBmap[(int)hit->GetFEB()]*96+hit->GetCh())] > hit->GetHG_pe())
        lowestHit_array[(int)(FEBmap[(int)hit->GetFEB()]*96+hit->GetCh())] = hit->GetHG_pe();

      //build arrays of lowest HG value for each channel, one for each fibre length
      if (hit->GetFEB() == 0 || hit->GetFEB() == 16){ //longest fibres
        if (lowestHit_array1[(int)(FEBmap1[(int)hit->GetFEB()]*96+hit->GetCh())] > hit->GetHG_pe())
          lowestHit_array1[(int)(FEBmap1[(int)hit->GetFEB()]*96+hit->GetCh())] = hit->GetHG_pe();
      }
      else if (hit->GetFEB() == 1 || hit->GetFEB() == 2 || hit->GetFEB() == 17 || hit->GetFEB() == 24){ //long fibres
        if (lowestHit_array2[(int)(FEBmap2[(int)hit->GetFEB()]*96+hit->GetCh())] > hit->GetHG_pe())
          lowestHit_array2[(int)(FEBmap2[(int)hit->GetFEB()]*96+hit->GetCh())] = hit->GetHG_pe();
      }
      else{ //short fibres
        if (lowestHit_array3[(int)(FEBmap3[(int)hit->GetFEB()]*96+hit->GetCh())] > hit->GetHG_pe())
          lowestHit_array3[(int)(FEBmap3[(int)hit->GetFEB()]*96+hit->GetCh())] = hit->GetHG_pe();
      }

      //build array of lowest HG value per channel per spill
      if (lowestHit_spill_array[(int)((FEBmap[(int)hit->GetFEB()]*96+hit->GetCh())*349+hit->GetSpillTag()-2)] > hit->GetHG_pe())
        lowestHit_spill_array[(int)((FEBmap[(int)hit->GetFEB()]*96+hit->GetCh())*349+hit->GetSpillTag()-2)] = hit->GetHG_pe();
    }
    // prevEventTime = unpackEvent->GetFEB12hitTfromSpill();
    //release memory used for hit information
    for (size_t i=0;i<mppc_hits.size();i++) delete mppc_hits[i];
  }
  cout << "_Processing events...100% done" << endl;

  TCanvas *cc = new TCanvas("cc","cc",0,0,800,630);
  TDirectory *HGplots = FileOutput->mkdir("HGplots");
  TDirectory *lowHGplots = FileOutput->mkdir("lowHGplots");
  HGplots->cd();
  cc->cd();
  
  //build histograms of lowest HG values, one for each MPPC type
  for (int ch=0; ch<1728; ch++){
    cc->Clear();
    int FEBind = (int)(ch/96);
    lowestHit[ch]->Draw();
    lowestHit[ch]->Write();
    lowestHit_all->Fill(lowestHit_array[ch]);
    // cout << lowestHit[it]->FindFirstBinAbove(0)*binSize-binSize << endl;
    if (FEBs[FEBind] == 8 || FEBs[FEBind] == 11){ //type 3 mppcs
      // lowVchannel3->Fill(it,lowestHit[it]->FindFirstBinAbove(0)*binSize-binSize);
      hitDist3->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
      hitDist_fine3->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
    }
    else if (FEBs[FEBind] == 3 || FEBs[FEBind] == 4 || FEBs[FEBind] == 9 || FEBs[FEBind] == 10){ //type 2 mppcs
      hitDist2->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
      hitDist_fine2->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
      // lowVchannel2->Fill(it,lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
    }
    else{ //type 1 mppcs
      hitDist1->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
      hitDist_fine1->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
      // lowVchannel1->Fill(it,lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
    }
    hitDist->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
    hitDist_fine->Fill(lowestHit[ch]->FindFirstBinAbove(0)*binSize-binSize);
    for (int sp=0; sp<349; sp++)
      lowestHit_spill[ch]->Fill(lowestHit_spill_array[ch*349+sp]);
  }

  lowestHit_all->Draw();
  lowestHit_all->Write();

  lowHGplots->cd();
  for (int ch=0; ch<1728; ch++){
    int FEBind = (int)(ch/96);
    lowestHit_spill[ch]->Draw();
    lowestHit_spill[ch]->Write();
    if (lowestHit_spill[ch]->GetStdDev() > 1) continue;
    if (FEBs[FEBind] == 0 || FEBs[FEBind] == 16){ //longest fibres
      lowestHit_peak_array1[FEBmap1[FEBs[FEBind]]*96+ch%96] = lowestHit_spill[ch]->GetMean();
      lowestHit_peak_errors1[FEBmap1[FEBs[FEBind]]*96+ch%96] = lowestHit_spill[ch]->GetStdDev();
    }
    else if (FEBs[FEBind] == 1 || FEBs[FEBind] == 2 || FEBs[FEBind] == 17 || FEBs[FEBind] == 24){ //long fibres
      lowestHit_peak_array2[FEBmap2[FEBs[FEBind]]*96+ch%96] = lowestHit_spill[ch]->GetMean();
      lowestHit_peak_errors2[FEBmap2[FEBs[FEBind]]*96+ch%96] = lowestHit_spill[ch]->GetStdDev();
    }
    else{ //short fibres
      lowestHit_peak_array3[FEBmap3[FEBs[FEBind]]*96+ch%96] = lowestHit_spill[ch]->GetMean();
      lowestHit_peak_errors3[FEBmap3[FEBs[FEBind]]*96+ch%96] = lowestHit_spill[ch]->GetStdDev();
    }
  }
  
  // Plotting Events__________________________________________________//______________________________

  TGraph *lowVchannel1 = new TGraph(channels1.size(), &channels1[0], &lowestHit_array1[0]);
  TGraph *lowVchannel2 = new TGraph(channels2.size(), &channels2[0], &lowestHit_array2[0]);
  TGraph *lowVchannel3 = new TGraph(channels3.size(), &channels3[0], &lowestHit_array3[0]);
  TGraphErrors *lowVchannelPeak1 = new TGraphErrors(channels1.size(), &channels1[0], &lowestHit_peak_array1[0], 0, &lowestHit_peak_errors1[0]);
  TGraphErrors *lowVchannelPeak2 = new TGraphErrors(channels2.size(), &channels2[0], &lowestHit_peak_array2[0], 0, &lowestHit_peak_errors2[0]);
  TGraphErrors *lowVchannelPeak3 = new TGraphErrors(channels3.size(), &channels3[0], &lowestHit_peak_array3[0], 0, &lowestHit_peak_errors3[0]);

  FileOutput->cd();
  TCanvas *cc1 = new TCanvas("cc1","cc1",0,0,800,630);
  cc1->cd();
  hitDist->Draw("hist");
  hitDist->GetXaxis()->SetRangeUser(0,6);
  cc1->Write("hitDist.pdf");
  cc1->SaveAs("hitDist.C");
  cc1->Clear();
  hitDist1->SetFillColorAlpha(kYellow, 0.35);
  hitDist1->SetLineColor(kBlack);
  hitDist2->SetLineColor(kBlack);
  hitDist3->SetLineColor(kBlack);
  hitDist2->SetFillColorAlpha(kBlue, 0.35);
  hitDist3->SetFillColorAlpha(kRed, 0.35);
  hitDist1->GetXaxis()->SetRangeUser(0,6);
  hitDist1->Draw("hist");
  hitDist2->Draw("same");
  hitDist3->Draw("same");
  TLegend *legendDist = new TLegend(0.71,0.7,0.99,0.95);
  legendDist->AddEntry(hitDist1, "Type 1 MPPCs", "f");
  legendDist->AddEntry(hitDist2, "Type 2 MPPCs", "f");
  legendDist->AddEntry(hitDist3, "Type 3 MPPCs", "f");
  legendDist->Draw("same");
  cc1->Update();
  cc1->Write("hitDist_types.pdf");
  cc1->SaveAs("hitDist_types.C");
  cc1->Clear();
  hitDist_fine->Draw("hist");
  hitDist_fine->GetXaxis()->SetRangeUser(0,6);
  cc1->Write("hitDist_fine.pdf");
  cc1->SaveAs("hitDist_fine.C");
  // cc1->WaitPrimitive();
  hitDist_fine1->SetFillColorAlpha(kYellow, 0.35);
  hitDist_fine1->SetLineColor(kBlack);
  hitDist_fine2->SetLineColor(kBlack);
  hitDist_fine3->SetLineColor(kBlack);
  hitDist_fine2->SetFillColorAlpha(kBlue, 0.35);
  hitDist_fine3->SetFillColorAlpha(kRed, 0.35);
  hitDist_fine1->GetXaxis()->SetRangeUser(0,6);
  hitDist_fine1->Draw("hist");
  hitDist_fine2->Draw("same");
  hitDist_fine3->Draw("same");
  TLegend *legendDist_fine = new TLegend(0.71,0.7,0.99,0.95);
  legendDist_fine->AddEntry(hitDist_fine1, "Type 1 MPPCs", "f");
  legendDist_fine->AddEntry(hitDist_fine2, "Type 2 MPPCs", "f");
  legendDist_fine->AddEntry(hitDist_fine3, "Type 3 MPPCs", "f");
  legendDist_fine->Draw("same");
  cc1->Update();
  cc1->Write("hitDist_types_fine.pdf");
  cc1->SaveAs("hitDist_types_fine.C");

  gStyle->SetOptStat(0);
  TCanvas *cc2 = new TCanvas("cc2","cc2",0,0,2000,630);
  cc2->cd();
  //gStyle->SetPalette(0);
  cc2->Update(); 
  lowVchannel3->Draw("ap");
  // lowVchannel1->SetTitle("; Channel; Lowest HG charge value [p.e.]");
  lowVchannel3->GetYaxis()->SetRangeUser(0.1,6.0);
  lowVchannel3->GetXaxis()->SetLimits(0,1728);
  lowVchannel3->GetXaxis()->SetTitle("Channel");
  lowVchannel3->GetYaxis()->SetTitle("Lowest HG charge value [p.e.]");
  lowVchannel3->SetTitle("");
  lowVchannel3->SetMarkerStyle(kPlus);
  lowVchannel3->SetMarkerSize(0.5);
  lowVchannel2->SetMarkerStyle(kMultiply);
  lowVchannel2->SetMarkerSize(0.5);
  lowVchannel1->SetMarkerStyle(kStar);
  lowVchannel1->SetMarkerSize(0.8);
  lowVchannelPeak1->SetMarkerSize(0.6);
  lowVchannelPeak2->SetMarkerSize(0.6);
  lowVchannelPeak3->SetMarkerSize(0.6);
  lowVchannelPeak1->SetLineColor(kRed);
  lowVchannelPeak2->SetLineColor(kRed);
  lowVchannelPeak3->SetLineColor(kRed);
  cc2->Update();
  TBox *type1_1 = new TBox(-0.5,0.1,287.5,6.0);
  type1_1->SetFillColorAlpha(kYellow, 0.35);
  TBox *type1_2 = new TBox(863.5,0.1,1727.5,6.0);
  type1_2->SetFillColorAlpha(kYellow, 0.35);
  TBox *type2_1 = new TBox(287.5,0.1,479.5,6.0);
  type2_1->SetFillColorAlpha(kBlue, 0.35);
  TBox *type2_2 = new TBox(575.5,0.1,767.5,6.0);
  type2_2->SetFillColorAlpha(kBlue, 0.35);
  TBox *type3_1 = new TBox(479.5,0.1,575.5,6.0);
  type3_1->SetFillColorAlpha(kRed, 0.35);
  TBox *type3_2 = new TBox(767.5,0.1,863.5,6.0);
  type3_2->SetFillColorAlpha(kRed, 0.35);
  TLegend *legend = new TLegend(0.91,0.7,0.99,0.95);
  legend->AddEntry(type1_1, "MPPC Type 1", "f");
  legend->AddEntry(type2_1, "MPPC Type 2", "f");
  legend->AddEntry(type3_1, "MPPC Type 3", "f");
  legend->AddEntry(lowVchannel1, "48 cm Fibres", "p");
  legend->AddEntry(lowVchannel3, "24 cm Fibres", "p");
  legend->AddEntry(lowVchannel3, "8 cm Fibres", "p");
  type1_1->Draw("same");
  type1_2->Draw("same");
  type2_1->Draw("same");
  type2_2->Draw("same");
  type3_1->Draw("same");
  type3_2->Draw("same");
  lowVchannel1->Draw("same p");
  lowVchannel2->Draw("same p");
  lowVchannel3->Draw("same p");
  TLine *line = new TLine(0,0,0,0);
  for (Int_t l=1; l<18; l++){
    line = new TLine(96*l-0.5, 0.1, 96*l-0.5, 6.0);
    line->SetLineStyle(9);
    line->Draw("same");
  }
  TGaxis *axis = new TGaxis(-0.5,6.0,1727.5,6.0,0,27,18,"-NMS");
  for (Int_t i=0; i<18; i++)
    axis->ChangeLabel(i+1,-1,-1,-1,-1,-1, to_string(FEBs[i]));
  axis->SetTitle("FEB Number");
  axis->SetTitleOffset(1.2);
  axis->SetLabelOffset(0.);
  axis->SetTickSize(-0.03);
  axis->Draw("same");
  legend->Draw("same");
  cc2->Modified(); cc2->Update();
  // cc2->Write("lowVChannel.pdf");
  cc2->SaveAs("lowVChannel.C");
  // cc2->WaitPrimitive();
  legend->AddEntry(lowVchannelPeak1, "Mean and std dev", "ep");
  legend->Draw("same");
  lowVchannelPeak1->Draw("same p");
  lowVchannelPeak2->Draw("same p");
  lowVchannelPeak3->Draw("same p");
  cc2->Modified(); cc2->Update();
  // cc2->Write("lowVChannel.pdf");
  cc2->SaveAs("lowVChannelPeak.C");

  // cc2->Update();

  cc->Close();
  cc1->Close();
  // cc2->Close();

  //cc->Clear();
  //timeSpec->Draw();

  //cc->Update();
  cc2->WaitPrimitive();

  //cc->SaveAs((foutDir+"/event_"+to_string(iev)+".png").c_str());

  // FileOutput->Close();
  // FileInput->Close();
  writeOutput();
  handleEndOfExecution();
}
