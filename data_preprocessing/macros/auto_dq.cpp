/* This file is part of BabyMINDdaq software package. This software
 * package is designed for internal use for the Baby MIND detector
 * collaboration and is tailored for this use primarily.
 *
 * BabyMINDdaq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BabyMINDdaq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BabyMINDdaq.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// ROOT libraries
#include "TROOT.h"
#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObject.h"
#include <TStyle.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TClass.h>
#include <TChain.h>
#include <TLegend.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Fit/Fitter.h>
#include <TRandom3.h>
#include <TEfficiency.h>

//C++ libraries
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <array>
#include <time.h>
#include <map>

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


#include <stdio.h>
#include <string.h>
#include <exception>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TMacro.h"
#include <TTree.h>
#include "MDfragmentBM.h"
#include "MDpartEventBM.h"
#include "MDargumentHandler.h"
#include "MDdataFile.h"

using namespace std;

string GetLocation(string str)
{

    int i = str.rfind("Slot");
    string way = str.substr(0,i);
    return way;
}

struct vectorsTree
{
  vector<double> FEBSN;
  vector<double> SpillNum;
  vector<double> SpillTime;
  vector<double> SpillTimeGTrig;
  vector<double> GTrigTag;
  vector<double> GTrigTime;
  vector<double> hitsChannel;
  vector<double> hitAmpl;
  vector<double> hitLGAmpl;
  vector<double> hitLeadTime;
  vector<double> hitTrailTime;
  vector<double> hitTimeDif;
  vector<double> hitTimefromSpill;
  vector<double> SpillTrailTime;
 // vector<double> SpillTemperature;
  vector<double> AsicTemperature;
  vector<double> FPGATemperature;
  vector<double> GlobalHV;
  vector<double> BoardTemperature;
  vector<double> BoardHumidity;
};



char *dataBuff;
uint32_t* dataPtr;

int auto_dq() {
  string sFileName;
  int maxSpills = -1;

  TH2F* h2D[30];
  TH1F* h1D[30];
  for (Int_t i=0;i<30;i++){
    h2D[i] = new TH2F("","Channel vs HG Counts; Channel number; HG Counts",400,0,2000,96,0,96);
    h1D[i] = new TH1F("","Hit time from Spill; Hit time from spill; Counts",100,0,4.E8);
  }
  int FEBs[24] = {0,1,2,3,8,9,10,11,12,16,17,18,19,20,24,25,26,27,56,57,58,59,60,61};

  vector<string> vFileNames;
  ifstream fList("febs_files_list.list");
  while (!fList.eof()) {
    fList >> sFileName;
    //cout << sFileName << endl;
    vFileNames.push_back(sFileName);
  }
  vFileNames.pop_back();
  for (int i = 0; i< vFileNames.size();i++){
      cout << vFileNames.at(i)<<endl;
  }
  int NumberOfFEB = 65; // large number will account for FEBs on US-Japan prototype if included
  vectorsTree FEB[NumberOfFEB];

  for (vector<string>::iterator itFileName=vFileNames.begin(); itFileName != vFileNames.end(); itFileName++) {
    sFileName = *itFileName;
    cout <<endl<< sFileName << endl;

    ifstream finData(sFileName.c_str());
    string stringBuf;
    string filepath;
    string filename;

    if (maxSpills < 0) cout << "Unpacking all spills" << endl;
    //if ( argh.GetValue("directory", stringBuf) != MDARGUMENT_STATUS_OK ) return -1;
    else cout << "Unpacking " << maxSpills << " spills" << endl;
    //filepath = stringBuf;
    //if ( argh.GetValue("file", stringBuf) != MDARGUMENT_STATUS_OK ) return -1;
    //filename = stringBuf;
    filename = sFileName;
    string rootFilename = sFileName;

    //TFile rfile("histos.root", "recreate");

    for (Int_t i=0;i<NumberOfFEB;i++){
      FEB[i].FEBSN.clear();
      FEB[i].SpillNum.clear();
      FEB[i].SpillTime.clear();
      FEB[i].SpillTimeGTrig.clear();
      FEB[i].hitsChannel.clear();
      FEB[i].hitAmpl.clear();
      FEB[i].hitLeadTime.clear();
      FEB[i].GTrigTag.clear();
      FEB[i].GTrigTime.clear();
      FEB[i].hitLGAmpl.clear();
      FEB[i].hitTrailTime.clear();
      FEB[i].hitTimeDif.clear();
      FEB[i].hitTimefromSpill.clear();
      FEB[i].SpillTrailTime.clear();
      //FEB[i].SpillTemperature.clear();
      FEB[i].AsicTemperature.clear();
      FEB[i].FPGATemperature.clear();
      FEB[i].GlobalHV.clear();
      FEB[i].BoardTemperature.clear();
      FEB[i].BoardHumidity.clear();
    }
    MDdateFile dfile(filename, filepath);
    // Open the file and loop over events.
    Int_t BordID;
    char *eventBuffer;
    bool _previousSpillTagExist = false;
    int _previousSpillTag = 0;
    int prevSpill = 0;
    int spillTag = 0;
    unsigned int _spillsProcessed = 0;
    if ( dfile.open() ) { // There is a valid file to unpack
      dfile.init();
      eventBuffer = dfile.GetNextEvent();
      int xEv(0);
      do { // Loop over all spills
        eventBuffer =  dfile.GetNextEvent();
        _spillsProcessed += 1;
        if (_spillsProcessed == maxSpills) break;
        try {
          MDfragmentBM   spill;

          spill.SetPreviousSpill(_previousSpillTagExist,_previousSpillTag);
          prevSpill = _previousSpillTag;
          spill.SetDataPtr(eventBuffer);
          spill.Init(spillTag);

          MDpartEventBM *event;
          int nTr = spill.GetNumOfTriggers();
          BordID = (Int_t)spill.GetBoardId();

          for (int i=0; i<nTr; ++i) {
            event = spill.GetTriggerEventPtr(i);
            //event->Dump();
            for (int ich=0; ich<BM_FEB_NCHANNELS; ++ich) {
              int nlHits = 0;
              int ntHits = 0;
              nlHits = event->GetNLeadingEdgeHits(ich);
              for (unsigned int ih=0; ih<nlHits; ++ih) {

                _previousSpillTag = spill.GetSpillTag();

		int ihere = -1;
		for (int loop = 0;loop< 24; loop++){
		  if(FEBs[loop] == spill.GetBoardId() ){
		    ihere = loop;
		    break;
		  }
		}
		h2D[ihere] -> Fill(ich, event->GetHitAmplitude(ich, 'h')); 

              }
            }
          }
        } 


        catch (MDexception & lExc) {
          std::cerr <<  lExc.GetDescription() << endl
                      << "Unpacking exception\n"
                      << -prevSpill+spillTag-1 << " spills missed!\n\n"
                      << "Here they should fill with negative numbers. \n\n";

          for (int misses=prevSpill; misses<spillTag;misses++){
            FEB[BordID].FEBSN.push_back(BordID);
            FEB[BordID].SpillNum.push_back(misses+1);
            FEB[BordID].SpillTime.push_back(-1);
            FEB[BordID].SpillTimeGTrig.push_back(-1);
            FEB[BordID].hitsChannel.push_back(-1);
            FEB[BordID].hitAmpl.push_back(-1);
            FEB[BordID].hitLeadTime.push_back(-1);
            FEB[BordID].GTrigTag.push_back(-1);
            FEB[BordID].GTrigTime.push_back(-1);
            FEB[BordID].hitLGAmpl.push_back(-1);
            FEB[BordID].hitTrailTime.push_back(-1);
            FEB[BordID].hitTimeDif.push_back(-1);
            FEB[BordID].hitTimefromSpill.push_back(-1);
            FEB[BordID].SpillTrailTime.push_back(-1);
            //FEB[BordID].SpillTemperature.push_back(-1);
            FEB[BordID].AsicTemperature.push_back(-1);
            FEB[BordID].FPGATemperature.push_back(-1);
            FEB[BordID].GlobalHV.push_back(-1);
            FEB[BordID].BoardTemperature.push_back(-1);
            FEB[BordID].BoardHumidity.push_back(-1);
          }

          _previousSpillTag = spillTag;

        } catch(std::exception & lExc) {
          std::cerr << lExc.what() << std::endl
                      << "Standard exception\n"
                      << "Spill skipped!\n\n";
        } catch(...) {
          std::cerr << "Unknown exception occurred...\n"
                      << "Spill skipped!\n\n";
        }
        //if (FEB[BordID].FEBSN.size()){
          cout<<"Number of events on FEB "<< BordID <<" is "<< FEB[BordID].FEBSN.size()<<endl;
          //FEBtree[BordID]->Fill();
          //FEBtree[BordID]-> Write("",TObject::kOverwrite);
          //FEBtree[BordID]-> Write();
          //FEBtree[BordID]->Delete();
          FEB[BordID].FEBSN.clear();
          FEB[BordID].SpillNum.clear();
          FEB[BordID].SpillTime.clear();
          FEB[BordID].SpillTimeGTrig.clear();
          FEB[BordID].hitsChannel.clear();
          FEB[BordID].hitAmpl.clear();
          FEB[BordID].hitLeadTime.clear();
          FEB[BordID].hitTrailTime.clear();
          FEB[BordID].hitTimeDif.clear();
          FEB[BordID].GTrigTag.clear();
          FEB[BordID].GTrigTime.clear();
          FEB[BordID].hitLGAmpl.clear();
          FEB[BordID].hitTimefromSpill.clear();
          FEB[BordID].SpillTrailTime.clear();
          //FEB[BordID].SpillTemperature.clear();
          FEB[BordID].AsicTemperature.clear();
          FEB[BordID].FPGATemperature.clear();
          FEB[BordID].GlobalHV.clear();
          FEB[BordID].BoardTemperature.clear();
          FEB[BordID].BoardHumidity.clear();
        // }

        ++xEv;
      //  } while (xEv < 5);
      } while ( eventBuffer );
    }

    //FEBtree[BordID]->Write("",TObject::kOverwrite);
    //FEBtree[BordID]->Delete();

    dfile.close();
    delete dataBuff;
  }
  //rfile.Close();

  TCanvas *cc = new TCanvas("cc", "cc", 800, 800);
  cc->Divide(5,6);
  for(int ii=0;ii<30;ii++){
    cc->cd(ii+1);
    h2D[ii]->Draw("colz");
  }
  cc->Update();
  cc->WaitPrimitive();
  
  fList.close();
  exit(0);
}

int main(int argc, char* argv[]) {
  TApplication theApp("woo",&argc, argv);
  auto_dq();
  theApp.Run();
  exit(0);
  return 0;
}

