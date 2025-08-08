//
//
// César Jesús-Valls (cjesus@ifae.es)
// 10-10-2019
//
// template file to show how to construct a loop that iterates
// over MC/data events, applies or not some reconstruction and does some analysis
//
//

#define THIS_NAME templateAna
#define NOINTERACTIVE_OUTPUT           
#define OVERRIDE_OPTIONS 

#include "../src/tools/global_header.hh"

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool    batch          = false;
bool    IsMC           = true; 
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "~/Desktop/templateAna_Out.root";
TString fileIn         = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
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

void templateAna() {
    
    /// ------------START--------------

    parseArguments();
    linkFilesToTTrees();

    std::vector<ND280SFGDHit*> mppc_hits;

    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;

        mppc_hits = getEventMPPChits(iev); 
        
        //________ANALYSIS_OF_EVT_STARTS_____//
        //
        //  ---> Write here your analysis.

        //e.g:
        if(IsMC){
            cout << "summary of true information: [without merging]" << endl;
            cout << "# of voxels: "   << inputEvent->GetVoxels().size() << endl;
            cout << "# of hits:   "   << inputEvent->GetHits().size()   << endl;
            cout << "# of tracks: "   << inputEvent->GetTrueTracks().size() << endl;
        }

        inputEvent->DrawHitsAndVoxels(true,"event");
        inputEvent->DrawVoxelsTruePE(true,"truePE");

        //
        //________ANALYSIS_OF_EVT_ENDS________//

        dataOut->Fill();
        selEvents++;
    }

    writeOutput();
    handleEndOfExecution();
}
