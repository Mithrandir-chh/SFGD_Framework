//
//
// César Jesús-Valls (cjesus@ifae.es)
// 10-10-2019
//
// template file to show how to construct a loop that iterates
// over MC/data events, applies or not some reconstruction and does some analysis
//
//

#define THIS_NAME crosstalkAna
#define NOINTERACTIVE_OUTPUT           
#define OVERRIDE_OPTIONS 

#include "../src/tools/global_header.hh"

//____DEFINE_GLOBAL_SCOPE_VARIABLES____: this needs to be declared before global_tools!
bool    batch          = false;
bool    IsMC           = false; 
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    RM_CROSSTALK   = false;
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "../../../results/crosstalk/crosstalkAna.root";
TString fileIn         = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
int IsPrototype        = false;
int SFGD_X             = 24;
int SFGD_Y             = 8;
int SFGD_Z             = 48;
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

void crosstalkAna() {
    
    /// ------------START--------------

    int VERBOSE = 2;

    if(!IsMC) fileIn = "/home/cjesus/Work/Data/SFGD_prototype/DATA/26August_ALL/26AugustAll.root";

    parseArguments();
    linkFilesToTTrees();

    std::vector<ND280SFGDHit*> mppc_hits;

    for (int iev=0; iev<nEvents; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;
        
        //________ANALYSIS_OF_EVT_STARTS_____//
        //

        cout << "Event" << iev << endl;
        mppc_hits = getEventMPPChits(iev); 
        cout << "# hits: " << mppc_hits.size() << endl;
        if (VERBOSE>1) for (auto h:mppc_hits) cout << "V-XYZ: " << h->GetView() << "," << h->GetX() << "," << h->GetY() << "," << h->GetZ() << endl; 

        //
        //________ANALYSIS_OF_EVT_ENDS________//

        selEvents++;
    }

    handleEndOfExecution();
}
