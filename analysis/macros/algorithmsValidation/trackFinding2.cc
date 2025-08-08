#define THIS_NAME trackFinding2
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
bool useClean          = true;
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

void trackFinding2() {
    
    /// ------------START--------------

    parseArguments();
    linkFilesToTTrees();

    std::vector<ND280SFGDHit*> mppc_hits;

    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;

        mppc_hits = getEventMPPChits(iev); 
        recoEvent->SetHits(MergeHits(mppc_hits,0));
        recoEvent->SetVoxels(HitsToVoxels(recoEvent->GetHits(),0));
        convertCoordinates(inputEvent->GetVoxels());
        FillTrueInformationInRecoEvent(inputEvent,recoEvent,0);

        int topology = 0;
        for(auto t:inputEvent->GetTrueTracks()) {if(t->GetParentID() == 0) topology++;}
        cout << "topology: " << topology << endl;


        std::vector<ND280SFGDVoxelSet*> clusters = FindClusters(recoEvent->GetVoxels(),0);
        ND280SFGDVoxelSet* largestClust = new ND280SFGDVoxelSet();
        for(auto cluster:clusters){
            if(cluster->GetVoxels().size() > largestClust->GetVoxels().size()) largestClust = cluster;
        }


        if(useNN){
            std::vector<ND280SFGDVoxel*> v_trk;
            v_trk = ApplyFakeNN(largestClust->GetVoxels());
            largestClust->SetVoxels(v_trk);
        }

        if(useClean){
            std::vector<ND280SFGDVoxel*> v_trk;
            v_trk = cleanVoxels(largestClust->GetVoxels());
            largestClust->SetVoxels(v_trk);
        }

        if(largestClust->GetVoxels().size()<6) continue;


        std::vector<ND280SFGDVoxel*> branchingPoints;
        branchingPoints = FindBranchingPoints(largestClust->GetVoxels(),0);

        cout << "number of bPoints: " << branchingPoints.size() << endl;


        
        //inputEvent->DrawHitsAndVoxels(true,"event");
        //recoEvent->DrawHitsAndVoxels(true,"event");
        largestClust->DrawTrueTracks(true,"trueTracks");
        //largestClust->DrawVoxelsTruePE(true,"event");
        selEvents++;
    }

    handleEndOfExecution();
}
