//
//
// César Jesús-Valls (cjesus@ifae.es)
// 19-10-2019
//
// macro to analyze the ReconstructVoxelsPE performance.
//
//

#define THIS_NAME lightReconstruction
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
TString fileIn         = "/Users/cjesus/Documents/PhD/SFGD/ROOT_DATA_v2/PI_MINS.root";
int IsPrototype        = false;
int SFGD_X             = 204;
int SFGD_Y             = 56;
int SFGD_Z             = 188;
int dTimeIni           = -100;
int dTimeFin           = -100;
int evtIni             = 0;
int evtFin             = 0;
bool useNN             = false;
bool useClean          = true;

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

void lightReconstruction() {

    TH1F* h_PE_resolution = new TH1F("PE resolution","PE resolution",100,-100,100);
    TH2F* h_PE_migration  = new TH2F("PE migration","PE migration",20,0,400,20,0,400);

    parseArguments();
    linkFilesToTTrees();

    std::vector<ND280SFGDHit*> mppc_hits;

    cout << "Initial Event: " << evtIni << endl;
    cout << "Final Event: " << evtFin << endl;
    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;

        mppc_hits = getEventMPPChits(iev); 

        // fill voxels, and link its information to the true event, to be able to plot true vs reconstructed distrubutions.
        recoEvent->SetHits(MergeHits(mppc_hits,0,true));
        recoEvent->SetVoxels(HitsToVoxels(recoEvent->GetHits(),0));
        convertCoordinates(inputEvent->GetVoxels());

        FillTrueInformationInRecoEvent(inputEvent,recoEvent,0);
        ReconstructVoxelsPE(recoEvent->GetHits(),recoEvent->GetVoxels(),0);

        for(auto voxel:recoEvent->GetVoxels()){
            h_PE_migration->Fill(voxel->GetTruePE(),voxel->GetRecoPE());
            h_PE_resolution->Fill(100.*(voxel->GetTruePE()-voxel->GetRecoPE())/voxel->GetTruePE());
        }
        selEvents++;
    }

    TCanvas c;
    c.Divide(2,1);
    c.cd(1);
    h_PE_migration->GetXaxis()->SetTitle("(True - Reco)/True [#pe]");
    h_PE_resolution->Draw("HIST");

    c.cd(2);
    h_PE_migration->GetXaxis()->SetTitle("True #pe");
    h_PE_migration->GetYaxis()->SetTitle("Reco #pe");
    h_PE_migration->Draw("COLZ");

    c.Update();
    c.WaitPrimitive();

    handleEndOfExecution();
}
