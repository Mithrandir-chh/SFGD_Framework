//
//
// César Jesús-Valls (cjesus@ifae.es)
// 19-10-2019
//
// macro to store MC data to csv files to perform the neuralNet studies for SFGD.

#define THIS_NAME inputNN
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
bool useNN             = false;
bool useClean          = true;
int evtIni             = 0;
int evtFin             = 0;
int dTimeIni           = -100;
int dTimeFin           = -100;

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

void inputNN() {
    
    /// ------------START--------------

    bool DEBUG = true;

    parseArguments();
    linkFilesToTTrees();

    std::vector<ND280SFGDHit*> mppc_hits;

    std::ofstream outCSVfile; 
    outCSVfile.open ("NN.csv",std::ofstream::out);
    
    //print the column headers
    outCSVfile<<"eventID"<<", "<<"X"<<", "<<"Y"<<", "<<"Z"<<", "
                <<"Qxy"<<", "<<"Qxz"<<", "<<"Qyz"<<", "<<"Mxy"<<", "<<"Mxz"<<", "<<"Myz"<<", "
                <<"QxyCorr"<<", "<<"QxzCorr"<<", "<<"QyzCorr"<<", "<<"Chi2Corr"<<", "
                <<"PullXCorr"<<", "<<"PullYCorr"<<", "<<"PullZCorr"<<", "<<"label"<< "\n";

    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;

        mppc_hits = getEventMPPChits(iev); 
        recoEvent->SetHits(MergeHits(mppc_hits,0));
        recoEvent->SetVoxels(HitsToVoxels(recoEvent->GetHits(),0));
        convertCoordinates(inputEvent->GetVoxels());
        FillTrueInformationInRecoEvent(inputEvent,recoEvent,0);

        if(DEBUG){
            
            int Ntype [3]    = {0,0,0};   // number of each true type of voxels
            int NemptyTruePE = 0;         // voxels that have 0 true PE but are not ghosts
            int NfakeGhosts  = 0;         // voxels that are true ghosts but have multiplicity < 2 in some of its fibers.
            int NisolatedXtalk   = 0;     // voxels that are true crosstalk but have no true track voxels around.
            for (auto n:recoEvent->GetVoxels()){
                // count how many voxels of each type there is in every event.
                Ntype[n->GetTrueType()]++;
                // check that any voxel has true Q == 0 and is not a ghost.
                if (!n->GetTruePE() && n->GetTrueType() !=2) NemptyTruePE++;
                // check that all voxels have at least multiplicity 2 in all its fibers.
                if (n->GetTrueType() == 2)
                    if (n->GetHits()[0]->GetMultiplicity()<2 || n->GetHits()[1]->GetMultiplicity()<2 || n->GetHits()[2]->GetMultiplicity()<2)
                        {NfakeGhosts++; cout << "trueType: " << n->GetTrueType() << " -- Mxy|Mxz|Myz: " << n->GetHits()[0]->GetMultiplicity() << "," << n->GetHits()[1]->GetMultiplicity() << "," << n->GetHits()[2]->GetMultiplicity() << " -- PExy|PExz|PEyz: " << n->GetHits()[0]->GetPE() << "," << n->GetHits()[1]->GetPE() << "," << n->GetHits()[2]->GetPE() << endl;}
                // check that all crosstalk voxels have a true voxel at a distance of 1cm
                
                if(n->GetTrueType() ==1){
                    bool found = false;
                    for (auto n_other:recoEvent->GetVoxels()) if (n != n_other && n->DistToVoxel(n_other) == 1 ) {found = true; break;}
                    if(!found) NisolatedXtalk++;
                }
            }
            
            cout << " --- Event " << selEvents << endl;
            cout << "TOTAL     [track-xtalk-ghost]:                " << Ntype[0] << "," << Ntype[1] << "," << Ntype[2] << endl; 
            cout << "percent   [track-xtalk-ghost]:                " << setprecision(3) << 100.*Ntype[0]/recoEvent->GetVoxels().size() << "," << setprecision(3) << 100.*Ntype[1]/recoEvent->GetVoxels().size() << "," << setprecision(3) << 100.*Ntype[2]/recoEvent->GetVoxels().size() << endl; 
            cout << "voxels [no-ghost] with true PE = 0:           " << NemptyTruePE << endl;
            cout << "ghosts with at least 1 multiplicity 1 fiber:  " << NfakeGhosts << endl;
            cout << "crosstalk without true track deposits around: " << NisolatedXtalk << endl;
            cout << " ------------- " << endl;

        }

        DumpToCSVfile(outCSVfile, recoEvent, selEvents);
        selEvents++;
    }

    outCSVfile.close();

    handleEndOfExecution();
}
