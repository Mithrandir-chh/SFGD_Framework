//
//
// César Jesús-Valls (cjesus@ifae.es)
// 16-11-2019
//
//

#define THIS_NAME pionReco
#define NOINTERACTIVE_OUTPUT           
#define OVERRIDE_OPTIONS 

#include "../src/tools/global_header.hh"

bool    batch          = false;
bool    IsMC           = true; 
int     maxEvents      = std::numeric_limits<int>::max();;
int     maxSelEvents   = std::numeric_limits<int>::max();;
int     selEvents      = 0;
float   start_time     = clock();
bool    SHOW_TRUE      = true;
bool    SHOW_RECO      = true;
TString fileOut        = "~/Desktop/pionRecoDist.root";
TString fileIn         = "/home/cjesus/Work/Dev/SFGD/sfgd_framework/analysis/data/204_56_188/PI_MINS.root";
int IsPrototype        = false;
int SFGD_X             = 204;
int SFGD_Y             = 52;
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

void pionReco() {

    parseArguments();
    linkFilesToTTrees();

    int VERBOSE    = 0;  // Prints messages to help debuging
    int MODE       = 0;  // The MODE defines which case is studied based on true information:
                         //     0: one pion, any number of protons and no other tracks.
                         //     1: one pion, any number of protons, one muon and no other tracks. [decay]
                         //     2: two pions, any number of protons, and no other tracks.
                         //     3: more than two pions and any other number of tracks.
    int USE_EDEP   = 0;  // Reconstructs Edep and uses it to reconstruct pion showers. // TODO: needs further development.
    int STORE      = 1;  // Generates new true-reco fits and stores them in files.

    // Add extra custom command line argument options:
    for (int iarg=0; iarg<gApplication->Argc(); iarg++)
        if (string( gApplication->Argv(iarg))=="-opt"){
            cout << "ENTER VERBOSE: " << endl;
            cin  >> VERBOSE;
            cout << "ENTER MODE: " << endl;
            cin  >> MODE;
            cout << "ENTER USE_EDEP: " << endl;
            cin  >> USE_EDEP;
        }

    // In the future get the functions from dedicated input files.
    // pByRng contains the predicted Momentum given a reconstructed range, based on MC.
    // TODO: improve the fit mathematical description.
    TF1 *pByRng = new TF1("f1","[0]+sqrt([1]*x)",0,200);
    if( MODE == 0 or MODE == 1 ) pByRng->SetParameters(7.8789413,1294.4196);
    if( MODE == 2 or MODE == 3 ) pByRng->SetParameters(176.16961,632.6997);

    if(VERBOSE>0) cout << "Params: " << std::setprecision(8) << pByRng->GetParameter(0) << "," << pByRng->GetParameter(1) << endl;

    // pByEdp contains the predicted Momentum given a reconstructed edep, based on true MC.
    // NOTE: this is an option to be further investigated.
    TF1 *pByEdp = new TF1("f1","[0]+sqrt([1]*x)",0,200);
    pByEdp->SetParameters( -44.457695 ,241.58553);

    // NOTE: all counters are organized by # of pions in the final state.
    // Counters of total events, independent of the number of other particles.
    int c_OnePion = 0;
    int c_TwoPion = 0;
    int c_MltPion = 0;
    // Counters of events with 'clean' [no other true track] conditions.
    int c_OnePionClean = 0;
    int c_TwoPionClean = 0;
    int c_MltPionClean = 0;
    // Counters of events involving protons.
    int c_OnePionPrFound = 0;
    int c_TwoPionPrFound = 0;
    int c_MltPionPrFound = 0;
    // Counters of events involving decays.
    int c_OnePionMuFound = 0;
    int c_TwoPionMuFound = 0;
    int c_MltPionMuFound = 0;
    // Counters of events involving electrons.
    int c_OnePionElFound = 0;
    int c_TwoPionElFound = 0;
    int c_MltPionElFound = 0;
    // Counters of events involving positrons.
    int c_OnePionPsFound = 0;
    int c_TwoPionPsFound = 0;
    int c_MltPionPsFound = 0;
    // Counters of events involving gammas.
    int c_OnePionPhFound = 0;
    int c_TwoPionPhFound = 0;
    int c_MltPionPhFound = 0;
    // Counters of events involving decays.
    int c_OnePionMuFoundAndDet = 0;
    int c_TwoPionMuFoundAndDet = 0;
    int c_MltPionMuFoundAndDet = 0;

    // Results, organized in histograms:
    TH2F* h_EnergyCont   = new TH2F("h_EnergyCont"    ,"h_Cont",          50,0,500,50,0,500);
    TH2F* h_RngVSRng     = new TH2F("h_RngVSRng"      ,"h_RngVSRng",      20,0,100,20,0,100);
    TH2F* h_EdpVSEdp     = new TH2F("h_EdpVSEdp"      ,"h_EdpVSEdp",      20,0,100,20,0,100);
    TH1F* h_RngRes       = new TH1F("h_RngRes"        ,"h_RngRes",        100,-50,50);
    TH1F* h_EdpRes       = new TH1F("h_EdpRes"        ,"h_EdpRes",        100,-50,50);
    TH2F* h_MomVSRng     = new TH2F("h_MomVSRng"      ,"h_MomVSRng",      50,0,200,50,0,500);
    TH2F* h_MomVSEdp     = new TH2F("h_MomVSEdp"      ,"h_MomVSEdp",      50,0,500,50,0,500);
    TH2F* h_MomVSRecoRng = new TH2F("h_MomVSRecoRng"  ,"h_MomVSRecoRng",  20,0,200,20,0,500);
    TH2F* h_MomVSRecoEdp = new TH2F("h_MomVSRecoEdp"  ,"h_MomVSRecoEdp",  20,0,200,20,0,500);
    TH1F* h_MomResRng    = new TH1F("h_MomResRng"     ,"h_MomResRng",     100,-50,200);
    TH1F* h_MomResEdp    = new TH1F("h_MomResEdp"     ,"h_MomResEdp",     100,-50,200);


    TH1F* h_MomResRng_split [5];
    h_MomResRng_split [0]    = new TH1F("h_MomResRng_0"     ,"h_MomResRng_0",     100,-50,200);
    h_MomResRng_split [1]    = new TH1F("h_MomResRng_1"     ,"h_MomResRng_1",     100,-50,200);
    h_MomResRng_split [2]    = new TH1F("h_MomResRng_2"     ,"h_MomResRng_2",     100,-50,200);
    h_MomResRng_split [3]    = new TH1F("h_MomResRng_3"     ,"h_MomResRng_3",     100,-50,200);
    h_MomResRng_split [4]    = new TH1F("h_MomResRng_4"     ,"h_MomResRng_4",     100,-50,200);

    // To extract tables: [mom vs range] & [mom vs edep]
    TH2F* h_pVSr  = new TH2F("h_pVSr","h_pVSr",50,0,200,50,0,500);
    TH2F* h_pVSe  = new TH2F("h_pVSe","h_pVSe",50,0,200,50,0,500);

    // To print accumulated hits:
    TH2F *hXY = new TH2F("XY","XY", SFGD_X, 0.5, SFGD_X+0.5, SFGD_Y, 0.5, SFGD_Y+0.5);      
    TH2F *hXZ = new TH2F("XZ","XZ", SFGD_X, 0.5, SFGD_X+0.5, SFGD_Z, 0.5, SFGD_Z+0.5);
    TH2F *hYZ = new TH2F("YZ","YZ", SFGD_Y, 0.5, SFGD_Y+0.5, SFGD_Z, 0.5, SFGD_Z+0.5);

    // To store 2D mppc_hits of each event:
    std::vector<ND280SFGDHit*> mppc_hits;

    // Loop over all events:
    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;
        if(!VERBOSE) if(iev%100 ==0) cout << "\33[2K\rProcessed: \033[1;33m" << std::setprecision(2) << (100.*iev/(evtFin-evtIni)) << "%\t\033[0m" << std::flush;

        mppc_hits = getEventMPPChits(iev); 
        recoEvent->SetHits(MergeHits(mppc_hits,0));
        recoEvent->SetVoxels(HitsToVoxels(recoEvent->GetHits(),0));
        convertCoordinates(inputEvent->GetVoxels());
        FillTrueInformationInRecoEvent(inputEvent,recoEvent,0);

        // Safety check: Find if there is a true pion. [Sometimes it decays into muon before entering SFGD if pgun is OutFV].
        ND280SFGDTrack* truePion = nullptr;
        for(auto tp:inputEvent->GetTrueTracks()) if(!tp->GetParentID()) truePion = tp;
        if(!truePion){
            if (VERBOSE>1) cout << "There is no true pion!" << endl;
            continue;
        }

        // Check if the pion is contained in SFGD:
        // TODO: use pgun starting within SFGD FV.
        bool OutFV = false;
        for(auto v:inputEvent->GetVoxels()){
            if(v->GetX()<2 or (SFGD_X-v->GetX())<2) OutFV = true;
            if(v->GetY()<2 or (SFGD_Y-v->GetY())<2) OutFV = true;
            if((SFGD_Z-v->GetZ())<2) OutFV = true;
        }
        if(OutFV) continue;
        // Apparently some pions are scattered back [TODO: check this] 
        if(((truePion->GetRange())/10-25.1)>150) continue;
        selEvents++;
        
        // Check the true tracks and classify them by true PDG:
        int Npions   = 0;
        int Nprotons = 0;
        int Nelectrs = 0;
        int Npositrs = 0;
        int Nmuons   = 0;
        int Nphot    = 0;
        int Nneut    = 0;

        for(auto t:inputEvent->GetTrueTracks()){
            if(abs(t->GetPDG()) ==  211) Npions++;
            else if(abs(t->GetPDG()) ==  13)  Nmuons++;
            else if(t->GetPDG() == 2212)      Nprotons++;
            else if(t->GetPDG() == 11)        Nelectrs++;
            else if(t->GetPDG() == -11)       Npositrs++;
            else if(t->GetPDG() == 22)        Nphot++;
            else if(t->GetPDG() == 2112)      Nneut++;
        }

        if(Npions == 1) c_OnePion++;
        if(Npions == 2) c_TwoPion++;
        if(Npions >= 3) c_MltPion++;

        if(Npions == 1 && !Nmuons && !Nelectrs && !Npositrs && !Nphot) c_OnePionClean++;
        if(Npions == 2 && !Nmuons && !Nelectrs && !Npositrs && !Nphot) c_TwoPionClean++;
        if(Npions >= 3 && !Nmuons && !Nelectrs && !Npositrs && !Nphot) c_MltPionClean++;

        if(Npions == 1 && Nprotons) c_OnePionPrFound++;
        if(Npions == 2 && Nprotons) c_TwoPionPrFound++;
        if(Npions >= 3 && Nprotons) c_MltPionPrFound++;

        if(Npions == 1 && Nmuons)   c_OnePionMuFound++;
        if(Npions == 2 && Nmuons)   c_TwoPionMuFound++;
        if(Npions >= 3 && Nmuons)   c_MltPionMuFound++;

        if(Npions == 1 && Nelectrs) c_OnePionElFound++;
        if(Npions == 2 && Nelectrs) c_TwoPionElFound++;
        if(Npions >= 3 && Nelectrs) c_MltPionElFound++;

        if(Npions == 1 && Npositrs) c_OnePionPsFound++;
        if(Npions == 2 && Npositrs) c_TwoPionPsFound++;
        if(Npions >= 3 && Npositrs) c_MltPionPsFound++;

        if(Npions == 1 && Nphot)    c_OnePionPhFound++;
        if(Npions == 2 && Nphot)    c_TwoPionPhFound++;
        if(Npions >= 3 && Nphot)    c_MltPionPhFound++;
            
        std::vector<ND280SFGDVoxel*> inputVoxels = inputEvent->GetVoxels();
        bool delayedSignal = false;
        if (Nmuons) delayedSignal = std::any_of(inputVoxels.begin(), inputVoxels.end(),
                                    [&](ND280SFGDVoxel* v) {return v->GetTime()>0;}) ? true : false;

        if(Npions == 1 && delayedSignal) c_OnePionMuFoundAndDet++;
        if(Npions == 2 && delayedSignal) c_TwoPionMuFoundAndDet++;
        if(Npions >= 3 && delayedSignal) c_MltPionMuFoundAndDet++;

        // single pion w or w/o protons:
        if (VERBOSE >0 ) cout << "[Pi,Pr,Mu,El,Ps,Ph]: " << Npions << "," << Nprotons << "," << Nmuons << "," << Nelectrs << "," << Npositrs << "," << Nphot << endl;
        if(MODE == 0) if(!(Npions == 1 && !Nmuons && !Nelectrs && !Npositrs && !Nphot)) continue;
        if(MODE == 1) if(!(Npions == 1 &&  Nmuons && !Nelectrs && !Npositrs && !Nphot)) continue;
        if(MODE == 2) if(!(Npions == 2 && !Nmuons && !Nelectrs && !Npositrs && !Nphot)) continue;
        if(MODE == 3) if(!(Npions == 3 && !Nmuons && !Nelectrs && !Npositrs && !Nphot)) continue;

        // The edep reconstruction is yet computationally expensive. Do not use it for the time being.
        // TODO: Improve and speed up Edep reconstruction.
        if(USE_EDEP) ReconstructVoxelsPE(recoEvent->GetHits(),recoEvent->GetVoxels(),0);

        //NOTE: For the time being I am using PE/100 instead of MeV. 
        //TODO: Change this scale to quenched MeV.
        //NOTE: Right now this not includes any quenching correction.
        double trueEdp = 0;
        double recoEdp = 0;
        for(auto v:recoEvent->GetVoxels()){
            recoEdp += v->GetRecoPE();
            trueEdp += v->GetTruePE();
        }
        trueEdp /= 100;
        recoEdp /= 100;

        // To find the reconstructed range as we will do in the final detector, we create voxels out of 2D MPPC Hits.
        // Then we look for the largest cluster, since for the time being we will assume it corresponds to the true pion.
        // TODO: Not always assume that largest cluster is the pion.
        std::vector<ND280SFGDVoxelSet*> clusters = FindClusters(recoEvent->GetVoxels(),0);
        ND280SFGDVoxelSet* largestClust = new ND280SFGDVoxelSet();
        for(auto cluster:clusters){
            if(cluster->GetVoxels().size() > largestClust->GetVoxels().size()) largestClust = cluster;
        }

        // Fill the hits in the accumulated views:
        for(UInt_t i=0; i<inputEvent->GetHits().size(); i++){
            ND280SFGDHit* hit = inputEvent->GetHits()[i];
            if(hit->GetView() == 0){
                hXY->Fill(hit->GetX(),hit->GetY(),hit->GetPE());
            }
            if(hit->GetView() == 1){
                hXZ->Fill(hit->GetX(),hit->GetZ(),hit->GetPE());
            }
            if(hit->GetView() == 2){
                hYZ->Fill(hit->GetY(),hit->GetZ(),hit->GetPE());
            }
        }

        // We take the reconstructed range as the maximum euclidian distance between 2 voxels in the largest cluster in the event.
        double recoRng = largestClust->GetMaxEuclDist();
        double offset = 25.1;  // NOTE: the offset is because the offset in the range coming from particle gun starting outside the SFGD.
                               // This value was computed [previously] from the mean of a gaussian filled with true range - reco range.

        // true vs reco range 2D
        h_RngVSRng->Fill(1.*(truePion->GetRange())/10-offset,largestClust->GetMaxEuclDist());
        // range resolution
        h_RngRes->Fill( 100.*((1.*(truePion->GetRange())/10-offset)-recoRng)/(1.*(truePion->GetRange())/10-offset));
        // true momentum vs true range
        h_MomVSRng->Fill(1.*(truePion->GetRange())/10-offset,truePion->GetMomentum());
        // true momentum vs reco range
        h_MomVSRecoRng->Fill(recoRng,truePion->GetMomentum());
        // momentum resolution from range
        h_MomResRng->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum());

        if(truePion->GetMomentum() > 0 && truePion->GetMomentum() < 100){
            h_MomResRng_split[0]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum());
        }        
        else if(truePion->GetMomentum() > 100 && truePion->GetMomentum() < 200){
            h_MomResRng_split[1]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        }
        else if(truePion->GetMomentum() > 200 && truePion->GetMomentum() < 300){
            h_MomResRng_split[2]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        }
        else if(truePion->GetMomentum() > 300 && truePion->GetMomentum() < 400){
            h_MomResRng_split[3]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        }
        else if(truePion->GetMomentum() > 400 && truePion->GetMomentum() < 500){
            h_MomResRng_split[4]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        }


        // if(true){
        //     if(1.*(truePion->GetRange())/10-offset > 0 && 1.*(truePion->GetRange())/10-offset < 5){
        //         h_MomResRng_split[0]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum());
        //     }        
        //     else if(1.*(truePion->GetRange())/10-offset > 5 && 1.*(truePion->GetRange())/10-offset < 10){
        //         h_MomResRng_split[1]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        //     }
        //     else if(1.*(truePion->GetRange())/10-offset > 10 && 1.*(truePion->GetRange())/10-offset < 15){
        //         h_MomResRng_split[2]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        //     }
        //     else if(1.*(truePion->GetRange())/10-offset > 15 && 1.*(truePion->GetRange())/10-offset < 20){
        //         h_MomResRng_split[3]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        //     }
        //     else if(1.*(truePion->GetRange())/10-offset > 45 && 1.*(truePion->GetRange())/10-offset < 50){
        //         h_MomResRng_split[4]->Fill(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum()); 
        //         // if(  100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum() > 20){
        //         //     cout << "true: " << 1.*(truePion->GetRange())/10-offset << endl;
        //         //     cout << "reco: " << recoRng << endl << endl;
        //         //     recoEvent->DrawTrueTracks(false,"tracks");
        //         //     largestClust->DrawTrueTracks(true,"cluster");
        //         // }
        //     }
        // }


        // true vs reco Edep 2D
        h_EdpVSEdp->Fill(trueEdp,recoEdp);
        // edep resolution
        h_EdpRes->Fill((trueEdp-recoEdp/trueEdp));
        // true momentum vs true edep
        h_MomVSEdp->Fill(truePion->GetMomentum(),trueEdp);
        // true momentum vs reco range
        h_MomVSRecoEdp->Fill(truePion->GetMomentum(),recoEdp);
        // edep resolution
        h_EdpRes->Fill( 100.*(truePion->GetMomentum()-pByRng->Eval(recoRng))/truePion->GetMomentum());
        // momentum resolution from edep
        h_MomResEdp->Fill(100.*(truePion->GetMomentum()-pByEdp->Eval(recoEdp))/truePion->GetMomentum());

        // Take the initial momentum [MeV/c] as the maximum kinetic energy. [E^2 = p^2 + m^2. If we assume the mass is not detected E=p]
        double Ekin = inputEvent->GetTrueTracks()[0]->GetMomentum();
        // Take the reconstructedd energy as the sum of the true Edep in each voxel.
        double Edep = 0;
        for (auto v:inputEvent->GetVoxels()) {Edep+=v->GetTrueEdep();}
        // energy containment
        h_EnergyCont->Fill(Ekin,Edep);
    }
    if(!VERBOSE) cout << "\33[2K\rProcessed: " << "\033[1;32m100%\033[0m" << std::flush;

    // Lines to summaryze shower modes:
    cout << endl << endl << "\033[1;34m\t--- SUMMARY of EVENTS ---\033[0m" << endl;
    cout << "\tSelected:"   << selEvents << endl;
    cout << "\tOne Pi:        " << c_OnePion << "," << std::setprecision(2) << 100.*c_OnePion/selEvents << " %"<< endl;
    cout << "\tTwo Pi:        " << c_TwoPion << "," << std::setprecision(2) << 100.*c_TwoPion/selEvents << " %"<< endl;
    cout << "\tMlt Pi:        " << c_MltPion << "," << std::setprecision(2) << 100.*c_MltPion/selEvents << " %"<< endl;
    cout << "\t\033[34m--- One Pi Breakdown ---\033[0m" << endl;
    cout << "\tClean:              " << c_OnePionClean   << ", " << std::setprecision(2) << 100.*c_OnePionClean/c_OnePion   << " %"<< endl;
    cout << "\tWith Protons        " << c_OnePionPrFound << ","  << std::setprecision(2) << 100.*c_OnePionPrFound/c_OnePion << " %\t included in Clean sample." << endl;
    cout << "\twith Electrons:     " << c_OnePionElFound << ","  << std::setprecision(2) << 100.*c_OnePionElFound/c_OnePion << " %"<< endl;
    cout << "\twith Positrons:     " << c_OnePionPsFound << ","  << std::setprecision(2) << 100.*c_OnePionPsFound/c_OnePion << " %"<< endl;
    cout << "\twith Gammas:        " << c_OnePionPhFound << ","  << std::setprecision(2) << 100.*c_OnePionPhFound/c_OnePion << " %"<< endl;
    cout << "\tDecay [detected]:   " << c_OnePionMuFound << " [" << std::setprecision(2) << c_OnePionMuFoundAndDet << "], " << std::setprecision(3) << 100.*c_OnePionMuFoundAndDet/c_OnePionMuFound << " %"<< endl;
    cout << "\t\033[34m--- Two Pi Breakdown ---\033[0m" << endl;
    cout << "\tClean:              " << c_TwoPionClean   << ", " << std::setprecision(2) << 100.*c_TwoPionClean/c_TwoPion   << " %"<< endl;
    cout << "\tWith Protons        " << c_TwoPionPrFound << ","  << std::setprecision(2) << 100.*c_TwoPionPrFound/c_TwoPion << " %\t included in Clean sample." << endl;
    cout << "\twith Electrons:     " << c_TwoPionElFound << ","  << std::setprecision(2) << 100.*c_TwoPionElFound/c_TwoPion << " %"<< endl;
    cout << "\twith Positrons:     " << c_TwoPionPsFound << ","  << std::setprecision(2) << 100.*c_TwoPionPsFound/c_TwoPion << " %"<< endl;
    cout << "\twith Gammas:        " << c_TwoPionPhFound << ","  << std::setprecision(2) << 100.*c_TwoPionPhFound/c_TwoPion << " %"<< endl;
    cout << "\tDecay [detected]:   " << c_TwoPionMuFound << " [" << std::setprecision(2) << c_TwoPionMuFoundAndDet << "], " << std::setprecision(3) << 100.*c_TwoPionMuFoundAndDet/c_TwoPionMuFound << " %"<< endl;
    cout << "\t\033[34m--- Mlt Pi Breakdown ---\033[0m" << endl;
    cout << "\tClean:              " << c_MltPionClean   << ", " << std::setprecision(2) << 100.*c_MltPionClean/c_MltPion   << " %"<< endl;
    cout << "\tWith Protons        " << c_MltPionPrFound << ","  << std::setprecision(2) << 100.*c_MltPionPrFound/c_MltPion << " %\t included in Clean sample." << endl;
    cout << "\twith Electrons:     " << c_MltPionElFound << ","  << std::setprecision(2) << 100.*c_MltPionElFound/c_MltPion << " %"<< endl;
    cout << "\twith Positrons:     " << c_MltPionPsFound << ","  << std::setprecision(2) << 100.*c_MltPionPsFound/c_MltPion << " %"<< endl;
    cout << "\twith Gammas:        " << c_MltPionPhFound << ","  << std::setprecision(2) << 100.*c_MltPionPhFound/c_TwoPion << " %"<< endl;
    cout << "\tDecay [detected]:   " << c_MltPionMuFound << " [" << std::setprecision(2) << c_MltPionMuFoundAndDet << "], " << std::setprecision(3) << 100.*c_MltPionMuFoundAndDet/c_MltPionMuFound << " %"<< endl;
   
    // Fill histograms to compute most likely true momentum for a given reconstructed measurement.
    h_pVSr = (TH2F*) h_MomVSRecoRng->Clone();
    h_pVSe = (TH2F*) h_MomVSRecoEdp->Clone();
    // We add 2 because of the underflow and overflow bins.
    int NBinsX = h_pVSr->GetNbinsX()+2;
    int NBinsY = h_pVSr->GetNbinsY()+2;
    // keep only the median
    for(int i=0; i<NBinsX; i++){
        int maxBinCnt = 0;
        for(int j=0; j<NBinsY; j++) if (maxBinCnt < h_pVSr->GetBinContent(i,j)) maxBinCnt = h_pVSr->GetBinContent(i,j);
        if(maxBinCnt) for(int j=0; j<NBinsY; j++)  if (maxBinCnt > h_pVSr->GetBinContent(i,j))  h_pVSr->SetBinContent(i,j,0);
    }
    // keep only the median
    for(int i=0; i<NBinsX; i++){
        int maxBinCnt = 0;
        for(int j=0; j<NBinsY; j++) if (maxBinCnt < h_pVSe->GetBinContent(i,j)) maxBinCnt = h_pVSe->GetBinContent(i,j);
        if(maxBinCnt) for(int j=0; j<NBinsY; j++)  if (maxBinCnt > h_pVSe->GetBinContent(i,j))  h_pVSe->SetBinContent(i,j,0);
    }
    // remove low statistics
    for(int i=0; i<NBinsX*NBinsY; i++)
        if(h_pVSr->GetBinContent(i) < 3) h_pVSr->SetBinContent(i,0);
    // remove low statistics
    for(int i=0; i<NBinsX*NBinsY; i++)
        if(h_pVSe->GetBinContent(i) < 2) h_pVSe->SetBinContent(i,0);


    // Plot the results in different Canvas.
    TCanvas* r0 = new TCanvas("energy containment");
    r0->cd(1);
    h_EnergyCont->SetTitle("Energy containment");
    h_EnergyCont->GetXaxis()->SetTitle("Total Energy [MeV]");
    h_EnergyCont->GetYaxis()->SetTitle("Energy Deposit [MeV]");
    h_EnergyCont->Draw("COLZ");
    r0->Update();

    TCanvas* r1 = new TCanvas("range histograms");
    r1->Divide(2,2);
    r1->cd(1);
    h_RngVSRng->SetTitle("Range Migration");
    h_RngVSRng->GetXaxis()->SetTitle("True Range [cm]");
    h_RngVSRng->GetYaxis()->SetTitle("Reco Range [cm]");
    h_RngVSRng->Draw("COLZ");
    r1->cd(2);
    h_MomVSRng->SetTitle("True Range vs Mom");
    h_MomVSRng->GetXaxis()->SetTitle("True Momentum [cm]");
    h_MomVSRng->GetYaxis()->SetTitle("True Range [cm]");
    h_MomVSRng->Draw("COLZ");
    r1->cd(3);
    h_RngRes->SetTitle("Range Resolution");
    h_RngRes->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_RngRes->Draw("COLZ");
    r1->cd(4);
    h_MomResRng->SetTitle("Mom Resolution by Range");
    h_MomResRng->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResRng->Draw("COLZ");
    r1->Update();

    TCanvas* r2 = new TCanvas("edep histograms");
    r2->Divide(2,2);
    r2->cd(1);
    h_EdpVSEdp->SetTitle("Edep Migration");
    h_EdpVSEdp->GetXaxis()->SetTitle("True Edep");
    h_EdpVSEdp->GetYaxis()->SetTitle("Reco Edep");
    h_EdpVSEdp->Draw("COLZ");
    r2->cd(2);
    h_MomVSEdp->SetTitle("True Edep vs Mom");
    h_MomVSEdp->GetXaxis()->SetTitle("True Momentum [cm]");
    h_MomVSEdp->GetYaxis()->SetTitle("True Edep");
    h_MomVSEdp->Draw("COLZ");
    r2->cd(3);
    h_EdpRes->SetTitle("Edep Resolution");
    h_EdpRes->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_EdpRes->Draw("COLZ");
    r2->cd(4);
    h_MomResEdp->SetTitle("Mom Resolution by Edep");
    h_MomResEdp->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResEdp->Draw("COLZ");
    r2->Update();

    TCanvas* r3 = new TCanvas("accumulated hits");
    r3->Divide(2,2);
    r3->cd(1);
    hXY->Draw("COLZ");
    r3->cd(2);
    hXZ->Draw("COLZ");
    r3->cd(3);
    hYZ->Draw("COLZ");
    r3->Update();

    TCanvas* r4 = new TCanvas("median from true information");
    r4->Divide(2,2);
    r4->cd(1);
    h_pVSr->Draw("COLZ");
    TF1 *f1 = new TF1("f1","[0]+sqrt([1]*x)",0,200);
    f1->SetParameters(1,1);
    h_pVSr->Fit(f1,"Q");
    if(VERBOSE >0) cout << "Params: " << std::setprecision(8) << f1->GetParameter(0) << "," << f1->GetParameter(1) << endl;
    r4->cd(2);
    TF1 *f2 = new TF1("f1","[0]+sqrt([1]*x)",0,200);
    f2->SetParameters(1,1);
    h_pVSe->Draw("COLZ");
    if(USE_EDEP){
        h_pVSe->Fit(f2,"Q");
        cout << "Params: " << std::setprecision(8) << f2->GetParameter(0) << "," << f2->GetParameter(1) << endl;
    }

    r4->cd(3);
    h_MomVSRecoRng->SetTitle("Reco Range vs Mom");
    h_MomVSRecoRng->GetXaxis()->SetTitle("Reco Range [cm]");
    h_MomVSRecoRng->GetYaxis()->SetTitle("True Momentum [MeV]");
    h_MomVSRecoRng->Draw("COLZ");
    r4->cd(4);
    h_MomVSRecoEdp->SetTitle("Reco Edep vs Mom");
    h_MomVSRecoRng->GetXaxis()->SetTitle("Reco Edep");
    h_MomVSRecoRng->GetYaxis()->SetTitle("True Momentum [MeV]");
    h_MomVSRecoEdp->Draw("COLZ");
    r4->Update();

    TCanvas* r5 = new TCanvas("edep histograms");
    r5->Divide(3,2);
    r5->cd(1);
    h_MomResRng_split[0]->SetTitle("Mom Resolution by Range [0-100 MeV]");
    h_MomResRng_split[0]->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResRng_split[0]->Draw("HIST");
    r5->cd(2);
    h_MomResRng_split[1]->SetTitle("Mom Resolution by Range [100-200 MeV]");
    h_MomResRng_split[1]->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResRng_split[1]->Draw("HIST");
    r5->cd(3);
    h_MomResRng_split[2]->SetTitle("Mom Resolution by Range [200-300 MeV]");
    h_MomResRng_split[2]->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResRng_split[2]->Draw("HIST");
    r5->cd(4);
    h_MomResRng_split[3]->SetTitle("Mom Resolution by Range [300-400 MeV]");
    h_MomResRng_split[3]->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResRng_split[3]->Draw("HIST");
    r5->cd(5);
    h_MomResRng_split[4]->SetTitle("Mom Resolution by Range [400-500 MeV]");
    h_MomResRng_split[4]->GetXaxis()->SetTitle("[1-Reco/True] [%]");
    h_MomResRng_split[4]->Draw("HIST");
    r5->Update();


    if (STORE and FileOutput->IsOpen()){
        FileOutput->cd();
        h_pVSr->Write("",TObject::kOverwrite);
        FileOutput->Close();
        cout << endl << "\033[1;32mGenerated output: \033[0m" << fileOut.Data() << endl;
    }

    handleEndOfExecution();
}
