//
//
// César Jesús-Valls (cjesus@ifae.es)
// 10-10-2019
//
// template file to show how to construct a loop that iterates
// over MC/data events, applies or not some reconstruction and does some analysis
//
//

#define THIS_NAME trackFinding
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
bool    SHOW_TRUE      = false;
bool    SHOW_RECO      = false;
TString fileOut        = "~/Desktop/templateAna_Out.root";
TString fileIn         = "../../data/204_56_188/Reconstructed_SFGD_MC.root";
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

void trackFinding() {
    
    /// ------------START--------------

    bool DEBUG = false;
    bool DRAW  = true;
    //maxEvents = 200;

    int Ntotl [2] = {0,0};
    int Nmuon [2] = {0,0};
    int Nprot [2] = {0,0};
    int Npion [2] = {0,0};
    int Nelec [2] = {0,0};

    parseArguments();
    linkFilesToTTrees();

    std::vector<ND280SFGDHit*> mppc_hits;

    TH2F* h_migration = new TH2F("h_migration","h_migration",10,-0.5,9.5,10,-0.5,9.5);

    TH1F* h_RngRes [5];
    h_RngRes[0] = new TH1F("h_RngRes_All", "h_RngRes_All", 50,-50,50);
    h_RngRes[1] = new TH1F("h_RngRes_Muon","h_RngRes_Muon",50,-50,50);
    h_RngRes[2] = new TH1F("h_RngRes_Pion","h_RngRes_Pion",50,-50,50);
    h_RngRes[3] = new TH1F("h_RngRes_Prot","h_RngRes_Prot",50,-50,50);
    h_RngRes[4] = new TH1F("h_RngRes_Elec","h_RngRes_Elec",50,-50,50);

    // h_effRngRec && h_effRngTot   // 0: muon, 1: pion+-, 2: proton, 3: electron
    TEfficiency* h_effRng [5];
    h_effRng[0] = new TEfficiency("h_effRng_All", "h_effRng_All", 20,0,100);
    h_effRng[1] = new TEfficiency("h_effRng_Muon","h_effRng_Muon",20,0,100);
    h_effRng[2] = new TEfficiency("h_effRng_Pion","h_effRng_Pion",10,0,50);
    h_effRng[3] = new TEfficiency("h_effRng_Prot","h_effRng_Prot",10,0,50);
    h_effRng[4] = new TEfficiency("h_effRng_Elec","h_effRng_Elec",10,0,50);

    TEfficiency* h_effAng [5];
    h_effAng[0] = new TEfficiency("h_effAng_All", "h_effAng_All", 20,-1,1);
    h_effAng[1] = new TEfficiency("h_effAng_Muon","h_effAng_Muon",20,-1,1);
    h_effAng[2] = new TEfficiency("h_effAng_Pion","h_effAng_Pion",20,-1,1);
    h_effAng[3] = new TEfficiency("h_effAng_Prot","h_effAng_Prot",20,-1,1);
    h_effAng[4] = new TEfficiency("h_effAng_Elec","h_effAng_Elec",20,-1,1);

    TH1F* h_AngRes [5];
    h_AngRes[0] = new TH1F("h_AngRes_All", "h_AngRes_All", 50,-50,50);
    h_AngRes[1] = new TH1F("h_AngRes_Muon","h_AngRes_Muon",50,-50,50);
    h_AngRes[2] = new TH1F("h_AngRes_Pion","h_AngRes_Pion",50,-50,50);
    h_AngRes[3] = new TH1F("h_AngRes_Prot","h_AngRes_Prot",50,-50,50);
    h_AngRes[4] = new TH1F("h_AngRes_Elec","h_AngRes_Elec",50,-50,50);

    TEfficiency* h_effAng2track = new TEfficiency("h_effAng2track", "h_effAng2track", 20,0,1);


    cout << "Initial Event: " << evtIni << endl;
    cout << "Final Event: " << evtFin << endl;
    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;

        if( iev % 100 == 0) cout << 100*iev/nEvents << "%\n";

        mppc_hits = getEventMPPChits(iev); 

        // only take events with 2 tracks:
        if(inputEvent->GetTrueTracks().size() < 3) continue;

        recoEvent->SetHits(MergeHits(mppc_hits,0));
        recoEvent->SetVoxels(HitsToVoxels(recoEvent->GetHits(),0));
        convertCoordinates(inputEvent->GetVoxels());
        FillTrueInformationInRecoEvent(inputEvent,recoEvent,0);

        // if some of the tracks has less than 8 true voxels, continue.
/*        bool jumpToNextEvt = false;
        int vxl_cnt = 0; 
        for (auto t:inputEvent->GetTrueTracks()){
            auto id = t->GetTrackID();
            for (auto v:recoEvent->GetVoxels()){
                for (auto v_id:v->GetTrueTrackIDs()){
                    if (id == v_id) vxl_cnt++;
                }
            }
            if(vxl_cnt < 20) jumpToNextEvt = true;
            if(DEBUG) cout << "vxl_cnt: " << vxl_cnt << endl;
        }

        if(jumpToNextEvt) continue;*/

        if(DEBUG) for (auto t:inputEvent->GetTrueTracks())
            cout << "trackID,parentID,PDG,Range: " << t->GetTrackID() << "," << t->GetParentID() << "," << t->GetPDG() << endl;
        
        std::vector<ND280SFGDVoxelSet*> clusters = FindClusters(recoEvent->GetVoxels(),0);
        
        ND280SFGDVoxelSet* largestClust = new ND280SFGDVoxelSet();
        for(auto cluster:clusters){
            if(DEBUG) cout << "#voxels: " << cluster->GetVoxels().size() << endl;
            if(cluster->GetVoxels().size() > largestClust->GetVoxels().size()) largestClust = cluster;
        }

   
        // reject events with several clusters, they are potentially shower-like.
        double ratio = 1.*largestClust->GetVoxels().size()/recoEvent->GetVoxels().size();
        if (ratio < 0.99) continue;

        if(useNN){
            std::vector<ND280SFGDVoxel*> v_trk;
            if(DEBUG) cout << "original size: " << largestClust->GetVoxels().size() << endl;
            v_trk = ApplyFakeNN(largestClust->GetVoxels());
            largestClust->SetVoxels(v_trk);
            if(DEBUG) cout << "new size: " << largestClust->GetVoxels().size() << endl;
        }

        if(useClean){
            std::vector<ND280SFGDVoxel*> v_trk;
            if(DEBUG) cout << "original size: " << largestClust->GetVoxels().size() << endl;;
            v_trk = cleanVoxels(largestClust->GetVoxels());
            largestClust->SetVoxels(v_trk);
            if(DEBUG) cout << "new size: " << largestClust->GetVoxels().size() << endl;
        }

        std::vector<ND280SFGDVoxel*> sample;
        sample = largestClust->GetVoxels();
        std::vector<ND280SFGDTrack*> trackSegments = FindTrackSegments(sample,0);
        largestClust->SetVoxels(sample);




        if (DRAW) largestClust->DrawTrueTracks(false,"true");
        if (DRAW) largestClust->DrawRecoTracks(false,"segments");
        
        std::vector<ND280SFGDTrack*> tracks = MergeTrackSegments(trackSegments,0);

        MatchRecoToTrueTracks(largestClust->GetVoxels(),tracks,inputEvent->GetTrueTracks());

        for (auto trk:tracks) if (trk->IsReco()) trk->GetTrueTrack()->SetIsReco(true);

        int cnt_reco = 0;
        for (auto trk:inputEvent->GetTrueTracks()) {
            if(trk->IsReco())++cnt_reco;
            double trueRng = trk->GetMaxEuclDist();
            double trueAng = cos(trk->fMomVec.Eta());
            h_effRng[0]->Fill(trk->IsReco(),trueRng);
            Ntotl[trk->IsReco()]++;
            if(trk->GetPDG() == 13)                                { h_effRng[1]->Fill(trk->IsReco(),trueRng); h_effAng[1]->Fill(trk->IsReco(),trueAng); Nmuon[trk->IsReco()]++; }
            else if(trk->GetPDG() == 211 or trk->GetPDG() == -211) { h_effRng[2]->Fill(trk->IsReco(),trueRng); h_effAng[2]->Fill(trk->IsReco(),trueAng); Npion[trk->IsReco()]++; }
            else if(trk->GetPDG() == 2212)                         { h_effRng[3]->Fill(trk->IsReco(),trueRng); h_effAng[3]->Fill(trk->IsReco(),trueAng); Nprot[trk->IsReco()]++; }
            else if(trk->GetPDG() == 11)                           { h_effRng[4]->Fill(trk->IsReco(),trueRng); h_effAng[4]->Fill(trk->IsReco(),trueAng); Nelec[trk->IsReco()]++; }
        }

        
        for (auto trk:tracks){
            // If the track is reconstructed, look to its reco range, and to its true range, and fill plots depending on the true PDG.
            if (trk->IsReco()){
                ND280SFGDTrack* trueT = trk->GetTrueTrack();
                double trueRng = trueT->GetMaxEuclDist();
                if(trueT->GetPDG() == 13)                                { h_RngRes[1]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                else if(trueT->GetPDG() == 211 or trk->GetPDG() == -211) { h_RngRes[2]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                else if(trueT->GetPDG() == 2212)                         { h_RngRes[3]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                else if(trueT->GetPDG() == 11)                           { h_RngRes[4]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
            } 
        }

        //find number of tracks exiting the vertex:
        int NoutgoingTracks = 0;
        for (int v1 = 0; v1<(int) largestClust->GetVoxels().size(); v1++){
            std::vector<int> ids;
            ids.push_back(largestClust->GetVoxels()[v1]->GetRecoTrackIDs()[0]);
            for (int v2 = v1+1; v2<(int) largestClust->GetVoxels().size(); v2++){
                if(largestClust->GetVoxels()[v1]->DistToVoxel(largestClust->GetVoxels()[v2])<3) ids.push_back(largestClust->GetVoxels()[v2]->GetRecoTrackIDs()[0]);
            }
            sort( ids.begin(), ids.end() );
            ids.erase( unique( ids.begin(), ids.end() ), ids.end() );            
            if(NoutgoingTracks<ids.size()) NoutgoingTracks = ids.size();
        }

        cout << "TRUE TOPO: " << inputEvent->GetTrueTracks().size() << endl; 
        cout << "RECO TOPO: " << NoutgoingTracks << endl;
        h_migration->Fill(inputEvent->GetTrueTracks().size() ,NoutgoingTracks);
        // if (DRAW or SHOW_RECO){
        //     TCanvas* c_ang = new TCanvas("segment angles","");
        //     c_ang->Divide(3,1);
        //     TH1F* h_Tht = new TH1F("CosTheta","CosTheta",20,0,1);
        //     TH1F* h_Phi = new TH1F("CosPhi",  "CosPhi",  20,0,1);
        //     TH2F* h_Ang = new TH2F("2Dangle", "2Dangle", 20,0,1,20,0,1);
        //     for (auto t:tracks){
        //         if (t->GetVoxels().size()<5) continue;
        //         std::vector<double> fitResults =TrackFitter(t->GetVoxels(),0);
        //         if(fitResults.size()){
        //             XYZVector aux_vec(fitResults[4],fitResults[5],fitResults[6]);
        //             double cos_phi = abs(cos(aux_vec.Phi()));
        //             double cos_tht = abs(cos(aux_vec.Eta()));
        //             cout << "ang: " << cos_tht << "," << cos_phi << endl;
        //             h_Tht->Fill(cos_tht);
        //             h_Phi->Fill(cos_phi);
        //             h_Ang->Fill(cos_tht,cos_phi);
        //         }
        //          else { cout << "Fit failed!" << endl;}
        //     }
        //     cout << "Drawing angular information.\n";
        //     c_ang->cd(1);
        //     h_Tht->DrawCopy("HIST");
        //     c_ang->cd(2);
        //     h_Phi->DrawCopy("HIST");
        //     c_ang->cd(3);
        //     h_Ang->DrawCopy("COLZ");
        //     c_ang->Update();
        //     delete h_Tht;
        //     delete h_Phi;
        //     delete h_Ang;
        // }
        cout << "#Of Reco Tracks: " << tracks.size() << endl;
        cout << "True-Reco Tracks: " << inputEvent->GetTrueTracks().size() << "," << cnt_reco << endl;

        if(SHOW_RECO){
            if(inputEvent->GetTrueTracks().size() != cnt_reco){
                largestClust->DrawTrueTracks(false, "true");
                largestClust->DrawRecoTracks(true, "reco");  
            } 
        }

        if (DRAW) largestClust->DrawRecoTracks(true, "reco");
        //if (DRAW) recoEvent->DrawHitsAndVoxels(false, "hits");
        //if (DRAW) recoEvent->DrawVoxelsTruePE(true, "pe");
        selEvents++;
    }

    cout << "-- SUMARY --" << endl;
    cout << "Totl Tracks: " << Ntotl[0]+Ntotl[1] << "\tReconstructed: " <<  Ntotl[1] << "," << (1.*Ntotl[1])/(Ntotl[0]+Ntotl[1]) << " %" << endl;
    cout << "Muon Tracks: " << Nmuon[0]+Nmuon[1] << "\tReconstructed: " <<  Nmuon[1] << "," << (1.*Nmuon[1])/(Nmuon[0]+Nmuon[1]) << " %" << endl;
    cout << "Pion Tracks: " << Npion[0]+Npion[1] << "\tReconstructed: " <<  Npion[1] << "," << (1.*Npion[1])/(Npion[0]+Npion[1]) << " %" << endl;
    cout << "Prot Tracks: " << Nprot[0]+Nprot[1] << "\tReconstructed: " <<  Nprot[1] << "," << (1.*Nprot[1])/(Nprot[0]+Nprot[1]) << " %" << endl;
    cout << "Elec Tracks: " << Nelec[0]+Nelec[1] << "\tReconstructed: " <<  Nelec[1] << "," << (1.*Nelec[1])/(Nelec[0]+Nelec[1]) << " %" << endl;

    h_effRng[0]->SetLineColor(kBlack);
    h_effRng[1]->SetLineColor(kRed);
    h_effRng[2]->SetLineColor(kBlue);
    h_effRng[3]->SetLineColor(kGreen+1);
    h_effRng[4]->SetLineColor(kCyan);

    h_effAng[0]->SetLineColor(kBlack);
    h_effAng[1]->SetLineColor(kRed);
    h_effAng[2]->SetLineColor(kBlue);
    h_effAng[3]->SetLineColor(kGreen+1);
    h_effAng[4]->SetLineColor(kCyan);

    TCanvas c1;
    c1.Divide(4,2);
    c1.cd(1);
    h_effRng[1]->SetTitle("Muon Eff vs Range");
    h_effRng[1]->Draw();
    c1.cd(2);
    h_RngRes[1]->SetTitle("[true -reco]/true Muon Range");
    h_RngRes[1]->Draw("HIST");
    c1.cd(3);
    h_effAng[1]->SetTitle("Muon Eff vs Angle");
    h_effAng[1]->Draw();
    c1.cd(4);
    h_AngRes[1]->SetTitle("[true -reco]/true Muon Angle");
    h_AngRes[1]->Draw("HIST");
    c1.cd(5);
    h_effRng[3]->SetTitle("Prot Eff vs Range");
    h_effRng[3]->Draw();
    c1.cd(6);
    h_RngRes[3]->SetTitle("[true -reco]/true Prot Range");
    h_RngRes[3]->Draw("HIST");
    c1.cd(7);
    h_effAng[3]->SetTitle("Prot Eff vs Angle");
    h_effAng[3]->Draw();
    c1.cd(8);
    h_AngRes[3]->SetTitle("[true -reco]/true Prot Angle");
    h_AngRes[3]->Draw("HIST");
    c1.Update();

    TCanvas c2;
    c2.cd();
    h_migration->Draw("COLZ");
    c2.Update();
    c2.WaitPrimitive();

    cout << "Total Events: " << nEvents << ", selEvents: " << selEvents << endl;


    handleEndOfExecution();
}
