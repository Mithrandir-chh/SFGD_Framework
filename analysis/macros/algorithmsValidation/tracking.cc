#define THIS_NAME tracking
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

#include <THStack.h>
#include <TLegend.h>

void tracking() {
    
    /// ------------START--------------

    parseArguments();
    linkFilesToTTrees();

    cout << "Program starts. " << endl;

    THStack *h_topVSNtrk = new THStack("h_topVSNtrk","");
    TH1F *h_Ntrk[4];
    // 0: CCQE || 1: 2p2h || 2: CC1PI || 3: DIS
    for(int i=0; i<4; ++i){
        TString h_name = "Ntrk_";
        h_name += i;
        h_Ntrk[i] = new TH1F(h_name.Data(),"",10,-0.5,9.5);
        h_Ntrk[i]->SetFillColor(i+1);
        h_Ntrk[i]->SetFillStyle(3001);
        h_topVSNtrk->Add(h_Ntrk[i]);
    }

    THStack *h_topVSmom = new THStack("h_topVSmom","");
    TH1F *h_Vmom[4];
    // 0: CCQE || 1: 2p2h || 2: CC1PI || 3: DIS
    for(int i=0; i<4; ++i){
        TString h_name = "Vmom_";
        h_name += i;
        h_Vmom[i] = new TH1F(h_name.Data(),"",20,0,2);
        h_Vmom[i]->SetFillColor(i+1);
        h_Vmom[i]->SetFillStyle(3001);
        h_topVSmom->Add(h_Vmom[i]);
    }

    TH2F* h_migration = new TH2F("h_migration","h_migration",10,-0.5,9.5,10,-0.5,9.5);
    TH1F* h_vtxRes[2];
    h_vtxRes[0]  = new TH1F("kink","", 100,0,100);
    h_vtxRes[1]  = new TH1F("brch","", 100,0,100);

    TH1*  h_nuMom     = new TH1F("h_nuMom", "", 100,0,10);

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

    int VERBOSE = 0;
    int topology [5] = {0,0,0,0,0}; // [CCQE, CC1PI, CCOthers, NoMuon, Others];

    std::vector<ND280SFGDHit*> mppc_hits;
    for (int iev=evtIni; iev<evtFin; iev++){
        if(iev == maxEvents-1 or selEvents >= maxSelEvents) break;

        mppc_hits = getEventMPPChits(iev);
        if(!VERBOSE) if(iev%10 ==0) cout << "\33[2K\rProcessed: \033[1;33m" << std::setprecision(2) << (100.*iev/(evtFin-evtIni)) << "%\t\033[0m" << std::flush;

        int Nmuon = 0;
        int Npion = 0;
        int Nprot = 0;
        int Nneut = 0;
        int Nexot = 0;

        for(auto t:inputEvent->GetTrueTracks()) {
            if(t->GetParentID() != 0) continue;
            if(t->GetPDG() == 13)   Nmuon++;
            else if(t->GetPDG() == 211 ) Npion++;
            else if(t->GetPDG() == -211) Npion++;
            else if(t->GetPDG() == 2212) Nprot++;
            else if(t->GetPDG() == 2112) Nneut++;
            else Nexot++;
        }

        bool CCQE     = (Nmuon && Nprot <= 1 && !(Npion)   && !(Nexot));
        bool CC2p2h   = (Nmuon && Nprot == 2 && !(Npion)   && !(Nexot));  // variables can not start by a number, so I use CC2p2h instead.
        bool CC1PI    = (Nmuon && Nprot <= 1 && Npion == 1 && !(Nexot));
        bool NoMuon   = !Nmuon;
        bool Exotic   = Nexot;
        bool CCOthers = (!CCQE && !CC2p2h && !CC1PI && !NoMuon && !Exotic);

        if (CCQE)     topology[0]++;
        if (CC1PI)    topology[1]++;
        if (CCOthers) topology[2]++;
        if (NoMuon)   topology[3]++;
        if (Exotic)   topology[4]++;

        if(!CC2p2h)   continue;

        recoEvent->SetHits(MergeHits(mppc_hits,0));
        recoEvent->SetVoxels(HitsToVoxels(recoEvent->GetHits(),0));
        convertCoordinates(inputEvent->GetVoxels());
        convertVtxCoordiantes(inputEvent->GetTrueVertex());
        FillTrueInformationInRecoEvent(inputEvent,recoEvent,0);


        double nuMom = inputEvent->GetNuMom()/1000;
        //cout << "neutrino momentum: " << nuMom << endl;
        h_nuMom->Fill(nuMom);
        //cout << "vertex   position: " << inputEvent->GetTrueVertex()->GetX() << "," << inputEvent->GetTrueVertex()->GetY() << "," << inputEvent->GetTrueVertex()->GetZ() << endl;
        
        if(useClean){
            Voxels clean = cleanVoxels(recoEvent->GetVoxels());
            recoEvent->SetVoxels(clean);
        }

        if(recoEvent->GetVoxels().size() < 10) continue;
        if(VERBOSE >0) for(auto t:inputEvent->GetTrueTracks())  cout << "trkID-prtID-PDG: " << t->GetTrackID() << "," << t->GetParentID() << "," << t->GetPDG() << endl;
        VxlSets listOfGraphs = FindClusters(recoEvent->GetVoxels(),0);
        int Ngraphs = listOfGraphs.size();

        ND280SFGDVoxelSet* largestClust = new ND280SFGDVoxelSet();
        for(auto cluster:listOfGraphs){;
            if(cluster->GetVoxels().size() > largestClust->GetVoxels().size()) largestClust = cluster;
        }

        if(Ngraphs>1) continue;

        Voxels  listOfBranchings     [Ngraphs];
        Voxels  listOfKinks          [Ngraphs];
        VxlSets listOfTrkCandidates  [Ngraphs];
        Tracks  listOfRecoTracks;

        // Loop over graphs:
        bool doReconstruction = true;
        if(CCOthers or Exotic or NoMuon) doReconstruction = false;
        if(doReconstruction){
            int NrecoTracks = 0;
            for (int i=0; i<Ngraphs; ++i){
                auto g = listOfGraphs[i];
                FindBranchingPoints(g,listOfBranchings[i]); 
                BreakGraph(g,listOfBranchings[i],listOfTrkCandidates[i]);
                if(listOfTrkCandidates[i].size() == 1){     //If te event is likely to be a CCQE:
                        KinkFinder(listOfKinks[i],listOfTrkCandidates[i]);     
                }
                if(listOfBranchings[i].size() == 1){
                    if(VERBOSE>1) cout << "\nBranchDistToVrx: " << listOfBranchings[i][0]->DistToVoxel(inputEvent->GetTrueVertex());
                    h_vtxRes[1]->Fill(listOfBranchings[i][0]->DistToVoxel(inputEvent->GetTrueVertex()));
                }
                else if (listOfKinks[i].size() == 1 && !listOfBranchings[i].size()){
                    if(VERBOSE>1)cout << "\nKinkDistToVrx: "   << listOfKinks[i][0]->DistToVoxel(inputEvent->GetTrueVertex());
                    h_vtxRes[0]->Fill(listOfKinks[i][0]->DistToVoxel(inputEvent->GetTrueVertex()));
                }
                CreateTracks(listOfBranchings[i],listOfTrkCandidates[i], listOfRecoTracks, NrecoTracks);        // asociates voxel and kink voxels to several track candidates.    
            }
            MatchRecoToTrueTracks(recoEvent->GetVoxels(),listOfRecoTracks,inputEvent->GetTrueTracks());
            for (auto trk:listOfRecoTracks) if (trk->IsReco()) trk->GetTrueTrack()->SetIsReco(true);
            int Nprimaries = 0;
            int Nreco      = 0;
            for (auto trk:inputEvent->GetTrueTracks()) {
                if(trk->GetParentID()) continue;
                if(trk->GetPDG()==2112) continue;
                int Nvxls = 0;
                for(auto v:trk->GetVoxels()) if(!v->GetTrueType()) Nvxls++;
                if(Nvxls<4) continue;
                Nprimaries++;
                if(VERBOSE>1){
                    cout << "#voxls: " << Nvxls << endl;
                    cout << "Is reco: " << trk->IsReco() << endl;
                    cout << "PDG: " << trk->GetPDG() << endl;
                }
                if(trk->IsReco()) Nreco++;
                double trueRng = trk->GetMaxEuclDist();
                double trueAng = cos(trk->fMomVec.Eta());
                h_effRng[0]->Fill(trk->IsReco(),trueRng);
                if(trk->GetPDG() == 13)                                { h_effRng[1]->Fill(trk->IsReco(),trueRng); h_effAng[1]->Fill(trk->IsReco(),trueAng);}
                else if(trk->GetPDG() == 211 or trk->GetPDG() == -211) { h_effRng[2]->Fill(trk->IsReco(),trueRng); h_effAng[2]->Fill(trk->IsReco(),trueAng);}
                else if(trk->GetPDG() == 2212)                         { h_effRng[3]->Fill(trk->IsReco(),trueRng); h_effAng[3]->Fill(trk->IsReco(),trueAng);}
                else if(trk->GetPDG() == 11)                           { h_effRng[4]->Fill(trk->IsReco(),trueRng); h_effAng[4]->Fill(trk->IsReco(),trueAng);}
            }
            for (auto trk:listOfRecoTracks){
                if (trk->IsReco()){
                    ND280SFGDTrack* trueT = trk->GetTrueTrack();
                    double trueRng = trueT->GetMaxEuclDist();
                    if(trueT->GetPDG() == 13)                                { h_RngRes[1]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                    else if(trueT->GetPDG() == 211 or trk->GetPDG() == -211) { h_RngRes[2]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                    else if(trueT->GetPDG() == 2212)                         { h_RngRes[3]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                    else if(trueT->GetPDG() == 11)                           { h_RngRes[4]->Fill( 100*(trueRng-trk->GetMaxEuclDist())/trueRng ) ;}
                } 
            }
            if(CCQE)        h_Ntrk[0]->Fill(Nprimaries);
            if(CC1PI)       h_Ntrk[1]->Fill(Nprimaries);
            if(CC2p2h)      h_Ntrk[2]->Fill(Nprimaries);
            if(CCOthers)    h_Ntrk[3]->Fill(Nprimaries);
            if(CCQE)        h_Vmom[0]->Fill(nuMom);
            if(CC1PI)       h_Vmom[1]->Fill(nuMom);
            if(CC2p2h)      h_Vmom[2]->Fill(nuMom);
            if(CCOthers)    h_Vmom[3]->Fill(nuMom);
            h_migration->Fill(Nprimaries,listOfRecoTracks.size());
        }

        if(true){
            // recoEvent->DrawBreakingPoints(false,"bPoints");
            // recoEvent->DrawVoxelsTruePE(false,"PE");

            recoEvent->DrawTrueTracks(false,"trueTracks");
            recoEvent->DrawRecoTracks(true,"recoTracks");   
        }
        selEvents++;
    }
    if(!VERBOSE) cout << "\33[2K\rProcessed: " << "\033[1;32m100%\033[0m" << std::flush;

    cout << endl <<  "Total-Selected Events: " << evtFin -evtIni << "," << selEvents << endl; 
    cout << "CCQE:     " << topology[0] << endl;
    cout << "CC1PI:    " << topology[1] << endl;
    cout << "CCOthers: " << topology[2] << endl;
    cout << "NoMuon:   " << topology[3] << endl;
    cout << "Others:   " << topology[4] << endl;

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
    c1.Divide(4,3);
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
    c1.cd(9);
    h_effRng[2]->SetTitle("Pion Eff vs Range");
    h_effRng[2]->Draw();
    c1.cd(10);
    h_RngRes[2]->SetTitle("[true -reco]/true Pion Range");
    h_RngRes[2]->Draw("HIST");
    c1.cd(11);
    h_effAng[2]->SetTitle("Pion Eff vs Angle");
    h_effAng[2]->Draw();
    c1.cd(12);
    h_AngRes[2]->SetTitle("[true -reco]/true Pion Angle");
    h_AngRes[2]->Draw("HIST");
    c1.Update();


    TCanvas c2;
    c2.cd();
    h_migration->Draw("COLZ");
    c2.Update();

    auto legend1 = new TLegend(0.65,0.75,0.85,0.85);
    legend1->SetBorderSize(0);
    legend1->AddEntry(h_vtxRes[0],"Kink","l");
    legend1->AddEntry(h_vtxRes[1],"Branching","l");

    auto legend2 = new TLegend(0.75,0.65,0.85,0.85);
    legend2->SetBorderSize(0);
    legend2->AddEntry(h_Vmom[0],"CCQE","f");
    legend2->AddEntry(h_Vmom[1],"CC1PI","f");
    legend2->AddEntry(h_Vmom[2],"2p2h","f");
    legend2->AddEntry(h_Vmom[3],"DIS","f");

    TCanvas c3;
    c3.Divide(2,2);
    c3.cd(1);
    gStyle->SetOptStat(0);
    h_vtxRes[0]->SetLineColor(kRed);
    h_vtxRes[1]->SetLineColor(kBlue);
    h_vtxRes[0]->Draw("HIST");
    h_vtxRes[1]->Draw("HIST same");
    legend1->Draw("same");
    c3.cd(3);
    h_topVSNtrk->Draw();
    legend2->Draw("same");
    c3.cd(4);
    h_topVSmom->Draw();
    legend2->Draw("same");
    c3.Update();
    c3.WaitPrimitive();

    handleEndOfExecution();
}
