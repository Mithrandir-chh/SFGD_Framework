
/////// ALGORITHMS TO ORGANIZE AND FILL THE TRUE INFO OF THE RECONSTRUCTED EVENT //////////

//***********************************************************************************************
void FillTrueInformationInRecoEvent_version0(ND280SFGDEvent* inputEvent, ND280SFGDEvent* recoEvent, bool DEBUG){
//***********************************************************************************************

    // ___FILL_VOXELS_TRUE_INFORMATION___
    std::vector <ND280SFGDVoxel*> trueVoxels = inputEvent->GetVoxels();
    std::vector <ND280SFGDVoxel*> recoVoxels = recoEvent->GetVoxels();

    for(auto rv:recoVoxels) rv->SetTrueType(2); // by default we define reconstructed voxels as ghosts.

    int ID = 0;

    double minTime = 1E6;
    for (auto tv:trueVoxels){
        if(tv->GetTrueTime() < minTime) minTime = tv->GetTrueTime();
    }
    for (auto tv:trueVoxels) tv->SetTrueTime(tv->GetTrueTime()-minTime);
    for (auto tv:trueVoxels) if(tv->GetTrueTime()<100) tv->SetTime(0); else(tv->SetTime(1));

    for(auto rv:recoVoxels){
        rv->SetID(ID++);
        double trueEdep = 0;
        double truePE   = 0;
        for(auto tv:trueVoxels){
            if (rv->GetX() == tv->GetX() && rv->GetY() == tv->GetY() && rv->GetZ() == tv->GetZ()){
                if (tv->GetTrueType() == 1 && rv->GetTrueType()) rv->SetTrueType(1);
                else if (!tv->GetTrueType()) rv->SetTrueType(0); // if there is a deposit coming from a track, the voxels is defined as track. 
                trueEdep += tv->GetTrueEdep();
                truePE   += tv->GetTruePE();
                if(DEBUG) if(tv->GetTruePDGs().size() != 1 or tv->GetTrueTrackIDs().size() != 1 or tv->GetTrueParentIDs().size() != 1) {cout << "FillTrueInformationInRecoEvent_version0: size of true track vector information from input event should be 1." << endl; exit(1);}
                if(!rv->GetTrueTrackIDs().size()){
                    rv->AddTruePDG(tv->GetTruePDGs()[0]);
                    rv->AddTrueTrackID(tv->GetTrueTrackIDs()[0]);
                    rv->AddTrueParentID(tv->GetTrueParentIDs()[0]);
                }
                else{
                    bool found = false;
                    std::vector <int> trackIDs = rv->GetTrueTrackIDs();
                    for(uint it=0; it<trackIDs.size(); ++it) if(trackIDs[it] == tv->GetTrueTrackIDs()[0]) found = true;
                    if(!found){
                        rv->AddTruePDG(tv->GetTruePDGs()[0]);
                        rv->AddTrueTrackID(tv->GetTrueTrackIDs()[0]);
                        rv->AddTrueParentID(tv->GetTrueParentIDs()[0]);
                    }
                }
            }
        }
        rv->SetTrueEdep(trueEdep);
        rv->SetTruePE(truePE);
    }

    // ___SET_HITS_MULTIPLICITY___
    for(auto rv:recoVoxels) {rv->GetHits()[0]->SetMultiplicity(0); rv->GetHits()[1]->SetMultiplicity(0); rv->GetHits()[2]->SetMultiplicity(0);}
    for(auto rv:recoVoxels) {rv->GetHits()[0]->SetMultiplicity(rv->GetHits()[0]->GetMultiplicity()+1); rv->GetHits()[1]->SetMultiplicity(rv->GetHits()[1]->GetMultiplicity()+1); rv->GetHits()[2]->SetMultiplicity(rv->GetHits()[2]->GetMultiplicity()+1);}

    // ___FILL_VECTOR_OF_VOXEL_POINTERS_IN_MPPC_HITS___
    for(auto rv:recoVoxels) {rv->GetHits()[0]->AddVoxel(rv); rv->GetHits()[1]->AddVoxel(rv); rv->GetHits()[2]->AddVoxel(rv);}

    if(DEBUG){
        cout << "----- Summary of TrueType information -----" << endl;
        cout << "number of reconstructed voxels: " << recoVoxels.size() << endl;
        cout << "number of true voxel: " << trueVoxels.size() << endl;
        
        int voxelWithMissingHits = 0;
        for(auto tv:trueVoxels) if (!tv->GetHits()[0]->GetPE() or !tv->GetHits()[1]->GetPE() or !tv->GetHits()[2]->GetPE()) ++voxelWithMissingHits;
        cout << "number of true voxels with(out) 3 hits > 0 PE: " << trueVoxels.size() - voxelWithMissingHits << "," << voxelWithMissingHits << endl;

        int nTrack=0;
        int nXtalk=0;
        int nGhost=0;
        for(auto rv:recoVoxels){
            if (!rv->GetTrueType()) ++nTrack;
            else if (rv->GetTrueType() == 1) ++nXtalk;
            else if (rv->GetTrueType() == 2) ++nGhost;
        }
        cout << "number of true track, crosstalk, ghost: " << nTrack << "," << nXtalk << "," << nGhost << endl;
        cout << "number of voxels without true flag:     " << recoVoxels.size() - nTrack - nXtalk - nGhost << endl << endl;

        cout << "----- Summary of Track information in the Voxels -----" << endl;
        cout << "number of reconstructed voxels: " << recoVoxels.size() << endl;
        int nullTrueEdep = 0;
        int nullTruePE   = 0;
        int contributors [3] = {0,0,0}; // 0 contributors, 1 contributors, more than 1 contributor.
        int maxContributors = 0;
        for(auto rv:recoVoxels){
            if(!rv->GetTrueEdep()) ++nullTrueEdep;
            if(!rv->GetTruePE())   ++nullTruePE;
            int voxel_contributors = (int) rv->GetTrueTrackIDs().size();
            if(!voxel_contributors) ++contributors[0];
            else if(voxel_contributors == 1) ++contributors[1];
            else if(voxel_contributors > 1)  ++contributors[2];
            if (voxel_contributors > maxContributors)  maxContributors = voxel_contributors;
        }
        cout << "Voxels with NullTrueEdep,NullTruePE,Ghosts: " << nullTrueEdep << "," << nullTruePE << "," << nGhost << endl;
        cout << "Track contributors [0,1,>1]: " << contributors[0] << "," << contributors[1] << "," << contributors[2] << endl;
        cout << "Maximum number of contributors in a single voxel:  " << maxContributors << endl;
        cout << "----- Summary of Hits information -----" << endl;
        cout << "number of reconstructed voxels: " << recoVoxels.size() << endl;
        cout << "number of measured hits:        " << recoEvent->GetHits().size() << endl;
        printLightInformation(recoVoxels);
    }

    // ___FILL_VOXELS_IN_TRUE_TRACKS___
    std::vector <ND280SFGDTrack*> trueTracks = inputEvent->GetTrueTracks();
    // reco voxels are the ones we care about to be filled on the tracks according to trueTrackID information.
    for (auto v:recoVoxels) for (auto id:v->GetTrueTrackIDs()) for (auto t:trueTracks) if (t->GetTrackID() == id) t->AddVoxel(v);
    if(DEBUG){
        int missingIDs = 0;
        int voxelIDs [(int)recoVoxels.size()];
        for (int it=0; it<(int) recoVoxels.size(); ++it) voxelIDs[it] = 0;
        for (auto t:trueTracks) for (auto v:t->GetVoxels()) voxelIDs[v->GetID()]++;
        for (auto v:recoVoxels) if (v->GetTrueType() == 2) voxelIDs[v->GetID()]++;
        for (auto id:voxelIDs)  if (!id) ++missingIDs;
        cout << "missingIDs: " << missingIDs << endl; 
        if (missingIDs) exit(1);
    }

}


//***********************************************************************************************
void FillTrueInformationInRecoEvent_version1(ND280SFGDEvent* inputEvent, ND280SFGDEvent* recoEvent, bool DEBUG){
//***********************************************************************************************

    cout << "FillTrueInformationInRecoEvent version 1 is not implemented yet!" << endl;

    return;
}


//***********************************************************************************************
void FillTrueInformationInRecoEvent(ND280SFGDEvent* inputEvent, ND280SFGDEvent* recoEvent, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    vector <ND280SFGDHit*> outputHits;

    if( version == 0 ){
        // the goal of this algorithm is to fill the TrueType of the recononstructed voxels, comparing the list of true voxels
        // and the list of reconstructed voxels in HitsToVoxels.
        FillTrueInformationInRecoEvent_version0(inputEvent, recoEvent, DEBUG);
    }
    else if ( version == 1){
        FillTrueInformationInRecoEvent_version1(inputEvent, recoEvent, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in FillTrueInformationInRecoEvent function!" << endl;
        exit(1);
    }
    return;
}
