
// /////// ALGORITHMS TO MATCH THE RECO TRACK  //////////


//***********************************************************************************************
void MatchRecoToTrueTracks_verson0(vector <ND280SFGDVoxel*> voxelsInCluster, vector <ND280SFGDTrack*> recoTracks,  vector <ND280SFGDTrack*> trueTracks, bool DEBUG){
//***********************************************************************************************
    map<int,int> trueTrackVoxels;
    std::vector <int> trueTrkID;
    // find all true trackIDs in the cluster.
    for (auto v:voxelsInCluster) for (auto id:v->GetTrueTrackIDs()) if(!v->GetTrueType()) trueTrkID.push_back(id);
    for_each( trueTrkID.begin(), trueTrkID.end(), [&trueTrackVoxels]( int val ){ trueTrackVoxels[val]++; } );

    if(DEBUG){
        cout << "true track [id-#voxels]: " << endl;
        for( auto t:trueTrackVoxels ) cout << t.first << " " << t.second << endl;
    }

    // for all reco tracks, find if there is a reco track that succesfully matches a true track.
    for (auto rtrk:recoTracks){
        map<int,int> recoTrackVoxels;
        std::vector <int> recoTrkID;
        int trkTypeVxls = 0;
        // find all true trackIDs in the reco track (coming from voxels that are intersected by a track (type=0))
        for (auto v:rtrk->GetVoxels()){
            for (auto id:v->GetTrueTrackIDs()) if(!v->GetTrueType()) recoTrkID.push_back(id);
        }
        for_each( recoTrkID.begin(), recoTrkID.end(), [&recoTrackVoxels]( int val ){ recoTrackVoxels[val]++; } );
        double best_precision = 0;
        double best_recall    = 0;
        double best_f1score   = 0;
        int    best_trkID = -999;
        for( auto map_rtrk:recoTrackVoxels ) trkTypeVxls+=map_rtrk.second;
        for( auto map_rtrk:recoTrackVoxels ){
            rtrk->SetPrecision(1.*(map_rtrk.second)/trkTypeVxls);
            rtrk->SetRecall(1.*(map_rtrk.second)/trueTrackVoxels.find(map_rtrk.first)->second);
            if(rtrk->GetF1Score() > best_f1score){
                best_f1score   = rtrk->GetF1Score();
                best_precision = rtrk->GetPrecision();
                best_recall    = rtrk->GetRecall();
                best_trkID = map_rtrk.first;
            }
            if(DEBUG)cout <<  "[#recoVoxels,#trueVoxels,Precision,Recall,F1-score]: "  <<  map_rtrk.second << "," << trueTrackVoxels.find(map_rtrk.first)->second << ","<< rtrk->GetPrecision() << "," << rtrk->GetRecall() << "," << rtrk->GetF1Score() << endl << endl;
            if(rtrk->GetPrecision()>=0.5 && rtrk->GetRecall() > 0.5) {
                if(DEBUG)cout << "Is Reconstructed!!!!" << endl;
                rtrk->SetIsReco(true);
            }
        }
        rtrk->SetPrecision(best_precision);
        rtrk->SetRecall(best_recall);
        if(rtrk->IsReco()) for (auto ttrk:trueTracks) if (ttrk->GetTrackID() == best_trkID) rtrk->SetTrueTrack(ttrk);
    }
}


//***********************************************************************************************
void MatchRecoToTrueTracks(vector <ND280SFGDVoxel*> voxelsInCluster, vector <ND280SFGDTrack*> recoTracks,  vector <ND280SFGDTrack*> trueTracks,  Int_t version=0, bool DEBUG=false){
//***********************************************************************************************
    if( version == 0 ){
        MatchRecoToTrueTracks_verson0(voxelsInCluster, recoTracks, trueTracks, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in MatchRecoToTrueTracks function!" << endl;
        exit(1);
    }
}
