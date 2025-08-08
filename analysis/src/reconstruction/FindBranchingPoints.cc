
// /////// ALGORITHMS TO FIND BRANCHING POINTS //////////

//***********************************************************************************************
vector <ND280SFGDVoxel*> FindBranchingPoints_version0(vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG){
//***********************************************************************************************


    /// TODO: RENAME WITH INTUITIVE CLEAR NAMESS!!!!!!

    double MIN_DIST = 2.5;

    // vertex candidates:
    vector <ND280SFGDVoxel*> branchingVoxels;
    vector <ND280SFGDVoxel*> vtxCandidates;
    // for all the input voxels, find vtx candidates:
    for(auto v:inputVoxels){
        vector <ND280SFGDVoxel*> auxList;
        for (auto v_reset:inputVoxels) {v_reset->SetClusterID(-1); if(v_reset->DistToVoxel(v)>MIN_DIST) auxList.push_back(v_reset);}
            std::vector<ND280SFGDVoxelSet*> clustList =FindClusters(auxList,0,false,1,2);
            if(clustList.size()>2) {
                bool save = true;
                for(auto c:clustList) if(c->GetVoxels().size() < 4) save = false;
                if(save)branchingVoxels.push_back(v);
            }
    }

    for (auto b:branchingVoxels) b->SetClusterID(-1);
    std::vector<ND280SFGDVoxelSet*> bClust =FindClusters(branchingVoxels,0,false,1,2);


    vector <ND280SFGDVoxel*> singularBranchingVoxels;
    for(auto c:bClust){
        std::vector<double> X = c->GetX();
        std::vector<double> Y = c->GetY();
        std::vector<double> Z = c->GetZ();
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;
        for(auto x_it:X) sumX += x_it;
        for(auto y_it:Y) sumY += y_it;
        for(auto z_it:Z) sumZ += z_it;
        double ctrX = sumX/X.size();
        double ctrY = sumY/Y.size();
        double ctrZ = sumZ/Z.size();
        double minDist = 1E6;
        //cout << "ctr: " << ctrX << "," << ctrY << "," << ctrZ << endl;
        ND280SFGDVoxel* branchingPoint = nullptr;
        for(auto v:c->GetVoxels()){
            double distance = sqrt(pow(v->GetX()-ctrX,2)+pow(v->GetY()-ctrY,2)+pow(v->GetZ()-ctrZ,2));
            if (distance < minDist) {minDist = distance; branchingPoint = v;}
        }
        branchingPoint->SetTruePE(10000); 
        singularBranchingVoxels.push_back(branchingPoint);
    }

    ND280SFGDVoxelSet representation;
    std::vector<ND280SFGDVoxelSet*> tracks;
    std::vector<ND280SFGDVoxelSet*> vertices;
    if(singularBranchingVoxels.size()){
        vector <ND280SFGDVoxel*> brokenGraphs;
        vector <ND280SFGDVoxel*> vertex;
        for (auto v:inputVoxels) {
            v->SetClusterID(-1);
            bool found = false;
            for (auto s:singularBranchingVoxels) if((v->DistToVoxel(s)<=MIN_DIST)) found = true;
            if(found) vertex.push_back(v);
            else brokenGraphs.push_back(v); 
        }
        tracks    = FindClusters(brokenGraphs,0,false,1,2); 
        vertices  = FindClusters(vertex,0,false,1,2);  
        std::vector<ND280SFGDVoxel*> all_voxels;
        for( auto t:tracks){
            std::vector<ND280SFGDVoxel*> track_voxels = t->GetVoxels();
            all_voxels.insert(all_voxels.begin(),track_voxels.begin(),track_voxels.end());
        }
        for (auto vtx:vertices){
            std::vector<ND280SFGDVoxel*> vertex_voxels = vtx->GetVoxels();
            for(auto v:vertex_voxels){
                v->SetClusterID(v->GetClusterID()+tracks.size());
                if(v->GetTruePE()<5000){
                    double minDist = 1E6;
                    int IDcloser = -999;
                    ND280SFGDVoxelSet* tcloser = nullptr;
                    for(auto t:tracks) {
                        std::vector<ND280SFGDVoxel*> aux_vec;
                        aux_vec.push_back(v);
                        double dist = computeMinDistanceOfSets(aux_vec,t->GetVoxels());
                        if(dist<minDist){
                            minDist = dist;
                            IDcloser = t->GetVoxels()[0]->GetClusterID();
                            tcloser = t;
                        }
                    }
                    cout << "IDcloser: " << IDcloser << endl;
                    tcloser->AddVoxel(v);
                    v->SetClusterID(IDcloser);
                }
            all_voxels.push_back(v);
            }
        }
        representation.SetVoxels(all_voxels);
    } else {
        ND280SFGDVoxelSet* auxVxls = new ND280SFGDVoxelSet();
        auxVxls->SetVoxels(inputVoxels);
        tracks.push_back(auxVxls);
        for (auto i:inputVoxels) i->SetClusterID(-1);
        representation.SetVoxels(inputVoxels);
    }

    representation.DrawClusters(false,"afterBranching");

    int trackCnt = 0;
    for(auto t:tracks){
        std::vector<ND280SFGDTrack*> trackSegments = FindTrackSegments(t->GetVoxels(),0);
        std::vector<ND280SFGDTrack*> tracksFinal   = MergeTrackSegments(trackSegments,0);
        for(auto trk:tracksFinal) {for(auto vxl:trk->GetVoxels()) vxl->SetTruePE(trackCnt); trackCnt++;}
        cout << "trackCnt: " << trackCnt << endl;
        // debug print:
        // for(auto trk:tracksFinal){
        //     //cout << "TRK: " << trk->GetTrackID() << endl;
        //     /for (auto vxl:trk->GetVoxels()) cout << "VtrkID: " << vxl->GetTruePE() << endl;
        // }
    }

    for(auto vox:singularBranchingVoxels) {vox->SetTruePE(trackCnt); trackCnt++;} 

    for(auto vf:representation.GetVoxels()) cout << "vf: " << vf->GetTruePE() << endl;


    cout << "InVoxels: " << inputVoxels.size() << endl;
    cout << "bPoints:  " << singularBranchingVoxels.size() << endl;
    cout << "Repres:   " << representation.GetVoxels().size() << endl;
    representation.DrawVoxelsTruePE(false,"trackID");

        

        // MatchRecoToTrueTracks(largestClust->GetVoxels(),tracks,inputEvent->GetTrueTracks());

        // for (auto trk:tracks) if (trk->IsReco()) trk->GetTrueTrack()->SetIsReco(true);








    return singularBranchingVoxels;
}

//***********************************************************************************************
vector <ND280SFGDVoxel*> FindBranchingPoints(vector <ND280SFGDVoxel*> inputVoxels, Int_t version, bool DEBUG=false){
//***********************************************************************************************

    vector <ND280SFGDVoxel*> branchingVoxels;
    if( version == 0 ){
         branchingVoxels = FindBranchingPoints_version0(inputVoxels, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in FindBranchingPoints function!" << endl;
        exit(1);
    }
    return branchingVoxels;

}
