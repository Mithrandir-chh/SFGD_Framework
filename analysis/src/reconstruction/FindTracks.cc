
double MIN_DIST = 2.5; 

typedef std::vector<ND280SFGDVoxelSet*> VxlSets;
typedef std::vector<ND280SFGDHit*>      MPPCHits;
typedef std::vector<ND280SFGDVoxel*>    Voxels;
typedef std::vector<ND280SFGDTrack*>    Tracks;

void FindBranchingPoints(ND280SFGDVoxelSet* graph, std::vector<ND280SFGDVoxel*>  &bPoints){

    // branching candidates:
    vector <ND280SFGDVoxel*> bCandidates;
    // for all the input voxels, find branching candidates:
    for(auto v:graph->GetVoxels()){
        vector <ND280SFGDVoxel*> auxList;
        for (auto x:graph->GetVoxels()) {x->SetClusterID(-1); if(x->DistToVoxel(v)>MIN_DIST) auxList.push_back(x);}
            std::vector<ND280SFGDVoxelSet*> clustList =FindClusters(auxList,0,false,1,2);
            if(clustList.size()>2) {
                bool save = true;
                for(auto c:clustList) if(c->GetVoxels().size() < 5) save = false;
                if(save) bCandidates.push_back(v);
            }
    }

    for (auto b:bCandidates) b->SetClusterID(-1);
    std::vector<ND280SFGDVoxelSet*> bClust =FindClusters(bCandidates,0,false,1,2);

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
        ND280SFGDVoxel* bPoint = nullptr;
        for(auto v:c->GetVoxels()){
            double distance = sqrt(pow(v->GetX()-ctrX,2)+pow(v->GetY()-ctrY,2)+pow(v->GetZ()-ctrZ,2));
            if (distance < minDist) {minDist = distance; bPoint = v;}
        }
        bPoint->SetIsBranching(true);
        bPoints.push_back(bPoint);
    }
    return;
}

void KinkFinder(std::vector<ND280SFGDVoxel*> &kinks, VxlSets &trks){


    VxlSets newListOfTrks;
    for(auto trk:trks){
        Voxels vxls = trk->GetVoxels();
        if(vxls.size() < 6) {newListOfTrks.push_back(trk);continue;}
        ND280SFGDVoxel* s1 = nullptr;
        ND280SFGDVoxel* s2 = nullptr;
        ND280SFGDVoxel* s3 = nullptr;
        double maxDist = 0;
        for (auto v1:vxls){
            for (auto v2:vxls){
                for (auto v3:vxls){
                    if( v1 == v2 or v1 == v3 or v2 == v3) continue;
                    double totDist = v1->DistToVoxel(v2)+v1->DistToVoxel(v3)+v2->DistToVoxel(v3);
                    if(totDist>maxDist){
                        maxDist = totDist;
                        s1 = v1;
                        s2 = v2;
                        s3 = v3;
                    }
                }
            }
        }
        double d1 = s1->DistToVoxel(s2);
        double d2 = s1->DistToVoxel(s3);
        double d3 = s2->DistToVoxel(s3);
        double minDist = d1;
        if(d2 < minDist) minDist = d2;
        if(d3 < minDist) minDist = d3;
        // compute lines between the 3 most separated points
        XYZVector x1(s1->GetX(), s1->GetY(), s1->GetZ()); 
        XYZVector x2(s2->GetX(), s2->GetY(), s2->GetZ());  
        XYZVector x3(s3->GetX(), s3->GetY(), s3->GetZ());  
        XYZVector u1 = (x2-x3).Unit();
        XYZVector u2 = (x1-x3).Unit();
        XYZVector u3 = (x1-x2).Unit();
        double cnt1 = 0;
        double cnt2 = 0;
        double cnt3 = 0;
        // compute distance to each line. Count how many voxels have a line as the closest.
        for(auto v:vxls){
            XYZVector xp(v->GetX(),v->GetY(),v->GetZ()); 
            double sep1 = sqrt(((xp-x2).Cross(u1)).Mag2());
            double sep2 = sqrt(((xp-x1).Cross(u2)).Mag2());
            double sep3 = sqrt(((xp-x1).Cross(u3)).Mag2());
            if      (sep1 <= sep2 && sep1 <= sep3) cnt1++;
            else if (sep2 <= sep1 && sep2 <= sep3) cnt2++;
            else if (sep3 <= sep2 && sep3 <= sep1) cnt3++;
        }
        // normalize the counts by the distance.
        cnt1 /= d3;
        cnt2 /= d2;
        cnt3 /= d1;
        ND280SFGDVoxel* kink = nullptr;
        if(cnt1 < cnt2 && cnt1 < cnt3) kink = s1;
        if(cnt2 < cnt1 && cnt2 < cnt3) kink = s2;
        if(cnt3 < cnt1 && cnt3 < cnt2) kink = s3;
        if(minDist > 3){
            ND280SFGDVoxelSet* vs1 = new ND280SFGDVoxelSet();
            ND280SFGDVoxelSet* vs2 = new ND280SFGDVoxelSet();
            Voxels vlist1;
            Voxels vlist2;
            kink->SetIsKink(true);
            kinks.push_back(kink);
            XYZVector x0(kink->GetX(), kink->GetY(), kink->GetZ());
            XYZVector l1;
            XYZVector l2; 
            if(kink == s1){ 
                l1 = u2;
                l2 = u3;
            }
            else if(kink == s2){
                l1 = u1;
                l2 = u3;
            }
            else{
                l1 = u1;
                l2 = u2;
            }
            for(auto v:vxls){
                XYZVector xv(v->GetX(),v->GetY(),v->GetZ()); 
                if( sqrt(((xv-x0).Cross(l1)).Mag2()) < sqrt(((xv-x0).Cross(l2)).Mag2())) vlist1.push_back(v);
                else vlist2.push_back(v);
            }
            vs1->SetVoxels(vlist1);
            vs2->SetVoxels(vlist2);
            newListOfTrks.push_back(vs1);
            newListOfTrks.push_back(vs2);
        }
        else{
            newListOfTrks.push_back(trk);
        }
    }
    trks = newListOfTrks;
}

void BreakGraph(ND280SFGDVoxelSet* graph, std::vector<ND280SFGDVoxel*>  breakPoints, std::vector<ND280SFGDVoxelSet*> &subGraphs){

    // To break the graph we remove a group of voxels around the breaking points [branching or kink], then apply DBSCAN to form subGraphs.
    // To finish, the voxels that were removed are associated to the subGraphs, by setting them to belong to the closes subGraph.

    // Classify volxes in close or far from some branching point.
    Voxels farFromBreaking;
    Voxels closeToBreaking;
    for(auto v:graph->GetVoxels()){
        // Reset the cluster ID which will be used by DBSCAN
        v->SetClusterID(-1);
        bool found = false;
        for (auto b:breakPoints){
            if(v->DistToVoxel(b) <= MIN_DIST) found = true;
        }
        if(found) closeToBreaking.push_back(v); 
        else      farFromBreaking.push_back(v);
    }
    // Form the subclusters for the disconnected voxels
    subGraphs = FindClusters(farFromBreaking,0,false,1,2);

    // We need to first sort the nodes in closeToBreaking by distance to subGraph
    for(auto v:closeToBreaking){
        // Not associate the breaking point to any subGraph in particular. 
        // Later [in a different function] the breaking point is associated to all sourrounding tracks.
        if(v->IsBranching()) continue;
        double minDist = 1E6;
        ND280SFGDVoxelSet* closestSG = nullptr;
        for(auto s:subGraphs){
            Voxels aux_v;
            aux_v.push_back(v);
            double dist = computeMinDistanceOfSets(aux_v,s->GetVoxels());
            if(dist < minDist) {minDist = dist; closestSG = s;}
        }
        //minDistances.push_back(minDist);
        v->SetDistance(minDist);
        if(!closestSG) {cout << "BreakGraph ERROR: closestSG not found!" << endl; exit(1);}
    }

    // sort closeToBreaking by minDistances...
    sort(closeToBreaking.begin(), closeToBreaking.end(), compareDistance);

    // use the sorted vector to distribute the unasigned voxels.s
    for(auto v:closeToBreaking){
        // Not associate the breaking point to any subGraph in particular. 
        // Later [in a different function] the breaking point is associated to all sourrounding tracks.
        if(v->IsBranching() or v->IsKink()) continue;
        double minDist = 1E6;
        ND280SFGDVoxelSet* closestSG = nullptr;
        for(auto s:subGraphs){
            Voxels aux_v;
            aux_v.push_back(v);
            double dist = computeMinDistanceOfSets(aux_v,s->GetVoxels());
            if(dist < minDist) {minDist = dist; closestSG = s;}
        }
        if(!closestSG) {cout << "BreakGraph ERROR: closestSG not found!" << endl; exit(1);}
        v->SetClusterID(closestSG->GetVoxels()[0]->GetClusterID());
        closestSG->AddVoxel(v);
    }

    for (auto b:breakPoints) b->SetClusterID(-1);
}

void CreateTracks(std::vector<ND280SFGDVoxel*>  breakPoints,  std::vector<ND280SFGDVoxelSet*> &trackCandidates, Tracks &recoTracks, int &NrecoTracks){

    // First include breaking points into track candidates
    for(auto b:breakPoints){
        Voxels aux_v;
        aux_v.push_back(b);
        double minDist = 1E6;
        for(auto trk:trackCandidates){
            if(computeMinDistanceOfSets(aux_v,trk->GetVoxels()) <= 2) {/*cout << b->GetID() << " -> " << trk->GetVoxels()[0]->GetClusterID() << endl;*/ trk->AddVoxel(b);}
            if(minDist>computeMinDistanceOfSets(aux_v,trk->GetVoxels())) minDist = computeMinDistanceOfSets(aux_v,trk->GetVoxels());
        }
    }
    // Then create tracks
    for(auto trk:trackCandidates){
        ND280SFGDTrack* newTrack = new ND280SFGDTrack();
        newTrack->SetVoxels(trk->GetVoxels());
        for(auto v:newTrack->GetVoxels()) v->AddRecoTrackID(NrecoTracks);
        recoTracks.push_back(newTrack);
        if(trk->GetVoxels().size()>4){
            newTrack->SetTrackID(NrecoTracks);
            NrecoTracks++;
        }
        else{
            newTrack->SetTrackID(-1);
        }
    }
}
