
// /////// ALGORITHMS TO FIND TRACK SEGMENTS FOR A GIVEN CLUSTER //////////


//***********************************************************************************************
void DFS(ND280SFGDVoxel* origen, bool DEBUG=false){
//***********************************************************************************************

    std::vector<ND280SFGDVoxel*> stack;
    ND280SFGDVoxel* node;
    stack.push_back(origen);
    int iterations = 0;
    double accum_distance = 0;
    while(stack.size() and iterations < 10000){
        ++iterations;
        node = stack.back();
        if(DEBUG) cout << "back to: " << node->GetID() << endl;
        node->SetIsVisited(true);
        ND280SFGDVoxel* deep_parent = node;
        accum_distance = node->GetDistance();
        if(DEBUG) cout << "node: "<< node->GetID() << ", neighbors: " << node->GetNeighbors().size() << endl;
        for (auto n:node->GetNeighbors())
            if(!n->IsVisited()){
                ND280SFGDVoxel* parent = node;
                node = n;
                node->SetIsVisited(true);
                stack.push_back(node);
                accum_distance+=parent->DistToVoxel(node);
                node->SetDistance(accum_distance);
                if(DEBUG) cout << "visiting: " << n->GetID() << endl;
            }
            else if(DEBUG) cout << "visited: " << n->GetID() << endl;
        if(node == deep_parent){
            if(DEBUG) cout << "no neighbors for: " << node->GetID() << endl; 
            stack.pop_back();
        }
    }

}


//***********************************************************************************************
void DefineGraph(std::vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG){
//***********************************************************************************************

    if(!inputVoxels.size()) return;

    // clusterID is used during graph generation. We store original cluster ID to restore it at the end of this function.
    int originalClusterID = inputVoxels[0]->GetClusterID();
    // prepare the voxels to be clustered with DBSCAN.
    for (auto v:inputVoxels) v->SetClusterID(-1);
    // search all voxels that can be connected by jumps of 1cm:
    std::vector<ND280SFGDVoxelSet*> clustList =FindClusters(inputVoxels,0,false,1,1);
    
    // for all voxels connected at 1cm, establish bi-direction links. [This is a non-directed graph.]
    for (ND280SFGDVoxelSet* cluster:clustList){
        std::vector<ND280SFGDVoxel*> nodes = cluster->GetVoxels();
            for(int i=0; i<(int) nodes.size(); ++i)
                for(int j=i+1; j<(int) nodes.size(); ++j)
                    if (nodes[i]->DistToVoxel(nodes[j]) == 1) {nodes[i]->AddNeighbor(nodes[j]); nodes[j]->AddNeighbor(nodes[i]);}
    }

    if(DEBUG) cout << "number of subclusters: " << clustList.size() << endl;

    int iterations = 0;
    double min_dist = 0;
    // take all nodes, and take 2 that, being from different subclusters are as close as possible. Link them as neighbors and merge the clusters.
    while (min_dist < 1E6 && iterations < (int)clustList.size() && clustList.size()>1){
        ++iterations;
        min_dist = 1E6;
        ND280SFGDVoxel* a;
        ND280SFGDVoxel* b;
        for(int i=0; i<(int) inputVoxels.size(); ++i)
            for(int j=i+1; j<(int) inputVoxels.size(); ++j){
                double dist = inputVoxels[i]->DistToVoxel(inputVoxels[j]);
                if(inputVoxels[i]->GetClusterID() == inputVoxels[j]->GetClusterID()) continue;
                if( dist > 1 && dist < min_dist){
                    a = inputVoxels[i];
                    b = inputVoxels[j];
                    min_dist = dist;
                }            
            }   
        a->AddNeighbor(b);
        b->AddNeighbor(a);
        int old_clust_id = b->GetClusterID();
        for(int i=0; i<(int) inputVoxels.size(); ++i) if(inputVoxels[i]->GetClusterID() == old_clust_id) inputVoxels[i]->SetClusterID(a->GetClusterID());
        if (DEBUG) {cout << "min_dist: " << min_dist << endl; cout << "clusters remaining: " << clustList.size() - iterations << endl;}
    }

    // resotore original clusterID for all voxels:
    for (auto v:inputVoxels) v->SetClusterID(originalClusterID);
}


//***********************************************************************************************
std::vector<ND280SFGDTrack*> FindTrackSegments_version0(vector <ND280SFGDVoxel*> nodes, bool DEBUG){
//***********************************************************************************************

    if (DEBUG) cout << "Nodes in the cluster: " << nodes.size() << endl;
    std::vector<ND280SFGDTrack*> segments;
    if(!nodes.size()) return segments;

    //int originalClusterID = nodes[0]->GetClusterID();
    TRandom3* rndm = new TRandom3(0);

    std::vector <ND280SFGDVoxel*> tNodes;  // tracked nodes
    std::vector <ND280SFGDVoxel*> uNodes;  // untracked nodes
    std::vector <ND280SFGDVoxelSet*> subGraphStack; // list of untracked nodes agrupataed in stacks;

    uNodes = nodes;

    int iterations = 0;
    while (tNodes.size() < nodes.size() && iterations < 10000){
        if(DEBUG )cout << "iteration: " << iterations << ", uNodes: " << uNodes.size() << endl;

        int MIN_SEGMENT_SIZE = 30; // this is local distance, not extrapolable to cm.
        int STEP_SIZE = 4; // this is local distance, not extrapolable to cm.

        if(useNN or useClean) MIN_SEGMENT_SIZE = 5;

        // if one subgraph has been completed, take the next one. [if there are not sub-graphs this will be always false].
        if(!uNodes.size()){
            uNodes = subGraphStack.back()->GetVoxels();
            subGraphStack.pop_back();
        }

        // we build a graph for the untracked nodes. We have to do this each time, in order to not jump to already tracked voxels.
        DefineGraph(uNodes,DEBUG);
        if(!iterations){
            // start from random origin.
            DFS(uNodes[rndm->Integer(uNodes.size())],false);
            // find one voxel close to an 'end-point' of the 3D graph. 
            double maxDist = 0;
            ND280SFGDVoxel* furthermostNode = nullptr;
            for(auto n:uNodes) if (maxDist<n->GetDistance()) {maxDist = n->GetDistance(); furthermostNode = n;}
            for(auto n:uNodes) {n->SetIsVisited(false); n->SetDistance(0);}
            // redefine interal distance with respect to the endpoint.
            if (furthermostNode) DFS(furthermostNode,false);
            // order nodes by distance to the endpoint.
        }
        else {for(auto n:uNodes) {n->SetIsVisited(false); n->SetDistance(0);}; DFS(uNodes[0],false);}
        sort(uNodes.begin(), uNodes.end(), compareDistance);

        // prepare the voxels to be clustered with DBSCAN.
        for (auto v:uNodes) v->SetClusterID(-1);
        // Find how many subgraphs there are on the graph formed by the untracked nodes.
        double max_separation = 3;
        if (useNN or useClean) max_separation = 2;
        std::vector<ND280SFGDVoxelSet*> clustList = FindClusters(uNodes,0,false,1,max_separation);
        if(DEBUG) cout << "clustList.size(): " << clustList.size() << endl;
        if(DEBUG) cout << "[u-t]Nodes.size(): " << uNodes.size() <<","<< tNodes.size() << "," << nodes.size() << endl;
        // if the graph can be split in sub-graphs, add the new sub-graphs to the stack of sub-graphs.
        if(clustList.size()>1){
            subGraphStack.insert(subGraphStack.begin(),clustList.begin(),clustList.end());
            uNodes = subGraphStack.back()->GetVoxels();
            subGraphStack.pop_back();
        }
        // select a group of nodes forming a linear segment
        double steps = 0;
        std::vector <ND280SFGDVoxel*> tmp_tNodes;
        std::vector <ND280SFGDVoxel*> tmp_uNodes;
        double dtot = -1;
        double d1   = -1;
        double d2   = -1;
        while(steps < 1000){
            ++steps;
            std::vector <ND280SFGDVoxel*> new_tmp_tNodes;
            std::vector <ND280SFGDVoxel*> new_tmp_uNodes;
            for (auto n:uNodes){
                if(tmp_tNodes.size()<MIN_SEGMENT_SIZE+STEP_SIZE*(steps-1)){
                    tmp_tNodes.push_back(n);
                }
                else{
                    tmp_uNodes.push_back(n);
                }
            }
            // take some nodes that are untracked and use them to try to extend the segment
            for (auto n:tmp_uNodes){
                if( (tmp_tNodes.size()+new_tmp_tNodes.size())< MIN_SEGMENT_SIZE+STEP_SIZE*(steps)){
                    new_tmp_tNodes.push_back(n);
                }
                else{
                    new_tmp_uNodes.push_back(n);
                }
            }
            if (DEBUG) cout << "sizes: " << tmp_tNodes.size() << "," << tmp_uNodes.size() << "," << new_tmp_tNodes.size() << "," << new_tmp_uNodes.size() << endl;
            
            std::vector<ND280SFGDVoxel*> voxels; 
            voxels.insert(voxels.begin(),new_tmp_tNodes.begin(),new_tmp_tNodes.end());
            voxels.insert(voxels.begin(),tmp_tNodes.begin(),tmp_tNodes.end());
            ND280SFGDVoxelSet vs_tot;
            vs_tot.SetVoxels(voxels);
            ND280SFGDVoxelSet vs_1;
            vs_1.SetVoxels(new_tmp_tNodes);
            ND280SFGDVoxelSet vs_2;
            vs_2.SetVoxels(tmp_tNodes);
            // Computing MaxEuclDist is computationally expensive. Here a simple way to makes things a bit easier. (This can be improved!)
            if (dtot < 0){
                dtot =  vs_tot.GetMaxEuclDist();
                d1   =  vs_1.GetMaxEuclDist();
                d2   =  vs_2.GetMaxEuclDist();
            }
            else{   
                //the previously formed segment had total length dtot, and now it is vs_2, so:
                d2 = dtot;
                // the new voxels are small, so this is easy to compute. Actually GetMaxEuclDist has to do N*(N-1) operations, where N is the voxel size!
                d1   =  vs_1.GetMaxEuclDist();
                // to compuate new dtot compare only the new voxels with the old ones. If the EuclDist is not larger then dtot is unchanged.
                double newMaxEuclDist = computeMaxDistanceOfSets(vs_1.GetVoxels(),vs_2.GetVoxels());
                if (newMaxEuclDist > dtot) dtot = newMaxEuclDist;
            }
            double alignment = 1.*(d1+d2)/dtot;
            // if the new voxels we are adding are probably not aligned, end the segment.
            if(DEBUG) cout << d1 << "," << d2 << "," << dtot << "," <<  new_tmp_tNodes.size() << "," << tmp_tNodes.size() << "," << voxels.size() << endl;
            if(DEBUG) cout << "alignment: " << alignment << endl;
            if ( !d1 or !d2 ) break;
            if (!(alignment < 1.05 && alignment > 0.95)) break;
            if (useNN or useClean) if (abs(alignment-1)>0.) break;
            tmp_tNodes.clear();
            tmp_uNodes.clear();
        }
        // update the list of track and untracked nodes for the newly formed segment.
        tNodes.insert(tNodes.begin(),tmp_tNodes.begin(),tmp_tNodes.end());
        uNodes.clear();
        for (auto u:tmp_uNodes) uNodes.push_back(u);
        std::vector<double> fitParameters = TrackFitter(tmp_tNodes,0);     
        //create a new track:
        ND280SFGDTrack* Segment = new ND280SFGDTrack();
        Segment->SetVoxels(tmp_tNodes);;
        Segment->SetFitParams(fitParameters);
        for (auto v:Segment->GetVoxels()) v->AddRecoTrackID(segments.size());
        segments.push_back(Segment);
        if(segments.size()>1) segments[segments.size()-2]->SetNextSegment(segments[segments.size()-1]);
        ++iterations;
    }
    
    //complete the cicle linking last segment to itself:
    segments[segments.size()-1]->SetNextSegment(segments[segments.size()-1]);

    if (DEBUG) cout << "#of trackCandidates: " << segments.size() << endl;
    if (DEBUG) cout << "#of tNodes: " << tNodes.size() << "," << nodes.size() << endl << endl;

    delete rndm;
    return segments;
}


//***********************************************************************************************
std::vector<ND280SFGDTrack*> FindTrackSegments(vector <ND280SFGDVoxel*> inputVoxels, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    
    std::vector<ND280SFGDTrack*> segments;
    if( version == 0 ){
        segments = FindTrackSegments_version0(inputVoxels, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in FindTrackCandidates function!" << endl;
        exit(1);
    }
    return segments;
}
