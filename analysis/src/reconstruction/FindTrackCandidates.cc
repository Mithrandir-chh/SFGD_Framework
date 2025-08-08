
// /////// ALGORITHMS TO FIND CLUSTERS FOR A GIVEN COLLECTION OF VOXELS //////////

bool compareDistance(ND280SFGDVoxel* v1, ND280SFGDVoxel* v2) 
{ 
    return (v1->GetDistance() < v2->GetDistance()); 
} 


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

    size_t iterations = 0;
    double min_dist = 0;
    // take all nodes, and take 2 that, being from different subclusters are as close as possible. Link them as neighbors and merge the clusters.
    while (min_dist < 1E6 && iterations < clustList.size() && clustList.size()>1){
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
std::vector<ND280SFGDTrack*> FindTrackCandidates_version0(vector <ND280SFGDVoxel*> nodes, bool DEBUG){
//***********************************************************************************************

    const int MIN_TRACK_LENGH = 10;

    if (DEBUG) cout << "Nodes in the cluster: " << nodes.size() << endl;
    std::vector<ND280SFGDTrack*> trackCandidates;
    DefineGraph(nodes,DEBUG);
    TRandom3* rndm = new TRandom3(0);
    DFS(nodes[rndm->Integer(nodes.size())],false);
    ND280SFGDVoxelSet* VS = new ND280SFGDVoxelSet();
    VS->SetVoxels(nodes);

    double maxDist = 0;
    ND280SFGDVoxel* furthermostNode = nullptr;
    for(auto n:nodes) if (maxDist<n->GetDistance()) {maxDist = n->GetDistance(); furthermostNode = n;}
    for(auto n:nodes) {n->SetIsVisited(false); n->SetDistance(0);}
    if (furthermostNode) DFS(furthermostNode,false);

    // order nodes by distance to the endpoint.
    sort(nodes.begin(), nodes.end(), compareDistance); 

    // select nodes until a minimum size is reached 

    for(int its=0; its<20; its++){
        std::vector <ND280SFGDVoxel*> newTrackVoxels;
        std::vector <ND280SFGDVoxel*> remainingVoxels;
        for (auto n:nodes){
            if(n->GetDistance()>=MIN_TRACK_LENGH+2*its) remainingVoxels.push_back(n);
            else{
                newTrackVoxels.push_back(n);
            }
        }
        // double euclMax = 0;
        // for(int a=0; a<(int)newTrackVoxels.size(); a++)
        //     for(int b=a+1; b<(int)newTrackVoxels.size(); b++)
        //         if(newTrackVoxels[a]->DistToVoxel(newTrackVoxels[b])>euclMax) euclMax = newTrackVoxels[a]->DistToVoxel(newTrackVoxels[b]);
        // cout << euclMax << "," << MIN_TRACK_LENGH+2*its << "," << euclMax/(MIN_TRACK_LENGH+2*its) << endl;
        std::vector<double> fitParameters = TrackFitter(newTrackVoxels,0);  
        if(fitParameters.size()){
            double tot_err = sqrt( pow(fitParameters[10],2)+pow(fitParameters[11],2)+pow(fitParameters[12],2));
            std::cout << "Fit errors: " << fitParameters[10] << "," << fitParameters[11] << "," << fitParameters[12] << "\t tot: " << tot_err << std::endl;
        }
    }




    delete VS;
    delete rndm;
    return trackCandidates;
}


//***********************************************************************************************
std::vector<ND280SFGDTrack*> FindTrackCandidates_version1(vector <ND280SFGDVoxel*> nodes, bool DEBUG){
//***********************************************************************************************

    const int MIN_TRACK_LENGH = 30;

    if (DEBUG) cout << "Nodes in the cluster: " << nodes.size() << endl;
    std::vector<ND280SFGDTrack*> trackCandidates;

    TRandom3* rndm = new TRandom3(0);

    std::vector <ND280SFGDVoxel*> tNodes;  // tracked nodes
    std::vector <ND280SFGDVoxel*> uNodes;  // untracked nodes

    uNodes = nodes;

    int iterations = 0;
    while (tNodes.size() < nodes.size() && iterations < 10000){
        ++iterations;
        cout << "iteration: " << iterations << ", uNodes: " << uNodes.size() << endl;
        // we build a graph for the untracked nodes.
        DefineGraph(uNodes,DEBUG);
        // start from random origin.
        DFS(uNodes[rndm->Integer(uNodes.size())],false);

        // find one voxel in a track end-point
        double maxDist = 0;
        ND280SFGDVoxel* furthermostNode = nullptr;
        for(auto n:uNodes) if (maxDist<n->GetDistance()) {maxDist = n->GetDistance(); furthermostNode = n;}
        for(auto n:uNodes) {n->SetIsVisited(false); n->SetDistance(0);}
        // redefine interal distance with respect to the endpoint.
        if (furthermostNode) DFS(furthermostNode,false);
        // order nodes by distance to the endpoint.
        sort(uNodes.begin(), uNodes.end(), compareDistance); 

        // select a group of nodes that are likely to be a track.
        double tot_err = 0;
        double best_tot_err = 1E3;
        double steps = 0;
        std::vector <ND280SFGDVoxel*> tmp_tNodes;
        std::vector <ND280SFGDVoxel*> tmp_uNodes;
        while(tot_err < 10 && steps < 1000){
            ++steps;
            tot_err = 1E6;
            for (auto n:uNodes){
                if(n->GetDistance() < MIN_TRACK_LENGH+2*steps) tmp_tNodes.push_back(n);
                else{
                    tmp_uNodes.push_back(n);
                }
            }
            std::vector<double> fitParameters = TrackFitter(tmp_tNodes,0);
            if(fitParameters.size()) tot_err = sqrt( pow(fitParameters[10],2)+pow(fitParameters[11],2)+pow(fitParameters[12],2));
            cout << "steps: " << steps << ", tNodes: " << tmp_tNodes.size() << endl;
            cout << "tot_err: " << tot_err << ", best_tot_err: " << best_tot_err << endl;
            if (tot_err == best_tot_err) break;
            if (best_tot_err > tot_err) best_tot_err = tot_err;
            if (tot_err > 2*best_tot_err) break;
            tmp_tNodes.clear();
            tmp_uNodes.clear();
        }
        cout << "steps-best_tot_err: " << steps << "," << best_tot_err << endl;

        uNodes = tmp_uNodes;
        for (auto t:tmp_tNodes) tNodes.push_back(t); // use insert!

        //create a new track:
        ND280SFGDTrack* Track = new ND280SFGDTrack();
        Track->SetVoxels(tmp_tNodes);
        for (auto v:Track->GetVoxels()) v->AddRecoTrackID(trackCandidates.size());
        trackCandidates.push_back(Track);
    }
    int cnt_vxl = 0;
    for (auto Track:trackCandidates) for (auto v:Track->GetVoxels()) cnt_vxl++;
    cout << "#of trackCandidates: " << trackCandidates.size() << endl;
    cout << "#of tNodes: " << cnt_vxl << "," << nodes.size() << endl;

    delete rndm;
    return trackCandidates;
}


//***********************************************************************************************
std::vector<ND280SFGDTrack*> FindTrackCandidates(vector <ND280SFGDVoxel*> inputVoxels, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    
    std::vector<ND280SFGDTrack*> trackCandidates;
    if( version == 0 ){
        trackCandidates = FindTrackCandidates_version0(inputVoxels, DEBUG);
    }
    if( version == 1 ){
        trackCandidates = FindTrackCandidates_version1(inputVoxels, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in FindTrackCandidates function!" << endl;
        exit(1);
    }
    return trackCandidates;
}
