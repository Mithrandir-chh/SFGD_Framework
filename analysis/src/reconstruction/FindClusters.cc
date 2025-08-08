
// /////// ALGORITHMS TO FIND CLUSTERS FOR A GIVEN COLLECTION OF VOXELS //////////

//***********************************************************************************************
std::vector<ND280SFGDVoxelSet*> FindClusters_version0(vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG, double MIN_NODES, double MAX_DISTANCE){
//***********************************************************************************************

    //___DBSCAN_ALGORITHM_CLUSTERING___

    // vector that will be filled with the output.
    std::vector<ND280SFGDVoxelSet*> recoClusters;

    // list of candidate nodes to be added to the cluster. 
    std::vector <ND280SFGDVoxel*> nodesToCheck;
    std::vector <ND280SFGDVoxel*> clusterVoxels;
    int clusterID = -1;  // ID of the first cluster. New cluster IDs will be clusterID++.
    int clustered = 0;   // number of voxels assigned to a cluster.
    while (clustered < (int) inputVoxels.size()){
        // if it is the beginning of the cluster, update clusterID, and add a voxel with unassigned cluster to this new cluster.
        // also add this voxel as a new candidate, to allow the algorithm to find new candidates from it.
        if (!nodesToCheck.size()){
            ND280SFGDVoxelSet* vSet = new ND280SFGDVoxelSet();
            recoClusters.push_back(vSet);
            clusterID++;
            auto it = std::find_if(inputVoxels.begin(), inputVoxels.end(), [](ND280SFGDVoxel* v) {return v->GetClusterID() < 0;});
            if (it != inputVoxels.end()) {
                nodesToCheck.push_back(*it);
                (*it)->SetClusterID(clusterID);
                clustered++;
            }
        }

        // find how many unclustered nodes there are within a sphere of radius MAX_DISTANCE, without including the origin node.
        // if the node is not already in nodesToCheck, add it to tmpNodes list. 
        std::vector <ND280SFGDVoxel*> tmpNodes;
        while(nodesToCheck.size() > 0){
            int closeNodes = 0;
            ND280SFGDVoxel* origin = nodesToCheck[0];
            for (auto neighbor:inputVoxels){
                if(origin->GetID() == neighbor->GetID()) continue;
                if(neighbor->DistToVoxel(origin) <= MAX_DISTANCE){
                    closeNodes++;
                    if(neighbor->GetClusterID() < 0){
                        bool found = std::any_of(nodesToCheck.begin(), nodesToCheck.end(),
                                    [&](ND280SFGDVoxel* v) { return v->GetID() == neighbor->GetID();}) ? true : false;
                        if(!found) tmpNodes.push_back(neighbor);
                    }
                }
            }
            // if the number of nodes within the spehere of radius MAX_DISTANCE is >= that MIN_NODES, then add tmpNodes to nodesToCheck.
            // add the origin node to the current voxel and increase the number of clustered voxels by 1.
            if(closeNodes >= MIN_NODES){
                nodesToCheck.insert(nodesToCheck.end(), tmpNodes.begin(), tmpNodes.end());
                if(origin->GetClusterID()<0){
                    clustered++;
                    origin->SetClusterID(clusterID);
                }
            }
            // clean tmpNodes to start from scratch for the next candidate node.
            tmpNodes.clear();
            // remove the origin node from nodeTocheck.
            nodesToCheck.erase(nodesToCheck.begin());
        }
    }
    for(auto v:inputVoxels) recoClusters[v->GetClusterID()]->AddVoxel(v);
    if (DEBUG){
        cout << "number of clusters: " << recoClusters.size() << endl;
        for(auto c:recoClusters) cout << "number of voxels: " << c->GetVoxels().size() << endl;
    }
    return recoClusters;
}

//***********************************************************************************************
std::vector<ND280SFGDVoxelSet*> FindClusters(vector <ND280SFGDVoxel*> inputVoxels, Int_t version, bool DEBUG=false, double MIN_NODES=1, double MAX_DISTANCE=2){
//***********************************************************************************************
    
    std::vector<ND280SFGDVoxelSet*> recoClusters;
    if( version == 0 ){
        
        recoClusters = FindClusters_version0(inputVoxels, DEBUG, MIN_NODES, MAX_DISTANCE);
    }
    else{
        cerr << "Version" << version << "does not exist in FindClusters function!" << endl;
        exit(1);
    }
    return recoClusters;
}
