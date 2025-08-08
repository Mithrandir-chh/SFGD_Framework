
// /////// ALGORITHMS TO FIND CLUSTERS FOR A GIVEN COLLECTION OF VOXELS //////////

double computeThickness(std::vector <ND280SFGDVoxel*> voxels){
    double thickness = 0;
    double centroid[3] = {0,0,0};
    for (auto v:voxels){
        centroid[0]+=v->GetX(); 
        centroid[1]+=v->GetY(); 
        centroid[2]+=v->GetZ(); 
    }
    centroid[0]/=voxels.size();
    centroid[1]/=voxels.size();
    centroid[2]/=voxels.size();
    for (auto v:voxels) thickness += sqrt(pow(v->GetX()-centroid[0],2)+pow(v->GetY()-centroid[1],2)+pow(v->GetZ()-centroid[2],2));
    thickness/=voxels.size();
    return thickness;
}

//***********************************************************************************************
std::vector<ND280SFGDTrack*> MergeTrackSegments_version0(vector <ND280SFGDTrack*> inputSegments, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    
    int MIN_ANGLE = 12;
    int MIN_SIZE = 3;
    double REASONABLY_CLOSE = 1.1;
    double CLOSE = 1;
    double REASONABLE_ANGLE = 25;
    double MIN_ALIGNMENT  = 0.05;
    double MIN_UNIFORMITY = 0.05;
    double MIN_REGULARITY = 0.07;

    if(useNN or useClean){
        MIN_ANGLE = 10;
        MIN_SIZE = 2;
        REASONABLY_CLOSE = 0.9;
        CLOSE = 0.3;
        REASONABLE_ANGLE = 25;
        MIN_ALIGNMENT  = 0.04;
        MIN_UNIFORMITY = 0.04;
        MIN_REGULARITY = 0.05;
    }


    DEBUG = false;

    std::vector<ND280SFGDTrack*> outputSegments = inputSegments;
    bool merged = false;
    //cout << "# initial segments: " << inputSegments.size() << endl;
    for (int i=0; i<(int) inputSegments.size()-1; ++i){
        ND280SFGDTrack* s1 = inputSegments[i];
        ND280SFGDTrack* s2 = s1->GetNextSegment();
        if(s1==s2) continue;
        double costheta = 0;
        double theta    = 1E3;
        if (s1->GetFitParams().size() && s2->GetFitParams().size()){
            XYZVector u1(s1->GetFitParams()[4],s1->GetFitParams()[5],s1->GetFitParams()[6]);
            XYZVector u2(s2->GetFitParams()[4],s2->GetFitParams()[5],s2->GetFitParams()[6]);
            costheta = (u1.Unit()).Dot(u2.Unit());
            theta = acos(costheta)*360/(2*3.141592);
        }
        if (DEBUG) cout << "theta: " << theta << endl;
        if (true){//abs(theta) < 25 or abs(theta-180) <25){
            std::vector<ND280SFGDVoxel*> voxels; // voxels of both merged segments;
            std::vector<ND280SFGDVoxel*> v1 = s1->GetVoxels();
            std::vector<ND280SFGDVoxel*> v2 = s2->GetVoxels();
            voxels.insert(voxels.begin(),v1.begin(),v1.end());
            voxels.insert(voxels.begin(),v2.begin(),v2.end());
            std::vector<double> fitResults =TrackFitter(voxels,0);
            double tot_err = 1E6;
            double ave_dist = 1E6;
            if(fitResults.size()){
                tot_err = sqrt( pow(fitResults[10],2)+pow(fitResults[11],2)+pow(fitResults[12],2));
                ave_dist = fitResults[0];
                if(DEBUG) cout << "tot_err: " << tot_err << endl;
                if(DEBUG) cout << "aveDist: " << ave_dist << endl;
            }
            // the angle is quite small:
            bool small_angle = (abs(theta) < MIN_ANGLE or abs(theta-180) < MIN_ANGLE);
            // one segment is very short:
            
            ND280SFGDVoxelSet vs1;
            ND280SFGDVoxelSet vs2;
            
            vs1.SetVoxels(v1);
            vs2.SetVoxels(v2);
            double d1   =  vs1.GetMaxEuclDist();
            double d2   =  vs2.GetMaxEuclDist();

            bool small_size  =  d1 < MIN_SIZE or d2 < MIN_SIZE;

            // the segments are short and the angle is not too bad:
            bool reasonable  = ((abs(theta) < REASONABLE_ANGLE or abs(theta-180) < REASONABLE_ANGLE) && ave_dist < REASONABLY_CLOSE) or (ave_dist < CLOSE);

            // if the fit fails... use triangularity:
            bool aligned = false;
            bool uniform = false;
            bool regular = false;

            // if 
            // -> any of the initial 2 segments had no succesful fit
            // -> two segments are reasonable close and have reasonably small angle 
            // -> if the fit failed.
            // then, lets look into it from a geometrical point of view, measuring alignment, uniformity and regularity of the merged segment.
            if ( theta>360 or !fitResults.size() or ave_dist < REASONABLY_CLOSE){
                if (fitResults.size() and ave_dist > 3) continue;
                ND280SFGDVoxelSet vs;
                vs.SetVoxels(voxels);
                double dtot =  vs.GetMaxEuclDist();
                double triangularity = (d1+d2)/dtot;
                aligned = abs(1-triangularity) < MIN_ALIGNMENT  && d1 != dtot && d2 != dtot;
                double uniformity = (computeThickness(v1)/d1)/(computeThickness(v2)/d2);
                uniform = abs(1-uniformity) < MIN_UNIFORMITY;
                double regularity = sqrt(pow(abs(1-triangularity),2)+pow(abs(1-uniformity),2));
                regular = regularity < MIN_REGULARITY;
                if (DEBUG) cout << "triangularity: " << triangularity << endl;
                if (DEBUG) cout << "uniformity: " << uniformity << endl;
                if (DEBUG) cout << "regularity: " << regularity << endl;
            }

            if(d1 > 20 or d2 > 20) {aligned = false;}// reasonable = false;}
            if (theta<360) if (!reasonable) aligned = false;
            // avoid to merge a small segment to another that is far away.
            if (small_size) if (computeMinDistanceOfSets(v1,v2)>2) small_size = false;

            if(DEBUG) cout << "ave_dist: " << ave_dist << endl;
            if(DEBUG) cout << "d1,d2,v1,v2: " << d1 << "," << d2 << "," << v1.size() << "," << v2.size() << endl;
            if(DEBUG) cout << "bool [angle,size,reas,align,regu]: " << small_angle << "," << small_size <<  "," << reasonable << "," << aligned << "," << regular <<endl;
            if( small_angle or  small_size or reasonable or aligned or regular) {
                if(DEBUG) cout << "merged!" << endl;
                s1->SetVoxels(voxels);
                s1->SetFitParams(fitResults);
                s1->SetNextSegment(s2->GetNextSegment());
                std::vector<ND280SFGDTrack*> tmp_segments;
                for (auto S:inputSegments) if (S != s2) tmp_segments.push_back(S);
                outputSegments = tmp_segments;
                merged = true;
                break;
            }
        }
        if(merged) break;
    } 
    if(merged) outputSegments = MergeTrackSegments_version0(outputSegments,0);

    bool found;
    while(1){
        found = false;
        for(auto s1:outputSegments){
            double minDist = 1E6;
            ND280SFGDTrack* closest = nullptr;
            if (s1->GetVoxels().size()<4){
                for (auto s2:outputSegments){
                    if(s1 != s2){
                        double dist = computeMinDistanceOfSets(s1->GetVoxels(),s2->GetVoxels());
                        if(dist < minDist) {minDist = dist; closest = s2; found = true;}
                    }
                }
            }
            if(found) {
                for(auto v:closest->GetVoxels()) s1->AddVoxel(v); 
                std::remove(outputSegments.begin(), outputSegments.end(), closest);
                break;
            }
        }
        if(!found) break;
    }

    for(int i=0; i<(int)outputSegments.size(); ++i) for(auto v:outputSegments[i]->GetVoxels()) {v->ClearRecoTrackIDs();v->AddRecoTrackID(i);}
    if(DEBUG) cout << "\n";

    return outputSegments;
}


//***********************************************************************************************
std::vector<ND280SFGDTrack*> MergeTrackSegments_version1(vector <ND280SFGDTrack*> inputSegments, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    
    std::vector<ND280SFGDTrack*> outputSegments = inputSegments;
    bool merged = false;
    cout << "# initial segments: " << inputSegments.size() << endl;
    for (int i=0; i<(int) inputSegments.size(); ++i){
        ND280SFGDTrack* s1 = inputSegments[i];
        if (!s1->GetFitParams().size()) continue;
        XYZVector u1(s1->GetFitParams()[4],s1->GetFitParams()[5],s1->GetFitParams()[6]);
        for (int j=i+1; j<(int) inputSegments.size(); ++j){
            ND280SFGDTrack* s2 = inputSegments[j];
            if (!s2->GetFitParams().size()) continue;
            XYZVector u2(s2->GetFitParams()[4],s2->GetFitParams()[5],s2->GetFitParams()[6]);
            double costheta = (u1.Unit()).Dot(u2.Unit());
            double theta = acos(costheta)*360/(2*3.141592);
            if (DEBUG) cout << u1.X() << "," << u1.Y() << "," << u1.Z() << " || "  << u2.X() << "," << u2.Y() << "," << u2.Z() << "\n";
            if (DEBUG) cout << "costheta: " << costheta << ",theta: " << theta <<endl;
            if (theta < 10){
                std::vector<ND280SFGDVoxel*> voxels; // voxels of both merged segments;
                std::vector<ND280SFGDVoxel*> v1 = s1->GetVoxels();
                std::vector<ND280SFGDVoxel*> v2 = s2->GetVoxels();
                voxels.insert(voxels.begin(),v1.begin(),v1.end());
                voxels.insert(voxels.begin(),v2.begin(),v2.end());
                std::vector<double> fitResults =TrackFitter(voxels,0);
                double tot_err = 1E6;
                if(fitResults.size()){
                    tot_err = sqrt( pow(fitResults[10],2)+pow(fitResults[11],2)+pow(fitResults[12],2));
                    //cout << "tot_err: " << tot_err << endl;
                }
                if(tot_err<0.1){
                    s1->SetVoxels(voxels);
                    s1->SetFitParams(fitResults);
                    std::vector<ND280SFGDTrack*> tmp_segments;
                    for (auto S:inputSegments) if (S != s2) tmp_segments.push_back(S);
                    outputSegments = tmp_segments;
                    merged = true;
                    break;
                }
            }
        }
        if(merged) break;
    } 
    if(merged) {cout << "break!: " << outputSegments.size() << endl; outputSegments = MergeTrackSegments_version0(outputSegments,0);}
    


    for(int i=0; i<(int)outputSegments.size(); ++i) for(auto v:outputSegments[i]->GetVoxels()) {v->ClearRecoTrackIDs();v->AddRecoTrackID(i);}
    cout << "final number of segments: " << outputSegments.size() << endl;
    // for(int i=0; i<(int)outputSegments.size(); ++i) for(auto v:outputSegments[i]->GetVoxels()) v->GetRecoTrackIDs()[0];
    
    return outputSegments;
}


//***********************************************************************************************
std::vector<ND280SFGDTrack*> MergeTrackSegments(vector <ND280SFGDTrack*> inputSegments, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    
    std::vector<ND280SFGDTrack*> outputSegments;
    if( version == 0 ){
        outputSegments = MergeTrackSegments_version0(inputSegments, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in FindTrackCandidates function!" << endl;
        exit(1);
    }
    return outputSegments;
}
