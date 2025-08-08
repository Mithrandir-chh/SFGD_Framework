
// /////// ALGORITHMS TO RECONSTRUCT THE LIGHT IN THE CUBES FROM THE 2D HIT MEASUREMENTS //////////

//***********************************************************************************************
void ReconstructVoxelsPE_version0(vector <ND280SFGDHit*> inputHits, vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG){
//***********************************************************************************************
    
    int iterations   = 0;       // counter of number of iterations.
    double TotChange = 1E9;     // total change from one iteration to the next.
    double ratio     = 0;       // ratio between TotPEr / TotPEm. Measures how much charge has been already asigned.  
    double TotPEr    = 0;       // total number of reconstructed PE (from 2D hits).
    double TotPEm    = 0;       // total number of measured PE      (from 2D hits).
    double TotTruePE = 0;       // total number of true PE          (3D - only MC).
    double TotRecoPE = 0;       // total number of true PE          (3D - only MC).

    // map to count the amount of PE that has been assigned for each 2D hit. Replace it as a method in the hits. map.find is probably very slow.
    std::map<ND280SFGDHit*,double> PE_assigned;
    for(auto h:inputHits) PE_assigned[h] = 0;

    for (auto h:inputHits)              TotPEm     += h->GetPE();
    if  (IsMC) for (auto v:inputVoxels) TotTruePE  += v->GetTruePE();
    for (auto v:inputVoxels) v->SetRecoPE(0);

    while(iterations < 10000 && abs(1-ratio)> 0.01 && ratio < 1){
        TotChange = 0;
        for (auto vxl:inputVoxels){

            double dPE[3]  = {0,0,0};   // measures the contribution to the variation to the reconstructed PE in the voxel from each measurement.
            double Denom   = 0;         // denominator of the weight
            double Numer   = 0;         // numerator   of the weight
            double alpha   = 1E-2;      // growth factor.

            //cout << "truePE: " << 0.25*vxl->GetTruePE() << " [" << vxl->GetHits()[0]->GetMultiplicity() << "," << vxl->GetHits()[1]->GetMultiplicity() << "," << vxl->GetHits()[2]->GetMultiplicity() << "] " << "- [" << vxl->GetHits()[0]->GetPE() << "," << vxl->GetHits()[1]->GetPE() << "," << vxl->GetHits()[2]->GetPE() << "] " << "- [" << vxl->GetHits()[0]->GetPE()/att(vxl->DistToMPPC(0)) << "," << vxl->GetHits()[1]->GetPE()/att(vxl->DistToMPPC(1)) << "," << vxl->GetHits()[2]->GetPE()/att(vxl->DistToMPPC(2)) << "] " << endl;

            for (int x=0; x<3; ++x){

                Denom = 0;
                Numer = 0;

                // define indices [a,b] perpendicular to the current 2D hit fiber: 0: 1 & 2 || 1: 0 & 2 || 2: 1 & 0.
                int a = (x+1)%3;
                int b = (x+2)%3;

                // Estimation from perpendicular hits, considering all voxels in the perpendicular fibers.
                // For each hit travel to all the voxels that it intersects. For each of the neighbor voxels, look how much PE has been measured in the 2 perpendicular fibers
                // and how much of this PE have been already been asigned. Weight each PE_available =  PE_measured - PE_assigned by the number of neighbors in the fiber. Correct always the attenuation.
                for (auto v:vxl->GetHits()[x]->GetVoxels()){
                    Denom += (v->GetHits()[a]->GetPE()/att(v->DistToMPPC(a))-PE_assigned.find(v->GetHits()[a])->second/att(v->DistToMPPC(a)))/vxl->GetHits()[a]->GetVoxels().size() 
                           + (v->GetHits()[b]->GetPE()/att(v->DistToMPPC(b))-PE_assigned.find(v->GetHits()[b])->second/att(v->DistToMPPC(b)))/vxl->GetHits()[b]->GetVoxels().size();
                }
                
                // Estimation from perpendicular hits, considering only the voxel we are evaluating.
                Numer = (vxl->GetHits()[a]->GetPE()/att(vxl->DistToMPPC(a)) - PE_assigned.find(vxl->GetHits()[a])->second/att(vxl->DistToMPPC(a)))/vxl->GetHits()[a]->GetVoxels().size() 
                      + (vxl->GetHits()[b]->GetPE()/att(vxl->DistToMPPC(b)) - PE_assigned.find(vxl->GetHits()[b])->second/att(vxl->DistToMPPC(b)))/vxl->GetHits()[b]->GetVoxels().size();

                // if the only voxel contributing to the measured charge in the perpendicular views is itself, the Weight is 1. 
                // if the voxel has neighbors with larger (smaller) measurements on its perpendicular fibers, then the Weight is < 1 (>1).
                double Weight = Numer / Denom;

                // the variation of the PE in fiber x, is simply the product of the growth factor (which ensure convergence, and smoothens out ordering effects).
                // the available charge in the fiber (which ensures charge conservation)
                // and the weight (which ensures proportional distribution of charge among the voxels).
                dPE[x] = alpha*(vxl->GetHits()[x]->GetPE()/att(vxl->DistToMPPC(x))-PE_assigned.find(vxl->GetHits()[x])->second/att(vxl->DistToMPPC(x)))*Weight;
                //cout << "dPE[" << x << "]: " << dPE[x] << endl;
            }

            // compute the total variation as the average of the 3 hits.
            double dPE_TOT = (dPE[0]+dPE[1]+dPE[2])/3;

            // Impose a condition to ensure positive reconstructed charge on all voxels.
            if((vxl->GetRecoPE()+2*dPE_TOT)>0){
                // update the reconstructed PE in the voxel.
                vxl->SetRecoPE(vxl->GetRecoPE() + 2*dPE_TOT);
                // update the total variation in the iteration.
                TotChange += dPE_TOT;
                //update the assigned charge of each hit in the voxel.
                for(int x=0; x<3; ++x) PE_assigned[vxl->GetHits()[x]] += dPE_TOT/att(vxl->DistToMPPC(x));
            }
        }
        TotPEr = 0;
        for (auto h:inputHits) TotPEr += PE_assigned.find(h)->second;
        ratio = TotPEr/TotPEm;
        if(DEBUG && iterations%1000 == 0){
            TotRecoPE = 0;
            for(auto v:inputVoxels) TotRecoPE += 4*v->GetRecoPE();
            cout << "Charge measured: " << TotPEm << " || Charge given: " << TotPEr << " || ratio: " << ratio << endl;
            cout << "TotRecoPE: " << TotRecoPE << " || TotTruePE: " << TotTruePE << " || ratio: " <<  TotRecoPE/TotTruePE <<endl; 
            cout << "rr: " << ratio/(TotRecoPE/TotTruePE) << endl;
        }
        iterations++;
    }

    // correct 25% MPPC inefficiency:
    for(auto v:inputVoxels) v->SetRecoPE(4*v->GetRecoPE());

    if(DEBUG){
        TotPEr = 0;
        for (auto h:inputHits) TotPEr += PE_assigned.find(h)->second;
        TotRecoPE = 0;
        for(auto v:inputVoxels) TotRecoPE += v->GetRecoPE();
        cout << "SOLUTION: " << iterations << endl;
        cout << "Charge measured: " << TotPEm << " || Charge given: " << TotPEr << " || ratio: " << TotPEr/TotPEm << endl;
        cout << "TotRecoPE: " << TotRecoPE << " || TotTruePE: " << TotTruePE << " || ratio: " <<  TotRecoPE/TotTruePE <<endl; 
        cout << "rr: " << (TotPEr/TotPEm)/(TotRecoPE/TotTruePE) << endl;
        printLightInformation(inputVoxels);
    }
}


//***********************************************************************************************
void ReconstructVoxelsPE_version1(vector <ND280SFGDHit*> inputHits, vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG){
//***********************************************************************************************

    cout << "ReconstructVoxelPE version 1 is not implemented yet!" << endl;
}


//***********************************************************************************************
void ReconstructVoxelsPE(vector <ND280SFGDHit*> inputHits, vector <ND280SFGDVoxel*> inputVoxels, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    if( version == 0 ){
        // the goal of the algorithm is to estimate the original amount of light in each cube from the measured 2D hits.
        ReconstructVoxelsPE_version0(inputHits, inputVoxels, DEBUG);
    }
    else if ( version == 1){
        ReconstructVoxelsPE_version1(inputHits, inputVoxels, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in ReconstructVoxelPE function!" << endl;
        exit(1);
    }
}
