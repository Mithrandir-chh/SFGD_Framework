
/////// ALGORITHMS TO CREATE VOXELS FROM THE 2D HIT MEASUREMENTS //////////

//***********************************************************************************************
vector <ND280SFGDVoxel*> HitsToVoxels_version0(vector <ND280SFGDHit*> inputHits, int DEBUG_VERBOSE){
//***********************************************************************************************

    vector <ND280SFGDVoxel*> outputVoxels;

    struct pseudoVoxel{
        Int_t         View;
        Double_t      X;
        Double_t      Y;
        Double_t      Z;
        ND280SFGDHit* hit1;
        ND280SFGDHit* hit2;
        ND280SFGDHit* hit3;
        Int_t         found;
    };

    vector <pseudoVoxel> vXYZ;
    vXYZ.clear();

    vector <Int_t> analyzedHits(inputHits.size(),0);

    for(UInt_t ihit=0; ihit<inputHits.size(); ihit++){
        if(analyzedHits[ihit]) continue;
        else analyzedHits[ihit] = 1;
        if(DEBUG_VERBOSE>1){            
            cout << "ihit: " << ihit << endl;
            cout << inputHits[ihit]->GetX() << ", " << inputHits[ihit]->GetY() << ", " << inputHits[ihit]->GetZ() << endl;
        }
        for(UInt_t jhit=0; jhit<inputHits.size(); jhit++){
            if(ihit == jhit) continue;
            if(inputHits[ihit]->GetView() == inputHits[jhit]->GetView()) continue;

            if(inputHits[ihit]->GetView() == 0 && inputHits[jhit]->GetView() == 1 && inputHits[ihit]->GetX() == inputHits[jhit]->GetX()){
                pseudoVoxel Vox;
                Vox.X = inputHits[ihit]->GetX();
                Vox.Y = inputHits[ihit]->GetY();
                Vox.Z = inputHits[jhit]->GetZ();
                Vox.View = 0;
                Vox.hit1 = inputHits[ihit];
                Vox.hit2 = inputHits[jhit];
                vXYZ.push_back(Vox);
                if(DEBUG_VERBOSE>1){                    
                    cout << "jhit: " << jhit << endl;
                    cout << inputHits[jhit]->GetX() << ", " << inputHits[jhit]->GetY() << ", " << inputHits[jhit]->GetZ() << endl;
                }
            }
            if(inputHits[ihit]->GetView() == 0 && inputHits[jhit]->GetView() == 2  && inputHits[ihit]->GetY() == inputHits[jhit]->GetY()){
                pseudoVoxel Vox;
                Vox.X = inputHits[ihit]->GetX();
                Vox.Y = inputHits[ihit]->GetY();
                Vox.Z = inputHits[jhit]->GetZ();
                Vox.View = 1;
                Vox.hit1 = inputHits[ihit];
                Vox.hit3 = inputHits[jhit];
                vXYZ.push_back(Vox);
                // if(DEBUG_VERBOSE){
                //     cout << "jhit: " << jhit << endl;
                //     cout << inputHits[jhit]->GetX() << ", " << inputHits[jhit]->GetY() << ", " << inputHits[jhit]->GetZ() << endl;    
                // }

            }
            if(inputHits[ihit]->GetView() == 1 && inputHits[jhit]->GetView() == 2  && inputHits[ihit]->GetZ() == inputHits[jhit]->GetZ()){
                pseudoVoxel Vox;
                Vox.X = inputHits[ihit]->GetX();
                Vox.Y = inputHits[jhit]->GetY();
                Vox.Z = inputHits[jhit]->GetZ();
                Vox.View = 2;
                Vox.hit2 = inputHits[ihit];
                Vox.hit3 = inputHits[jhit];
                vXYZ.push_back(Vox);
                // if(DEBUG_VERBOSE){
                //     cout << "jhit: " << jhit << endl;
                //     cout << inputHits[jhit]->GetX() << ", " << inputHits[jhit]->GetY() << ", " << inputHits[jhit]->GetZ() << endl;
                // }
            }
        }
    }

    if(DEBUG_VERBOSE) cout << "# of 2D matches: " << vXYZ.size() << endl;

    if(DEBUG_VERBOSE){        
        Int_t CNT = 0;
        for(UInt_t ivox=0; ivox<vXYZ.size(); ivox++){
            for(UInt_t jvox=0; jvox<vXYZ.size(); jvox++){
                if(ivox == jvox) continue;
                if(vXYZ[ivox].X == vXYZ[jvox].X && vXYZ[ivox].Y == vXYZ[jvox].Y && vXYZ[ivox].Z == vXYZ[jvox].Z && vXYZ[ivox].View == vXYZ[jvox].View){
                    CNT++;
                }
            }     
        }
        cout << "# of duplicated 2D matches: " << CNT << endl;
    }

    if(DEBUG_VERBOSE>1){        
        for(UInt_t ivox=0; ivox<vXYZ.size(); ivox++){
            cout << "XYZ-View:  " << vXYZ[ivox].X << ", " << vXYZ[ivox].Y<< ", "  << vXYZ[ivox].Z << ", " << vXYZ[ivox].View << endl; 
        }
    }

    vector <Int_t> analyzedVoxels(vXYZ.size(),0);

    for(UInt_t ivox=0; ivox<vXYZ.size(); ivox++) vXYZ[ivox].found = 0;
    
    for(UInt_t ivox=0; ivox<vXYZ.size(); ivox++){
        Int_t found = 0;
        if(analyzedVoxels[ivox]) continue;
        analyzedVoxels[ivox] = 1;
        for(UInt_t jvox=0; jvox<vXYZ.size(); jvox++){
            if(analyzedVoxels[jvox]) continue;
            if(vXYZ[ivox].X == vXYZ[jvox].X && vXYZ[ivox].Y == vXYZ[jvox].Y && vXYZ[ivox].Z == vXYZ[jvox].Z && vXYZ[ivox].View != vXYZ[jvox].View){
                analyzedVoxels[jvox] = 1;
                found++;
                if(vXYZ[ivox].found){
                    vXYZ[ivox].found = found;
                    continue;
                }
                ND280SFGDVoxel* newVoxel = new ND280SFGDVoxel(vXYZ[ivox].X,vXYZ[ivox].Y,vXYZ[ivox].Z);
                vector <ND280SFGDHit*> hitsInVoxel;
                hitsInVoxel.resize(3);
                if(vXYZ[ivox].View == 0){              
                    hitsInVoxel[0] = vXYZ[ivox].hit1;
                    hitsInVoxel[1] = vXYZ[ivox].hit2;
                    hitsInVoxel[2] = vXYZ[jvox].hit3;
                }
                else if(vXYZ[ivox].View == 1){                    
                    hitsInVoxel[0] = vXYZ[ivox].hit1;
                    hitsInVoxel[1] = vXYZ[jvox].hit2;
                    hitsInVoxel[2] = vXYZ[ivox].hit3;
                }
                else if(vXYZ[ivox].View == 2){                    
                    hitsInVoxel[0] = vXYZ[jvox].hit1;
                    hitsInVoxel[1] = vXYZ[ivox].hit2;
                    hitsInVoxel[2] = vXYZ[ivox].hit3;
                }
                else{
                    cerr << "View must be 0,1 or 2!" << endl;
                    exit(1);
                }
                newVoxel->SetHits(hitsInVoxel);
                outputVoxels.push_back(newVoxel);
                vXYZ[ivox].found = found;
            }
        }
    }

    if(DEBUG_VERBOSE){        
        cout << "Voxels:" << endl;
        for(UInt_t ivox=0; ivox<outputVoxels.size(); ivox++){
            cout << ivox << "| XYZ: " << outputVoxels[ivox]->GetX() << ", " << outputVoxels[ivox]->GetY() << ", " <<  outputVoxels[ivox]->GetZ() << endl;   
        }
    }

    if(DEBUG_VERBOSE){        
        Int_t CNT = 0;
        for(UInt_t ivox=0; ivox<outputVoxels.size(); ivox++){
            for(UInt_t jvox=0; jvox<outputVoxels.size(); jvox++){
                if(ivox == jvox) continue;
                if(outputVoxels[ivox]->GetX() == outputVoxels[jvox]->GetX() && outputVoxels[ivox]->GetY() == outputVoxels[jvox]->GetY() && outputVoxels[ivox]->GetZ() == outputVoxels[jvox]->GetZ() ){
                    CNT++;
                }
            }     
        }
        cout << "# of duplicated voxels: " << CNT << endl;
    }

    return outputVoxels;
}



//***********************************************************************************************
vector <ND280SFGDVoxel*> HitsToVoxels_version1(vector <ND280SFGDHit*> inputHits, int DEBUG_VERBOSE){
//***********************************************************************************************

    vector <ND280SFGDVoxel*> outputVoxels;

    cout << "HitsToVoxels version 1 is not implemented yet!" << endl;

    return outputVoxels;
}


//***********************************************************************************************
vector <ND280SFGDVoxel*> HitsToVoxels(vector <ND280SFGDHit*> inputHits, Int_t version, int DEBUG_VERBOSE=0){
//***********************************************************************************************
    
    vector <ND280SFGDVoxel*> outputVoxels;

    if( version == 0 ){
        // the goal of the algorithm is to create voxels at X,Y,Z positions when 2 coordinates
        // matching produce a voxel in 2 different complementary 2D planes.
        outputVoxels = HitsToVoxels_version0(inputHits, DEBUG_VERBOSE);
    }
    else if ( version == 1){
        outputVoxels = HitsToVoxels_version1(inputHits, DEBUG_VERBOSE);
    }
    else{
        cerr << "Version" << version << "does not exist in HitsToVoxels function!" << endl;
        exit(1);
    }
    return outputVoxels;
}
