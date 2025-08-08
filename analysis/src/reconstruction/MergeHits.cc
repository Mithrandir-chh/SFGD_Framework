
/////// ALGORITHMS TO MERGE 2D HITS INTO 3D VOXELS //////////

//***********************************************************************************************
vector <ND280SFGDHit*> MergeHits_version0(vector <ND280SFGDHit*> inputHits, bool DEBUG){
//***********************************************************************************************

    vector <ND280SFGDHit*> outputHits;
    vector <int> analyzedHits(inputHits.size(),0);

    if(DEBUG) std::cout << "Original # of hits: " << inputHits.size() << std::endl;

    for(UInt_t ihit=0; ihit<inputHits.size(); ihit++) if(inputHits[ihit]->GetPE()<1) analyzedHits[ihit] = 1;

    for(UInt_t ihit=0; ihit<inputHits.size(); ihit++){
        if(analyzedHits[ihit]) continue;
        outputHits.push_back(inputHits[ihit]);
        analyzedHits[ihit]=1;
        for(UInt_t jhit=0; jhit<inputHits.size(); jhit++){
            if(analyzedHits[jhit]) continue;
            if(inputHits[ihit]->GetView() == inputHits[jhit]->GetView() && inputHits[ihit]->GetX() == inputHits[jhit]->GetX() && inputHits[ihit]->GetY() == inputHits[jhit]->GetY() && inputHits[ihit]->GetZ() == inputHits[jhit]->GetZ() ){
                analyzedHits[jhit] = 1;
                inputHits[ihit]->SetPE( inputHits[ihit]->GetPE() + inputHits[jhit]->GetPE() );
            }
        }
    }

    if(DEBUG) std::cout << "# of merged hits: " << outputHits.size() << std::endl;

    if(DEBUG){
        int iniPE = 0;
        int finPE= 0;
        for(auto h:inputHits)  iniPE += h->GetPE();
        for(auto h:outputHits) finPE += h->GetPE();
        std::cout << "Initial total PE: " << iniPE << std::endl;
        std::cout << "Final   total PE: " << finPE << std::endl;
    }
    return outputHits;
}


//***********************************************************************************************
vector <ND280SFGDHit*> MergeHits_version1(vector <ND280SFGDHit*> inputHits, bool DEBUG){
//***********************************************************************************************

    vector <ND280SFGDHit*> outputHits;

    cout << "MergeHits version 1 is not implemented yet!" << endl;

    return outputHits;
}


//***********************************************************************************************
vector <ND280SFGDHit*> MergeHits(vector <ND280SFGDHit*> inputHits, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    vector <ND280SFGDHit*> outputHits;

    if( version == 0 ){
        // the first version of this algorithm merges hits regardless of its time information.
        // all hits in the same fiber are merged into a single hit with #PE equal to the sum of all the merged hits.
        outputHits = MergeHits_version0(inputHits, DEBUG);
    }
    else if ( version == 1){
        outputHits = MergeHits_version1(inputHits, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in MergeHits function!" << endl;
        exit(1);
    }
    return outputHits;
}
