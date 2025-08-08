
// /////// ALGORITHMS TO FIT 3D LINES TO A GIVEN SET OF VOXELS //////////

//***********************************************************************************************
std::vector<double> TrackFitter_version0(vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG){
//***********************************************************************************************
    
    std::vector<double> fitParameters;

    // function Object to be minimized
    struct SumDistance2 {
      const vector <ND280SFGDVoxel*> fVoxels;

      SumDistance2(const vector <ND280SFGDVoxel*> g) : fVoxels(g) {}

      // calculate distance line-point
      double distance2(double x,double y,double z, const double *p) {
         // distance line point is D= | (xp-x0) cross  ux |
         // where ux is direction of line and x0 is a point in the line (like t = 0)
         XYZVector xp(x,y,z);
         XYZVector x0(p[0], p[1], p[2]);
         XYZVector x1(p[0] + p[3], p[1] + p[4], p[2]+ p[5]);
         XYZVector u = (x1-x0).Unit();
         double d2 = ((xp-x0).Cross(u)).Mag2();
         return d2;
      }

      // implementation of the function to be minimized
      double operator() (const double *par) {
         assert(fVoxels.size() != 0);
         double sum = 0;
         for (auto v:fVoxels) {
            double d = distance2(v->GetX(),v->GetY(),v->GetZ(),par);
            // possible improvement, weight by RecoQ: d*=v->GetRecoQ() ---> maybe?
            sum += d;
         }
         return sum;
      }

    };

    gStyle->SetOptStat(0);
    gStyle->SetOptFit();

    gStyle->SetCanvasColor(0);
    gStyle->SetMarkerStyle(21);
    gStyle->SetMarkerSize(1.05);

    ROOT::Fit::Fitter  fitter;

    // make the functor objet
    SumDistance2 sdist(inputVoxels);
    ROOT::Math::Functor fcn(sdist,6);
    // set the function and the initial parameter values
    double pStart[6] = {1,1,1,1,1,1};
    fitter.SetFCN(fcn,pStart);
    // set step sizes different than default ones (0.3 times parameter values)
    for (int i = 0; i < 6; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

    std::vector <double> p_result;

    bool ok = fitter.FitFCN();
    if (!ok) {
        //Error("line3Dfit","Line3D Fit failed");
        return fitParameters;
    }

    const ROOT::Fit::FitResult & result = fitter.Result();

    if(DEBUG){
        std::cout << "Total final distance square: " << result.MinFcnValue() << std::endl;
        std::cout << "Ave Dist: " << result.MinFcnValue() / (inputVoxels.size()) << std::endl;
        result.Print(std::cout);
    }

    // get fit parameters
    const double * parFit = result.GetParams();
    const double * parFitErr = result.GetErrors();

    double p0 = parFit[0];
    double p1 = parFit[1];
    double p2 = parFit[2];
    double p3 = parFit[3];
    double p4 = parFit[4];
    double p5 = parFit[5];


    fitParameters.push_back(result.MinFcnValue() / (inputVoxels.size()));  // ave dist
    fitParameters.push_back(p0); // point x
    fitParameters.push_back(p1); // point y
    fitParameters.push_back(p2); // point z
    fitParameters.push_back(p3); // vector x
    fitParameters.push_back(p4); // vector y
    fitParameters.push_back(p5); // vector z
    fitParameters.push_back(parFitErr[0]);
    fitParameters.push_back(parFitErr[1]);
    fitParameters.push_back(parFitErr[2]);
    fitParameters.push_back(parFitErr[3]);
    fitParameters.push_back(parFitErr[4]);
    fitParameters.push_back(parFitErr[5]);

    // draw the fitted line
    int n = 1000;
    double t0 = 0;
    double unitVecNorm = pow(pow(p3,2)+pow(p4,2)+pow(p5,2),0.5);
    double dt = 50/unitVecNorm;
    TPolyLine3D *l = new TPolyLine3D(n);

    int x0 = 0;
    int y0 = 0;
    int z0 = 0;

    for(auto v:inputVoxels){
        x0+=v->GetX();
        y0+=v->GetY();
        z0+=v->GetZ();
    }

    x0/=inputVoxels.size();
    y0/=inputVoxels.size();
    z0/=inputVoxels.size();

    // fitParameters.push_back(x0);
    // fitParameters.push_back(y0);
    // fitParameters.push_back(z0);
    // fitParameters.push_back(unitVecNorm);

    if(DEBUG){
        for (int i = 0; i <n;++i) {
            double t = t0+ dt*i/n-dt/(2);
            double x,y,z;
            line(t,parFit,x,y,z,x0,y0,z0);
            l->SetPoint(i,x,y,z);
        }
        l->SetLineColor(kRed);
        std::cout << "drawing " << inputVoxels.size() << "voxels." << std::endl;
        TNtuple *event3D = new TNtuple("event3D", "event3D", "x:y:z:c");
        for (auto v:inputVoxels){
          event3D->Fill(v->GetZ(),v->GetY(),v->GetX(),v->GetTruePE());
        }
        TCanvas *canv = new TCanvas("canv", "canv", 800, 600, 800, 600);
        canv->cd(1);
        event3D->Draw("x:y:z:c","","box2");
        TH3F *htemp = (TH3F*)gPad->GetPrimitive("htemp");
        htemp->GetXaxis()->SetLimits(0,SFGD_X);
        htemp->GetYaxis()->SetLimits(0,SFGD_Y);
        htemp->GetZaxis()->SetLimits(0,SFGD_Z);
        htemp->SetTitle("");
        l->Draw("same");
        canv->Update();
        canv->WaitPrimitive();
        delete htemp;
        delete canv;
        delete event3D;
        cout << "Fitt completed: solution is: " << p0 << "," << p1 << "," << p2 << "," << p3 << "," << p4 << "," << p5 << endl;
    }

    return fitParameters;

}

//***********************************************************************************************
std::vector<double> TrackFitter_version1(vector <ND280SFGDVoxel*> inputVoxels, bool DEBUG){
//***********************************************************************************************
    std::vector<double> fitParameters;
    cout << "TrackFitter version 1 is not implemented yet!" << endl;
    return fitParameters;
}


//***********************************************************************************************
std::vector<double> TrackFitter(vector <ND280SFGDVoxel*> inputVoxels, Int_t version, bool DEBUG=false){
//***********************************************************************************************
    
    std::vector<double> fitParameters;

    if( version == 0 ){
        // Find fit estimation for a track.
        fitParameters = TrackFitter_version0(inputVoxels, DEBUG);
    }
    else if ( version == 1){
        fitParameters = TrackFitter_version1(inputVoxels, DEBUG);
    }
    else{
        cerr << "Version" << version << "does not exist in TrackFitter function!" << endl;
        exit(1);
    }
    return fitParameters;
}
