// SingleTrackFilter.cpp - Filter for identifying single track events from EventDisplayLite output
// Compatible with the SFGD framework build system

#define THIS_NAME SingleTrackFilter
#define OVERRIDE_OPTIONS

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TKey.h"
#include "TMath.h"
#include "TApplication.h"

// Include the framework headers FIRST (before defining any globals)
#include "../../src/tools/global_header.hh"

using namespace std;

// Detector geometry constants
const int DETECTOR_X = 24;
const int DETECTOR_Y = 8;
const int DETECTOR_Z = 48;

// Filter criteria thresholds
const double MAX_X_RANGE = 12.0;
const double MAX_Y_RANGE = 4.0;
const double MAX_Z_RANGE = 12.0;
const double MAX_CHI2_NDF = 2.0;
const double MIN_R_SQUARED = 0.8;
const int MIN_HITS_PER_VIEW = 3;
const double AUTO_PASS_X_RANGE = 3.0;  
const double AUTO_PASS_Y_RANGE = 2.0;  
const double AUTO_PASS_Z_RANGE = 7.0;  

// Structure to hold fit quality information
struct FitQuality {
    bool isGoodFit;
    double chi2_ndf;
    double rSquared;
    int nHits;
    double slope;
    double intercept;
    
    FitQuality() : isGoodFit(false), chi2_ndf(999.0), rSquared(0.0), nHits(0), slope(0.0), intercept(0.0) {}
};

// Structure to hold event data
struct EventData {
    int eventNumber;
    TH2F* hXY;
    TH2F* hXZ;
    TH2F* hZY;
    
    EventData() : eventNumber(-1), hXY(nullptr), hXZ(nullptr), hZY(nullptr) {}
    
    ~EventData() {
        if (hXY) delete hXY;
        if (hXZ) delete hXZ;
        if (hZY) delete hZY;
    }
};

// Function to extract hit positions from a 2D histogram
TGraphErrors* extractHitPositions(TH2F* hist, string& viewType) {
    vector<double> xPos, yPos;
    
    // Extract non-zero bins
    for (int ix = 1; ix <= hist->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hist->GetNbinsY(); iy++) {
            if (hist->GetBinContent(ix, iy) > 0) {
                // Get bin centers (0-indexed for detector coordinates)
                double x = hist->GetXaxis()->GetBinCenter(ix);
                double y = hist->GetYaxis()->GetBinCenter(iy);
                xPos.push_back(x);
                yPos.push_back(y);
            }
        }
    }
    
    int nHits = xPos.size();
    if (nHits < MIN_HITS_PER_VIEW) return nullptr;
    
    // Create TGraphErrors with hit positions
    TGraphErrors* graph = new TGraphErrors(nHits);
    for (int i = 0; i < nHits; i++) {
        graph->SetPoint(i, xPos[i], yPos[i]);
        graph->SetPointError(i, 0.5, 0.5); // Bin width uncertainty
    }
    
    return graph;
}

// Function to calculate R-squared value
double calculateRSquared(TGraphErrors* graph, TF1* fit) {
    int nPoints = graph->GetN();
    double sumY = 0, sumY2 = 0;
    
    // Calculate mean of y values
    for (int i = 0; i < nPoints; i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        sumY += y;
        sumY2 += y * y;
    }
    double meanY = sumY / nPoints;
    
    // Calculate SSres and SStot
    double SSres = 0, SStot = 0;
    for (int i = 0; i < nPoints; i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        double yFit = fit->Eval(x);
        SSres += pow(y - yFit, 2);
        SStot += pow(y - meanY, 2);
    }
    
    // Avoid division by zero
    if (SStot == 0) return 0;
    
    return 1.0 - (SSres / SStot);
}

// Function to evaluate linear fit quality for a view
FitQuality evaluateLinearFit(TH2F* hist, const string& viewName) {
    FitQuality result;
    
    string viewType = viewName;
    TGraphErrors* graph = extractHitPositions(hist, viewType);
    if (!graph) {
        return result; // Not enough hits
    }
    
    result.nHits = graph->GetN();
    
    // Calculate hit ranges
    double xMin = 999, xMax = -999, yMin = 999, yMax = -999;
    for (int i = 0; i < result.nHits; i++) {
        double x, y;
        graph->GetPoint(i, x, y);
        xMin = min(xMin, x); xMax = max(xMax, x);
        yMin = min(yMin, y); yMax = max(yMax, y);
    }
    
    double xRange = xMax - xMin;
    double yRange = yMax - yMin;
    
    // Determine auto-pass thresholds based on view and axis
    bool autoPass = false;
    string autoPassReason = "";
    
    if (viewName == "XY") {
        if (xRange < AUTO_PASS_X_RANGE) {
            autoPass = true;
            autoPassReason = "X-range < " + to_string(AUTO_PASS_X_RANGE);
        } else if (yRange < AUTO_PASS_Y_RANGE) {
            autoPass = true;
            autoPassReason = "Y-range < " + to_string(AUTO_PASS_Y_RANGE);
        }
    } else if (viewName == "XZ") {
        if (xRange < AUTO_PASS_X_RANGE) {
            autoPass = true;
            autoPassReason = "X-range < " + to_string(AUTO_PASS_X_RANGE);
        } else if (yRange < AUTO_PASS_Z_RANGE) { // y-axis represents Z in XZ view
            autoPass = true;
            autoPassReason = "Z-range < " + to_string(AUTO_PASS_Z_RANGE);
        }
    } else if (viewName == "ZY") {
        if (xRange < AUTO_PASS_Z_RANGE) { // x-axis represents Z in ZY view
            autoPass = true;
            autoPassReason = "Z-range < " + to_string(AUTO_PASS_Z_RANGE);
        } else if (yRange < AUTO_PASS_Y_RANGE) {
            autoPass = true;
            autoPassReason = "Y-range < " + to_string(AUTO_PASS_Y_RANGE);
        }
    }
    
    if (autoPass) {
        // Automatically pass fit quality tests
        result.isGoodFit = true;
        result.chi2_ndf = 0.1;  // Excellent chi2/ndf
        result.rSquared = 0.99; // Excellent R-squared
        result.slope = (yRange > xRange) ? 999.0 : 0.0;
        result.intercept = (yRange > xRange) ? xMin : yMin;
        
    } else {
        // Standard linear fit for tracks with significant spread in both dimensions
        TF1* fitFunc = new TF1("linearFit", "pol1", xMin - 1, xMax + 1);
        graph->Fit(fitFunc, "Q");
        
        result.slope = fitFunc->GetParameter(1);
        result.intercept = fitFunc->GetParameter(0);
        
        double chi2 = fitFunc->GetChisquare();
        int ndf = fitFunc->GetNDF();
        result.chi2_ndf = (ndf > 0) ? chi2 / ndf : 999.0;
        result.rSquared = calculateRSquared(graph, fitFunc);
        
        result.isGoodFit = (result.chi2_ndf < MAX_CHI2_NDF && 
                           result.rSquared > MIN_R_SQUARED && 
                           result.nHits >= MIN_HITS_PER_VIEW);

        delete fitFunc;
    }
    
    delete graph;
    return result;
}

// Function to check geometric range
bool passesRangeCheck(TH2F* hXY, TH2F* hXZ, TH2F* hZY, 
                     double& xRange, double& yRange, double& zRange) {
    double xMin = 999, xMax = -999;
    double yMin = 999, yMax = -999;
    double zMin = 999, zMax = -999;
    
    // Scan XY view
    for (int ix = 1; ix <= hXY->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hXY->GetNbinsY(); iy++) {
            if (hXY->GetBinContent(ix, iy) > 0) {
                double x = hXY->GetXaxis()->GetBinCenter(ix);
                double y = hXY->GetYaxis()->GetBinCenter(iy);
                xMin = min(xMin, x);
                xMax = max(xMax, x);
                yMin = min(yMin, y);
                yMax = max(yMax, y);
            }
        }
    }
    
    // Scan XZ view
    for (int ix = 1; ix <= hXZ->GetNbinsX(); ix++) {
        for (int iz = 1; iz <= hXZ->GetNbinsY(); iz++) {
            if (hXZ->GetBinContent(ix, iz) > 0) {
                double x = hXZ->GetXaxis()->GetBinCenter(ix);
                double z = hXZ->GetYaxis()->GetBinCenter(iz);
                xMin = min(xMin, x);
                xMax = max(xMax, x);
                zMin = min(zMin, z);
                zMax = max(zMax, z);
            }
        }
    }
    
    // Scan ZY view
    for (int iz = 1; iz <= hZY->GetNbinsX(); iz++) {
        for (int iy = 1; iy <= hZY->GetNbinsY(); iy++) {
            if (hZY->GetBinContent(iz, iy) > 0) {
                double z = hZY->GetXaxis()->GetBinCenter(iz);
                double y = hZY->GetYaxis()->GetBinCenter(iy);
                zMin = min(zMin, z);
                zMax = max(zMax, z);
                yMin = min(yMin, y);
                yMax = max(yMax, y);
            }
        }
    }
    
    xRange = xMax - xMin;
    yRange = yMax - yMin;
    zRange = zMax - zMin;
    bool xPass = xRange < MAX_X_RANGE;
    bool yPass = yRange < MAX_Y_RANGE;
    bool zPass = zRange < MAX_Z_RANGE;

    // return ((xPass && yPass) || (yPass && zPass) || (zPass && xPass));
    return (xPass && yPass && zPass);
}

// Helper function to check if a hit is at a detector corner
bool isCornerHit(double x, double y, double z, const string& viewType) {
    if (viewType == "XY") {
        // XY view corners: (0,0), (0,7), (23,0), (23,7)
        return ((x <= 0.5 || x >= 23.5) && (y <= 0.5 || y >= 7.5));
    } else if (viewType == "XZ") {
        // XZ view corners: (0,0), (0,47), (23,0), (23,47)
        return ((x <= 0.5 || x >= 23.5) && (z <= 0.5 || z >= 47.5));
    } else if (viewType == "ZY") {
        // ZY view corners: (0,0), (0,7), (47,0), (47,7)
        return ((z <= 0.5 || z >= 47.5) && (y <= 0.5 || y >= 7.5));
    }
    return false;
}

// Modified boundary check that excludes corner hits
bool hasBoundaryHits(TH2F* hXY, TH2F* hXZ, TH2F* hZY) {
    bool xyBoundary = false;
    bool xzBoundary = false;
    bool zyBoundary = false;
    
    // Check XY view
    vector<pair<double, double>> xyHits;
    bool hasCornerHit = false;
    
    for (int ix = 1; ix <= hXY->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hXY->GetNbinsY(); iy++) {
            if (hXY->GetBinContent(ix, iy) > 0) {
                double x = hXY->GetXaxis()->GetBinCenter(ix);
                double y = hXY->GetYaxis()->GetBinCenter(iy);
                
                // Check if this hit is at a corner
                if (isCornerHit(x, y, 0, "XY")) {
                    hasCornerHit = true;
                    break;  // Immediately reject if corner hit found
                }
                
                xyHits.push_back({x, y});
            }
        }
        if (hasCornerHit) break;
    }
    
    // If corner hit found, return false immediately
    if (hasCornerHit) return false;
    
    // Continue with normal boundary check for XY view
    vector<int> xBoundaryHits, yBoundaryHits;
    for (size_t i = 0; i < xyHits.size(); i++) {
        double x = xyHits[i].first;
        double y = xyHits[i].second;
        
        // Check X boundaries (but not corners since we already excluded them)
        if (x <= 0.5 || x >= 23.5) {
            xBoundaryHits.push_back(i);
        }
        
        // Check Y boundaries
        if (y <= 0.5 || y >= 7.5) {
            yBoundaryHits.push_back(i);
        }
    }
    
    // Check if we have both types from different hits
    if (!xBoundaryHits.empty() && !yBoundaryHits.empty()) {
        for (int xIdx : xBoundaryHits) {
            for (int yIdx : yBoundaryHits) {
                if (xIdx != yIdx) {
                    xyBoundary = true;
                    break;
                }
            }
            if (xyBoundary) break;
        }
    }
    
    // Check XZ view
    vector<pair<double, double>> xzHits;
    hasCornerHit = false;
    
    for (int ix = 1; ix <= hXZ->GetNbinsX(); ix++) {
        for (int iz = 1; iz <= hXZ->GetNbinsY(); iz++) {
            if (hXZ->GetBinContent(ix, iz) > 0) {
                double x = hXZ->GetXaxis()->GetBinCenter(ix);
                double z = hXZ->GetYaxis()->GetBinCenter(iz);
                
                if (isCornerHit(x, 0, z, "XZ")) {
                    hasCornerHit = true;
                    break;
                }
                
                xzHits.push_back({x, z});
            }
        }
        if (hasCornerHit) break;
    }
    
    if (hasCornerHit) return false;
    
    vector<int> xBoundaryHitsXZ, zBoundaryHits;
    for (size_t i = 0; i < xzHits.size(); i++) {
        double x = xzHits[i].first;
        double z = xzHits[i].second;
        
        if (x <= 0.5 || x >= 23.5) {
            xBoundaryHitsXZ.push_back(i);
        }
        if (z <= 0.5 || z >= 47.5) {
            zBoundaryHits.push_back(i);
        }
    }
    
    if (!xBoundaryHitsXZ.empty() && !zBoundaryHits.empty()) {
        for (int xIdx : xBoundaryHitsXZ) {
            for (int zIdx : zBoundaryHits) {
                if (xIdx != zIdx) {
                    xzBoundary = true;
                    break;
                }
            }
            if (xzBoundary) break;
        }
    }
    
    // Check ZY view
    vector<pair<double, double>> zyHits;
    hasCornerHit = false;
    
    for (int iz = 1; iz <= hZY->GetNbinsX(); iz++) {
        for (int iy = 1; iy <= hZY->GetNbinsY(); iy++) {
            if (hZY->GetBinContent(iz, iy) > 0) {
                double z = hZY->GetXaxis()->GetBinCenter(iz);
                double y = hZY->GetYaxis()->GetBinCenter(iy);
                
                if (isCornerHit(0, y, z, "ZY")) {
                    hasCornerHit = true;
                    break;
                }
                
                zyHits.push_back({z, y});
            }
        }
        if (hasCornerHit) break;
    }
    
    if (hasCornerHit) return false;
    
    vector<int> zBoundaryHitsZY, yBoundaryHitsZY;
    for (size_t i = 0; i < zyHits.size(); i++) {
        double z = zyHits[i].first;
        double y = zyHits[i].second;
        
        if (z <= 0.5 || z >= 47.5) {
            zBoundaryHitsZY.push_back(i);
        }
        if (y <= 0.5 || y >= 7.5) {
            yBoundaryHitsZY.push_back(i);
        }
    }
    
    if (!zBoundaryHitsZY.empty() && !yBoundaryHitsZY.empty()) {
        for (int zIdx : zBoundaryHitsZY) {
            for (int yIdx : yBoundaryHitsZY) {
                if (zIdx != yIdx) {
                    zyBoundary = true;
                    break;
                }
            }
            if (zyBoundary) break;
        }
    }
    
    return (xyBoundary || xzBoundary || zyBoundary);
}

// Function to determine if event contains a single track
bool isSingleTrackEvent(const FitQuality& xy, const FitQuality& xz, const FitQuality& zy) {
    int goodFits = 0;
    if (xy.isGoodFit) goodFits++;
    if (xz.isGoodFit) goodFits++;
    if (zy.isGoodFit) goodFits++;
    
    return (goodFits >= 2);
}

// Main function following framework convention
void SingleTrackFilter() {
    
    // Default parameters
    TString inputFile = "";
    TString outputFile = "single_track_filtered.root";
    TString acceptDir = "./AcceptedEvents/";
    TString rejectDir = "./RejectedEvents/";
    bool saveImages = true;
    bool verbose = false;
    
    // Get command line arguments from the global application
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Parse command line arguments
    for (int iarg = 0; iarg < argc; iarg++) {
        string arg = argv[iarg];
        
        if (arg == "-h" || arg == "--help") {
            cout << "Single Track Filter for SFGD Event Display\n";
            cout << "Usage: " << argv[0] << " -i input.root [options]\n";
            cout << "Options:\n";
            cout << "  -i FILE              Input ROOT file from EventDisplayLite\n";
            cout << "  -o FILE              Output ROOT file (default: single_track_filtered.root)\n";
            cout << "  --accept-dir DIR     Directory for accepted event images (default: ./AcceptedEvents/)\n";
            cout << "  --reject-dir DIR     Directory for rejected event images (default: ./RejectedEvents/)\n";
            cout << "  --no-images          Don't save PNG images\n";
            cout << "  -v, --verbose        Verbose output\n";
            cout << "  -h, --help           Show this help\n";
            cout << "\nFilter Criteria:\n";
            cout << "  - Compact events (X<" << MAX_X_RANGE << ", Y<" << MAX_Y_RANGE 
                 << ", Z<" << MAX_Z_RANGE << ") must have boundary hits\n";
            cout << "  - All events must have good linear fits in at least 2 views\n";
            cout << "  - Linear fit criteria: chi2/ndf<" << MAX_CHI2_NDF 
                 << ", R2>" << MIN_R_SQUARED << ", hits>=" << MIN_HITS_PER_VIEW << "\n";
            return;
        }
        else if ((arg == "-i" || arg == "--input") && iarg + 1 < argc) {
            iarg++;
            inputFile = argv[iarg];
        }
        else if ((arg == "-o" || arg == "--output") && iarg + 1 < argc) {
            iarg++;
            outputFile = argv[iarg];
        }
        else if (arg == "--accept-dir" && iarg + 1 < argc) {
            iarg++;
            acceptDir = argv[iarg];
            if (!acceptDir.EndsWith("/")) acceptDir += "/";
        }
        else if (arg == "--reject-dir" && iarg + 1 < argc) {
            iarg++;
            rejectDir = argv[iarg];
            if (!rejectDir.EndsWith("/")) rejectDir += "/";
        }
        else if (arg == "--no-images") {
            saveImages = false;
        }
        else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
    }
    
    if (inputFile == "") {
        cerr << "Error: Input file required. Use -h for help.\n";
        return;
    }
    
    TString inputDir = gSystem->DirName(inputFile);
    TString inputBaseName = gSystem->BaseName(inputFile);
    inputBaseName.ReplaceAll(".root", "");
    
    // If no custom output specified, create in same directory as input
    if (outputFile == "single_track_filtered.root") {
        outputFile = inputDir + "/" + inputBaseName + "_singletrack.root";
    }
    
    // Set default directories relative to input file location if not specified by user
    if (acceptDir == "./AcceptedEvents/") {
        acceptDir = inputDir + "/" + "AcceptedEvents/";
    }
    if (rejectDir == "./RejectedEvents/") {
        rejectDir = inputDir + "/" + "RejectedEvents/";
    }
    
    // Ensure directories end with slash
    if (!acceptDir.EndsWith("/")) acceptDir += "/";
    if (!rejectDir.EndsWith("/")) rejectDir += "/";
    
    // Set style
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    
    // Create output directories
    if (saveImages) {
        gSystem->mkdir(acceptDir.Data(), kTRUE);
        gSystem->mkdir(rejectDir.Data(), kTRUE);
    }
    
    // Open input file
    TFile* fInput = TFile::Open(inputFile.Data(), "READ");
    if (!fInput || fInput->IsZombie()) {
        cerr << "Error: Cannot open input file: " << inputFile << endl;
        return;
    }
    
    // Open output file
    TFile* fOutput = new TFile(outputFile.Data(), "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        cerr << "Error: Cannot create output file: " << outputFile << endl;
        fInput->Close();
        return;
    }
    
    // Create diagnostic histograms
    TH1F* hRangeX = new TH1F("hRangeX", "X Range Distribution;X Range [units];Events", 50, 0, 25);
    TH1F* hRangeY = new TH1F("hRangeY", "Y Range Distribution;Y Range [units];Events", 50, 0, 10);
    TH1F* hRangeZ = new TH1F("hRangeZ", "Z Range Distribution;Z Range [units];Events", 50, 0, 50);
    
    TH1F* hChi2XY = new TH1F("hChi2XY", "XY View #chi^{2}/ndf;#chi^{2}/ndf;Events", 100, 0, 10);
    TH1F* hChi2XZ = new TH1F("hChi2XZ", "XZ View #chi^{2}/ndf;#chi^{2}/ndf;Events", 100, 0, 10);
    TH1F* hChi2ZY = new TH1F("hChi2ZY", "ZY View #chi^{2}/ndf;#chi^{2}/ndf;Events", 100, 0, 10);
    
    TH1F* hR2XY = new TH1F("hR2XY", "XY View R^{2};R^{2};Events", 100, 0, 1);
    TH1F* hR2XZ = new TH1F("hR2XZ", "XZ View R^{2};R^{2};Events", 100, 0, 1);
    TH1F* hR2ZY = new TH1F("hR2ZY", "ZY View R^{2};R^{2};Events", 100, 0, 1);
    
    TH1F* hSlopeXY = new TH1F("hSlopeXY", "XY Track Slopes;Slope;Events", 100, -5, 5);
    TH1F* hSlopeXZ = new TH1F("hSlopeXZ", "XZ Track Slopes;Slope;Events", 100, -5, 5);
    TH1F* hSlopeZY = new TH1F("hSlopeZY", "ZY Track Slopes;Slope;Events", 100, -5, 5);
    
    TH1F* hFilterSteps = new TH1F("hFilterSteps", "Events Passing Each Filter Step", 5, 0, 5);
    hFilterSteps->GetXaxis()->SetBinLabel(1, "Total Events");
    hFilterSteps->GetXaxis()->SetBinLabel(2, "Events Analyzed");
    hFilterSteps->GetXaxis()->SetBinLabel(3, "Pass Boundary/Range");
    hFilterSteps->GetXaxis()->SetBinLabel(4, "Pass Linear Fit");
    hFilterSteps->GetXaxis()->SetBinLabel(5, "Single Track Events");
    
    // Scan input file for event histograms
    map<int, EventData*> events;
    TIter next(fInput->GetListOfKeys());
    TKey* key;
    
    cout << "Scanning input file for event histograms..." << endl;
    
    while ((key = (TKey*)next())) {
        string keyName = key->GetName();
        if (keyName.find("Event") != 0) continue;
        
        // Extract event number and view type
        int eventNum;
        char viewType[10];
        if (sscanf(keyName.c_str(), "Event%d_%[^_]_Filtered", &eventNum, viewType) == 2) {
            // Create event entry if doesn't exist
            if (events.find(eventNum) == events.end()) {
                events[eventNum] = new EventData();
                events[eventNum]->eventNumber = eventNum;
            }
            
            // Load histogram
            TH2F* hist = (TH2F*)key->ReadObj();
            string view = viewType;
            
            if (view == "XY") {
                events[eventNum]->hXY = (TH2F*)hist->Clone();
            } else if (view == "XZ") {
                events[eventNum]->hXZ = (TH2F*)hist->Clone();
            } else if (view == "ZY") {
                events[eventNum]->hZY = (TH2F*)hist->Clone();
            }
        }
    }
    
    cout << "Found " << events.size() << " events to process" << endl;
    
    // Canvas for saving images
    TCanvas* canvas = nullptr;
    if (saveImages) {
        canvas = new TCanvas("SingleTrackDisplay", "Single Track Event Display", 1200, 900);
        canvas->Divide(2, 2);
    }
    
    // Process each event
    int nTotal = 0;
    int nPassRange = 0;
    int nPassBoundary = 0;
    int nPassLinearFit = 0;
    int nSingleTrack = 0;
    
    for (auto& eventPair : events) {
        EventData* event = eventPair.second;
        nTotal++;
        
        if (nTotal % 10 == 0) {
            cout << "Processing event " << nTotal << "/" << events.size() << "\r" << flush;
        }
        
        // Check that all views are present
        if (!event->hXY || !event->hXZ || !event->hZY) {
            if (verbose) {
                cout << "Event " << event->eventNumber << ": Missing views, skipping" << endl;
            }
            continue;
        }
        
        // Step 1: Calculate geometric ranges
        double xRange, yRange, zRange;
        bool isCompactEvent = passesRangeCheck(event->hXY, event->hXZ, event->hZY, xRange, yRange, zRange);
        
        hRangeX->Fill(xRange);
        hRangeY->Fill(yRange);
        hRangeZ->Fill(zRange);
        
        // Step 2: Apply boundary check ONLY for compact events
        bool needsBoundaryCheck = isCompactEvent;  // Events with small ranges need boundary validation
        bool passBoundary = true;  // Default to true for events that don't need boundary check
        
        if (needsBoundaryCheck) {
            passBoundary = hasBoundaryHits(event->hXY, event->hXZ, event->hZY);
            
            if (!passBoundary) {
                if (verbose) {
                    cout << "Event " << event->eventNumber 
                         << ": Compact event (X:" << xRange << ", Y:" << yRange 
                         << ", Z:" << zRange << ") failed boundary check" << endl;
                }
                // Save to rejected folder
                if (saveImages && canvas) {
                    canvas->cd(1); event->hXY->Draw("COLZ"); gPad->SetRightMargin(0.15);
                    canvas->cd(2); event->hXZ->Draw("COLZ"); gPad->SetRightMargin(0.15);
                    canvas->cd(3); event->hZY->Draw("COLZ"); gPad->SetRightMargin(0.15);
                    canvas->cd(4);
                    TPaveText* info = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
                    info->AddText(Form("Event: %d", event->eventNumber));
                    info->AddText("REJECTED: Compact Event Failed Boundary Check");
                    info->AddText(Form("Ranges: X=%.1f, Y=%.1f, Z=%.1f", xRange, yRange, zRange));
                    info->AddText("Compact events must have boundary hits");
                    info->SetTextColor(kRed);
                    info->Draw();
                    canvas->SaveAs((rejectDir + Form("event_%06d_boundary_fail.png", event->eventNumber)).Data());
                    delete info;
                }
                continue;
            }
            nPassBoundary++;
        } else {
            // Large range events automatically pass boundary check
            nPassBoundary++;
            if (verbose) {
                cout << "Event " << event->eventNumber 
                     << ": Large range event (X:" << xRange << ", Y:" << yRange 
                     << ", Z:" << zRange << ") - boundary check not required" << endl;
            }
        }
        
        nPassRange++;  // Count all events that were processed (not filtered by range)
        
        // Step 3: Linear fit quality assessment
        FitQuality fitXY = evaluateLinearFit(event->hXY, "XY");
        FitQuality fitXZ = evaluateLinearFit(event->hXZ, "XZ");
        FitQuality fitZY = evaluateLinearFit(event->hZY, "ZY");
        
        // Fill fit quality histograms
        if (fitXY.nHits >= MIN_HITS_PER_VIEW) {
            hChi2XY->Fill(fitXY.chi2_ndf);
            hR2XY->Fill(fitXY.rSquared);
            hSlopeXY->Fill(fitXY.slope);
        }
        if (fitXZ.nHits >= MIN_HITS_PER_VIEW) {
            hChi2XZ->Fill(fitXZ.chi2_ndf);
            hR2XZ->Fill(fitXZ.rSquared);
            hSlopeXZ->Fill(fitXZ.slope);
        }
        if (fitZY.nHits >= MIN_HITS_PER_VIEW) {
            hChi2ZY->Fill(fitZY.chi2_ndf);
            hR2ZY->Fill(fitZY.rSquared);
            hSlopeZY->Fill(fitZY.slope);
        }
        
        // Step 4: Multi-view consistency check
        bool isSingleTrack = isSingleTrackEvent(fitXY, fitXZ, fitZY);
        
        if (!isSingleTrack) {
            if (verbose) {
                cout << "Event " << event->eventNumber << ": Failed single track criteria" << endl;
            }
            // Save to rejected folder
            if (saveImages && canvas) {
                canvas->cd(1); event->hXY->Draw("COLZ"); gPad->SetRightMargin(0.15);
                canvas->cd(2); event->hXZ->Draw("COLZ"); gPad->SetRightMargin(0.15);
                canvas->cd(3); event->hZY->Draw("COLZ"); gPad->SetRightMargin(0.15);
                canvas->cd(4);
                TPaveText* info = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
                info->AddText(Form("Event: %d", event->eventNumber));
                info->AddText("REJECTED: Failed Linear Fit");
                info->AddText(Form("XY: #chi^{2}/ndf=%.2f, R^{2}=%.2f %s", 
                    fitXY.chi2_ndf, fitXY.rSquared, fitXY.isGoodFit ? "(GOOD)" : "(BAD)"));
                info->AddText(Form("XZ: #chi^{2}/ndf=%.2f, R^{2}=%.2f %s", 
                    fitXZ.chi2_ndf, fitXZ.rSquared, fitXZ.isGoodFit ? "(GOOD)" : "(BAD)"));
                info->AddText(Form("ZY: #chi^{2}/ndf=%.2f, R^{2}=%.2f %s", 
                    fitZY.chi2_ndf, fitZY.rSquared, fitZY.isGoodFit ? "(GOOD)" : "(BAD)"));
                info->SetTextColor(kRed);
                info->Draw();
                canvas->SaveAs((rejectDir + Form("event_%06d_fit_fail.png", event->eventNumber)).Data());
                delete info;
            }
            continue;
        }
        
        nPassLinearFit++;
        nSingleTrack++;
        
        // Event passed all filters - save to output and accepted folder
        fOutput->cd();
        event->hXY->Write(Form("Event%06d_XY_SingleTrack", event->eventNumber));
        event->hXZ->Write(Form("Event%06d_XZ_SingleTrack", event->eventNumber));
        event->hZY->Write(Form("Event%06d_ZY_SingleTrack", event->eventNumber));
        
        if (saveImages && canvas) {
            canvas->cd(1); event->hXY->Draw("COLZ"); gPad->SetRightMargin(0.15);
            canvas->cd(2); event->hXZ->Draw("COLZ"); gPad->SetRightMargin(0.15);
            canvas->cd(3); event->hZY->Draw("COLZ"); gPad->SetRightMargin(0.15);
            canvas->cd(4);
            TPaveText* info = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
            info->AddText(Form("Event: %d", event->eventNumber));
            info->AddText("ACCEPTED: Single Track Event");
            info->AddText(Form("Ranges: X=%.1f, Y=%.1f, Z=%.1f", xRange, yRange, zRange));
            info->AddText(Form("XY: #chi^{2}/ndf=%.2f, R^{2}=%.2f", fitXY.chi2_ndf, fitXY.rSquared));
            info->AddText(Form("XZ: #chi^{2}/ndf=%.2f, R^{2}=%.2f", fitXZ.chi2_ndf, fitXZ.rSquared));
            info->AddText(Form("ZY: #chi^{2}/ndf=%.2f, R^{2}=%.2f", fitZY.chi2_ndf, fitZY.rSquared));
            info->AddText(needsBoundaryCheck ? "Compact event - boundary check passed" : "Large range event");
            info->SetTextColor(kGreen+2);
            info->Draw();
            canvas->SaveAs((acceptDir + Form("event_%06d_accepted.png", event->eventNumber)).Data());
            delete info;
        }
        
        if (verbose) {
            cout << "Event " << event->eventNumber << ": ACCEPTED as single track" << endl;
        }
    }
    
    cout << "\nProcessing complete!" << endl;
    
    // Fill filter statistics
    hFilterSteps->SetBinContent(1, nTotal);
    hFilterSteps->SetBinContent(2, nPassRange);
    hFilterSteps->SetBinContent(3, nPassBoundary);
    hFilterSteps->SetBinContent(4, nPassLinearFit);
    hFilterSteps->SetBinContent(5, nSingleTrack);
    
    // Save diagnostic histograms
    fOutput->cd();
    hRangeX->Write();
    hRangeY->Write();
    hRangeZ->Write();
    hChi2XY->Write();
    hChi2XZ->Write();
    hChi2ZY->Write();
    hR2XY->Write();
    hR2XZ->Write();
    hR2ZY->Write();
    hSlopeXY->Write();
    hSlopeXZ->Write();
    hSlopeZY->Write();
    hFilterSteps->Write();
    
    // Create and save summary canvas
    if (saveImages) {
        TCanvas* summaryCanvas = new TCanvas("FilterSummary", "Filter Summary", 1400, 1000);
        summaryCanvas->Divide(3, 3);
        
        summaryCanvas->cd(1); hRangeX->Draw();
        summaryCanvas->cd(2); hRangeY->Draw();
        summaryCanvas->cd(3); hRangeZ->Draw();
        
        summaryCanvas->cd(4); hChi2XY->Draw();
        summaryCanvas->cd(5); hChi2XZ->Draw();
        summaryCanvas->cd(6); hChi2ZY->Draw();
        
        summaryCanvas->cd(7); hR2XY->Draw();
        summaryCanvas->cd(8); hR2XZ->Draw();
        summaryCanvas->cd(9); hFilterSteps->Draw();
        
        summaryCanvas->SaveAs("filter_summary.png");
        delete summaryCanvas;
    }
    
    // Print summary statistics
    cout << "\n=== FILTERING SUMMARY ===" << endl;
    cout << "Total events processed: " << nTotal << endl;
    cout << "Events analyzed (not pre-filtered by range): " << nPassRange << " (" 
         << (100.0 * nPassRange / nTotal) << "%)" << endl;
    cout << "Events passing boundary check (or not requiring it): " << nPassBoundary << " (" 
         << (100.0 * nPassBoundary / nTotal) << "%)" << endl;
    cout << "Events passing linear fit: " << nPassLinearFit << " (" 
         << (100.0 * nPassLinearFit / nTotal) << "%)" << endl;
    cout << "Single track events: " << nSingleTrack << " (" 
         << (100.0 * nSingleTrack / nTotal) << "%)" << endl;
    cout << "\nNote: Boundary check only applied to compact events (X<" << MAX_X_RANGE 
         << ", Y<" << MAX_Y_RANGE << ", Z<" << MAX_Z_RANGE << ")" << endl;
    cout << "\nOutput saved to: " << outputFile << endl;
    if (saveImages) {
        cout << "Accepted events saved to: " << acceptDir << endl;
        cout << "Rejected events saved to: " << rejectDir << endl;
    }
    
    // Cleanup
    for (auto& eventPair : events) {
        delete eventPair.second;
    }
    
    delete hRangeX;
    delete hRangeY;
    delete hRangeZ;
    delete hChi2XY;
    delete hChi2XZ;
    delete hChi2ZY;
    delete hR2XY;
    delete hR2XZ;
    delete hR2ZY;
    delete hSlopeXY;
    delete hSlopeXZ;
    delete hSlopeZY;
    delete hFilterSteps;
    
    if (canvas) delete canvas;
    
    fOutput->Close();
    delete fOutput;
    
    fInput->Close();
    delete fInput;
    
    exit(0);
}