// Calculates distances from high-confidence voxels to fitted trajectory line
#define THIS_NAME DistanceAnalysis_3D
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TKey.h"
#include "TPolyLine3D.h"
#include "TF1.h"

// Include the framework headers
#include "../../src/tools/global_header.hh"

using namespace std;

// Function to create directory if it doesn't exist
TString createDistanceOutputDirectory(const TString& inputFile) {
    // Extract directory path from input file
    TString inputDir = gSystem->DirName(inputFile);
    
    // Go one level up from input directory
    TString parentDir = gSystem->DirName(inputDir);
    
    // Create the FitDistance_PNG folder path
    TString outputDir = parentDir + "/FitDistance_PNG";
    
    // Check if directory exists, create if not
    if (gSystem->AccessPathName(outputDir)) {
        std::cout << "Creating distance output directory: " << outputDir << std::endl;
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            std::cerr << "Using current directory instead." << std::endl;
            return "./FitDistance_PNG";
        }
    } else {
        std::cout << "Using existing distance output directory: " << outputDir << std::endl;
    }
    
    return outputDir;
}

// Function to calculate distance from point to 3D line
double calculatePointToLineDistance(double px, double py, double pz, 
                                   double x0, double y0, double z0,
                                   double vx, double vy, double vz) {
    // Vector from line point to target point
    double dx = px - x0;
    double dy = py - y0;
    double dz = pz - z0;
    
    // Cross product: (P - P0) × d
    double crossX = dy * vz - dz * vy;
    double crossY = dz * vx - dx * vz;
    double crossZ = dx * vy - dy * vx;
    
    // Distance = |cross product| / |direction vector|
    double crossMagnitude = sqrt(crossX*crossX + crossY*crossY + crossZ*crossZ);
    
    return crossMagnitude;
}

// Main entry point function that matches THIS_NAME
void DistanceAnalysis_3D() {
    gROOT->SetBatch(kTRUE);
    
    // Get command line arguments from the global application
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Check arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_3D_file.root>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /path/to/file_3D_event42_Trj_GroupResSide.root" << std::endl;
        return;
    }
    
    TString inputFile = argv[1];
    
    std::cout << "Processing file: " << inputFile << std::endl;
    
    // Open input file
    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << inputFile << std::endl;
        return;
    }
    
    // Load trajectory parameters
    TTree* trajTree = (TTree*)file->Get("Trajectory");
    if (!trajTree) {
        std::cerr << "Error: No trajectory tree found in file" << std::endl;
        file->Close();
        return;
    }
    
    Double_t x0, y0, z0, vx, vy, vz;
    Int_t valid;
    
    trajTree->SetBranchAddress("x0", &x0);
    trajTree->SetBranchAddress("y0", &y0);
    trajTree->SetBranchAddress("z0", &z0);
    trajTree->SetBranchAddress("vx", &vx);
    trajTree->SetBranchAddress("vy", &vy);
    trajTree->SetBranchAddress("vz", &vz);
    trajTree->SetBranchAddress("valid", &valid);
    
    if (trajTree->GetEntries() == 0) {
        std::cerr << "Error: No trajectory data found" << std::endl;
        file->Close();
        return;
    }
    
    trajTree->GetEntry(0);
    
    if (!valid) {
        std::cerr << "Error: Invalid trajectory data" << std::endl;
        file->Close();
        return;
    }
    
    std::cout << "Trajectory loaded: Entry=(" << x0 << "," << y0 << "," << z0 << ")" << std::endl;
    std::cout << "Direction=(" << vx << "," << vy << "," << vz << ")" << std::endl;
    
    // Load h3D_Alpha histogram which contains confidence levels
    TH3F* h3D_Alpha = (TH3F*)file->Get("h3D_Alpha");
    if (!h3D_Alpha) {
        std::cerr << "Error: Cannot find h3D_Alpha histogram" << std::endl;
        file->Close();
        return;
    }
    
    std::cout << "Processing h3D_Alpha histogram with confidence levels..." << std::endl;
    
    // Create distance distribution histogram with 1 cm bins
    TH1F* h_distances = new TH1F("h_distances", 
                                "Distance from High-Confidence Voxels to Fitted Line;Distance [cm];Number of Voxels",
                                100, 0, 5);  // 100 bins, 0-20 cm range = 0.2 cm bin size
    
    int highConfidenceCount = 0;
    int totalVoxels = 0;
    
    // Loop through all bins in h3D_Alpha
    for (int binX = 1; binX <= h3D_Alpha->GetNbinsX(); binX++) {
        for (int binY = 1; binY <= h3D_Alpha->GetNbinsY(); binY++) {
            for (int binZ = 1; binZ <= h3D_Alpha->GetNbinsZ(); binZ++) {
                
                double confidence = h3D_Alpha->GetBinContent(binX, binY, binZ);
                
                if (confidence > 0) {  // Non-empty voxel
                    totalVoxels++;
                    
                    if (confidence >= 0.5) {  // High-confidence voxel
                        // Get voxel center coordinates
                        double voxel_z = h3D_Alpha->GetXaxis()->GetBinCenter(binX);
                        double voxel_y = h3D_Alpha->GetYaxis()->GetBinCenter(binY);
                        double voxel_x = h3D_Alpha->GetZaxis()->GetBinCenter(binZ);
                        
                        // Calculate 3D distance to fitted line
                        double distance = calculatePointToLineDistance(voxel_x, voxel_y, voxel_z, 
                                                                      x0, y0, z0, vx, vy, vz);
                        
                        // Fill distance distribution
                        h_distances->Fill(distance);
                        
                        highConfidenceCount++;
                    }
                }
            }
        }
    }
    
    std::cout << "Total voxels processed: " << totalVoxels << std::endl;
    std::cout << "Distances calculated for " << highConfidenceCount << " high-confidence voxels" << std::endl;
   
     // Extract event number from filename (must be exactly "event" followed by digits)
    TString eventNumber = "";
    TString baseName = gSystem->BaseName(inputFile);
    
    // Use regex-like approach to find "event" followed by digits
    Int_t eventPos = baseName.Index("event");
    while (eventPos != kNPOS) {
        // Check if this is exactly "event" (not part of "events" or other words)
        Bool_t isExactMatch = kTRUE;
        
        // Check character before "event" (should be non-alphabetic or start of string)
        if (eventPos > 0) {
            char prevChar = baseName[eventPos - 1];
            if ((prevChar >= 'a' && prevChar <= 'z') || (prevChar >= 'A' && prevChar <= 'Z')) {
                isExactMatch = kFALSE;
            }
        }
        
        // Check character after "event" (should be non-alphabetic)
        if (isExactMatch && eventPos + 5 < baseName.Length()) {
            char nextChar = baseName[eventPos + 5];
            if ((nextChar >= 'a' && nextChar <= 'z') || (nextChar >= 'A' && nextChar <= 'Z')) {
                isExactMatch = kFALSE;
            }
        }
        
        if (isExactMatch) {
            // Extract substring starting after "event"
            TString afterEvent = baseName(eventPos + 5, baseName.Length() - eventPos - 5);
            
            // Find the end of the number (next non-digit character)
            Int_t numberEnd = 0;
            for (Int_t i = 0; i < afterEvent.Length(); i++) {
                if (afterEvent[i] >= '0' && afterEvent[i] <= '9') {
                    numberEnd = i + 1;
                } else {
                    break;
                }
            }
            
            if (numberEnd > 0) {
                eventNumber = afterEvent(0, numberEnd);
                // Pad with zeros to make it 6 digits
                while (eventNumber.Length() < 6) {
                    eventNumber = "0" + eventNumber;
                }
                break; // Found valid event number, exit loop
            }
        }
        
        // Look for next occurrence if this wasn't a match
        eventPos = baseName.Index("event", eventPos + 1);
    }

    // Create and save distance distribution plot
    TCanvas* c1 = new TCanvas("c1", "Distance Distribution", 800, 600);
    h_distances->Draw();
    h_distances->SetLineColor(kBlue);
    h_distances->SetLineWidth(2);
    h_distances->SetFillColor(kBlue);
    h_distances->SetFillStyle(3001);
    
    // Add event number to plot title if available
    if (eventNumber != "") {
        TString plotTitle = Form("Distance Distribution - Event %s", eventNumber.Data());
        h_distances->SetTitle(plotTitle);
    }
    
    // Add statistics
    double mean = h_distances->GetMean();
    double rms = h_distances->GetRMS();
    
    std::cout << "\nDistance statistics for high-confidence voxels:" << std::endl;
    std::cout << "  Mean distance: " << std::fixed << std::setprecision(2) << mean << " cm" << std::endl;
    std::cout << "  RMS: " << std::fixed << std::setprecision(2) << rms << " cm" << std::endl;
    
    // Save plot to PNG following the framework pattern
    TString pngDir = createDistanceOutputDirectory(inputFile);
    
    // Create PNG filename with event number
    TString pngFile;
    if (eventNumber != "") {
        pngFile = Form("%s/event_%s_distanceGood.png", pngDir.Data(), eventNumber.Data());
    } else {
        // Fallback to original naming if event number not found
        pngFile = Form("%s/%s_distance_distribution.png", pngDir.Data(), baseName.Data());
    }
    
    c1->SaveAs(pngFile);
    std::cout << "Distance distribution plot saved to: " << pngFile << std::endl;
    
    file->Close();
    
    std::cout << "\nDistance analysis complete!" << std::endl;
    std::cout << "PNG plot saved to: " << pngFile << std::endl;
    
    exit(0);
}