// Saves reconstructed hits to TH3 histograms for offline viewing
#define THIS_NAME EventDisplays_3D
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"

// Include the framework headers
#include "../../src/tools/global_header.hh"

using namespace std;

// Structure to hold 3D hit information
struct Hit3D {
    double x, y, z;
    double charge;
    double timeXY, timeXZ, timeZY;
    int nViews;
    std::string matchType;
};

// Structure to organize hits by view
struct ViewHits {
    std::vector<int> indices;
    std::vector<std::pair<int, int>> positions;
    std::vector<double> charges;
    std::vector<double> times;
};

// Function to create directory if it doesn't exist
TString createOutputDirectory(const TString& inputFile) {
    // Extract directory path from input file
    TString inputDir = gSystem->DirName(inputFile);
    
    // Create the subfolder path
    TString outputDir = inputDir + "/3D_Display_root";
    
    // Check if directory exists, create if not
    if (gSystem->AccessPathName(outputDir)) {
        std::cout << "Creating output directory: " << outputDir << std::endl;
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            std::cerr << "Using current directory instead." << std::endl;
            return "./3D_Display_root";
        }
    } else {
        std::cout << "Using existing output directory: " << outputDir << std::endl;
    }
    
    return outputDir;
}

// Function to reconstruct 3D hits from 2D projections
std::vector<Hit3D> reconstructHits(Event* event) {
    std::vector<Hit3D> hits3D;
    
    // Organize hits by view
    ViewHits viewXY, viewXZ, viewZY;
    
    TClonesArray* hits = event->GetHits();
    if (!hits) return hits3D;
    
    // Sort hits into views
    for (int i = 0; i < hits->GetEntries(); i++) {
        Hit* hit = (Hit*)hits->At(i);
        if (!hit) continue;
        
        if (hit->GetView() == 0) {  // XY view
            viewXY.indices.push_back(i);
            viewXY.positions.push_back({hit->GetX(), hit->GetY()});
            viewXY.charges.push_back(hit->GetPE());
            viewXY.times.push_back(hit->GetTfromSpill());
        }
        else if (hit->GetView() == 1) {  // XZ view
            viewXZ.indices.push_back(i);
            viewXZ.positions.push_back({hit->GetX(), hit->GetZ()});
            viewXZ.charges.push_back(hit->GetPE());
            viewXZ.times.push_back(hit->GetTfromSpill());
        }
        else if (hit->GetView() == 2) {  // ZY view
            viewZY.indices.push_back(i);
            viewZY.positions.push_back({hit->GetZ(), hit->GetY()});
            viewZY.charges.push_back(hit->GetPE());
            viewZY.times.push_back(hit->GetTfromSpill());
        }
    }
    
    std::cout << "Hits per view: XY=" << viewXY.indices.size() 
              << ", XZ=" << viewXZ.indices.size() 
              << ", ZY=" << viewZY.indices.size() << std::endl;
    
    // Track which hits have been used
    std::vector<bool> usedXY(viewXY.indices.size(), false);
    std::vector<bool> usedXZ(viewXZ.indices.size(), false);
    std::vector<bool> usedZY(viewZY.indices.size(), false);
    
    // First pass: Find 3-view matches
    for (size_t ixy = 0; ixy < viewXY.positions.size(); ixy++) {
        if (usedXY[ixy]) continue;
        
        int x_xy = viewXY.positions[ixy].first;
        int y_xy = viewXY.positions[ixy].second;
        
        for (size_t ixz = 0; ixz < viewXZ.positions.size(); ixz++) {
            if (usedXZ[ixz]) continue;
            
            int x_xz = viewXZ.positions[ixz].first;
            int z_xz = viewXZ.positions[ixz].second;
            
            if (x_xy != x_xz) continue;
            
            for (size_t izy = 0; izy < viewZY.positions.size(); izy++) {
                if (usedZY[izy]) continue;
                
                int z_zy = viewZY.positions[izy].first;
                int y_zy = viewZY.positions[izy].second;
                
                if (y_xy == y_zy && z_xz == z_zy) {
                    Hit3D hit3d;
                    hit3d.x = x_xy;
                    hit3d.y = y_xy;
                    hit3d.z = z_xz;
                    hit3d.charge = (viewXY.charges[ixy] + viewXZ.charges[ixz] + viewZY.charges[izy]) / 3.0;
                    hit3d.timeXY = viewXY.times[ixy];
                    hit3d.timeXZ = viewXZ.times[ixz];
                    hit3d.timeZY = viewZY.times[izy];
                    hit3d.nViews = 3;
                    hit3d.matchType = "Full";
                    
                    hits3D.push_back(hit3d);
                    usedXY[ixy] = true;
                    usedXZ[ixz] = true;
                    usedZY[izy] = true;
                    break;
                }
            }
        }
    }
    
    // Second pass: Find 2-view matches
    // XY-XZ matches
    for (size_t ixy = 0; ixy < viewXY.positions.size(); ixy++) {
        if (usedXY[ixy]) continue;
        
        int x_xy = viewXY.positions[ixy].first;
        int y_xy = viewXY.positions[ixy].second;
        
        for (size_t ixz = 0; ixz < viewXZ.positions.size(); ixz++) {
            if (usedXZ[ixz]) continue;
            
            int x_xz = viewXZ.positions[ixz].first;
            int z_xz = viewXZ.positions[ixz].second;
            
            if (x_xy == x_xz) {
                Hit3D hit3d;
                hit3d.x = x_xy;
                hit3d.y = y_xy;
                hit3d.z = z_xz;
                hit3d.charge = (viewXY.charges[ixy] + viewXZ.charges[ixz]) / 2.0;
                hit3d.timeXY = viewXY.times[ixy];
                hit3d.timeXZ = viewXZ.times[ixz];
                hit3d.timeZY = -1;
                hit3d.nViews = 2;
                hit3d.matchType = "Partial-XY-XZ";
                
                hits3D.push_back(hit3d);
                usedXY[ixy] = true;
                usedXZ[ixz] = true;
                break;
            }
        }
    }
    
    // XY-ZY matches
    for (size_t ixy = 0; ixy < viewXY.positions.size(); ixy++) {
        if (usedXY[ixy]) continue;
        
        int x_xy = viewXY.positions[ixy].first;
        int y_xy = viewXY.positions[ixy].second;
        
        for (size_t izy = 0; izy < viewZY.positions.size(); izy++) {
            if (usedZY[izy]) continue;
            
            int z_zy = viewZY.positions[izy].first;
            int y_zy = viewZY.positions[izy].second;
            
            if (y_xy == y_zy) {
                Hit3D hit3d;
                hit3d.x = x_xy;
                hit3d.y = y_xy;
                hit3d.z = z_zy;
                hit3d.charge = (viewXY.charges[ixy] + viewZY.charges[izy]) / 2.0;
                hit3d.timeXY = viewXY.times[ixy];
                hit3d.timeXZ = -1;
                hit3d.timeZY = viewZY.times[izy];
                hit3d.nViews = 2;
                hit3d.matchType = "Partial-XY-ZY";
                
                hits3D.push_back(hit3d);
                usedXY[ixy] = true;
                usedZY[izy] = true;
                break;
            }
        }
    }
    
    // XZ-ZY matches
    for (size_t ixz = 0; ixz < viewXZ.positions.size(); ixz++) {
        if (usedXZ[ixz]) continue;
        
        int x_xz = viewXZ.positions[ixz].first;
        int z_xz = viewXZ.positions[ixz].second;
        
        for (size_t izy = 0; izy < viewZY.positions.size(); izy++) {
            if (usedZY[izy]) continue;
            
            int z_zy = viewZY.positions[izy].first;
            int y_zy = viewZY.positions[izy].second;
            
            if (z_xz == z_zy) {
                Hit3D hit3d;
                hit3d.x = x_xz;
                hit3d.y = y_zy;
                hit3d.z = z_xz;
                hit3d.charge = (viewXZ.charges[ixz] + viewZY.charges[izy]) / 2.0;
                hit3d.timeXY = -1;
                hit3d.timeXZ = viewXZ.times[ixz];
                hit3d.timeZY = viewZY.times[izy];
                hit3d.nViews = 2;
                hit3d.matchType = "Partial-XZ-ZY";
                
                hits3D.push_back(hit3d);
                usedXZ[ixz] = true;
                usedZY[izy] = true;
                break;
            }
        }
    }
    
    // Count unmatched hits
    int unmatchedXY = std::count(usedXY.begin(), usedXY.end(), false);
    int unmatchedXZ = std::count(usedXZ.begin(), usedXZ.end(), false);
    int unmatchedZY = std::count(usedZY.begin(), usedZY.end(), false);
    
    std::cout << "Unmatched hits: XY=" << unmatchedXY 
              << ", XZ=" << unmatchedXZ 
              << ", ZY=" << unmatchedZY << std::endl;
    
    return hits3D;
}

// Function to print summary
void printSummary(const std::vector<Hit3D>& hits3D, int eventID) {
    std::cout << "\n==================================================\n";
    std::cout << "3D Reconstruction Summary for Event " << eventID << "\n";
    std::cout << "==================================================\n";
    
    int nFull = 0, nPartial = 0;
    std::map<std::string, int> partialCounts;
    
    for (const auto& hit : hits3D) {
        if (hit.nViews == 3) nFull++;
        else {
            nPartial++;
            partialCounts[hit.matchType]++;
        }
    }
    
    std::cout << "\nTotal 3D hits reconstructed: " << hits3D.size() << std::endl;
    std::cout << "  - Full matches (3 views): " << nFull << std::endl;
    std::cout << "  - Partial matches (2 views): " << nPartial << std::endl;
    
    if (nPartial > 0) {
        std::cout << "\nPartial match breakdown:\n";
        for (const auto& p : partialCounts) {
            std::cout << "    " << p.first << ": " << p.second << std::endl;
        }
    }
    
    std::cout << "\nDetailed hit list:\n";
    std::cout << std::setw(5) << "Hit#" 
              << std::setw(15) << "Type"
              << std::setw(8) << "X" 
              << std::setw(8) << "Y" 
              << std::setw(8) << "Z"
              << std::setw(10) << "Charge"
              << std::setw(12) << "Time_XY"
              << std::setw(12) << "Time_XZ"
              << std::setw(12) << "Time_ZY"
              << std::endl;
    std::cout << std::string(95, '-') << std::endl;
    
    int hitNum = 0;
    for (const auto& hit : hits3D) {
        std::cout << std::setw(5) << hitNum++
                  << std::setw(15) << hit.matchType
                  << std::setw(8) << hit.x
                  << std::setw(8) << hit.y
                  << std::setw(8) << hit.z
                  << std::setw(10) << std::fixed << std::setprecision(1) << hit.charge;
        
        if (hit.timeXY >= 0)
            std::cout << std::setw(12) << std::fixed << std::setprecision(1) << hit.timeXY;
        else
            std::cout << std::setw(12) << "---";
            
        if (hit.timeXZ >= 0)
            std::cout << std::setw(12) << std::fixed << std::setprecision(1) << hit.timeXZ;
        else
            std::cout << std::setw(12) << "---";
            
        if (hit.timeZY >= 0)
            std::cout << std::setw(12) << std::fixed << std::setprecision(1) << hit.timeZY;
        else
            std::cout << std::setw(12) << "---";
            
        std::cout << std::endl;
    }
}

// Main entry point function that matches THIS_NAME
void EventDisplays_3D() {
    // Get command line arguments from the global application
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Check arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <event_number>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /path/to/events.root 42" << std::endl;
        return;
    }
    
    TString inputFile = argv[1];
    int eventNumber = atoi(argv[2]);
    
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Event number: " << eventNumber << std::endl;
    
    // Open input file
    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << inputFile << std::endl;
        return;
    }
    
    // Get the TimeGroupedEvents tree
    TTree* tree = (TTree*)file->Get("TimeGroupedEvents");
    if (!tree) {
        std::cerr << "Error: Cannot find TimeGroupedEvents tree" << std::endl;
        file->Close();
        return;
    }
    
    // Check event number
    Long64_t nEvents = tree->GetEntries();
    if (eventNumber < 0 || eventNumber >= nEvents) {
        std::cerr << "Error: Event number " << eventNumber 
                  << " out of range. File contains " << nEvents << " events." << std::endl;
        file->Close();
        return;
    }
    
    // Set up event reading
    Event* event = new Event();
    tree->SetBranchAddress("Event", &event);
    
    // Get the specified event
    tree->GetEntry(eventNumber);
    
    std::cout << "\nProcessing Event " << eventNumber 
              << " (Event ID: " << event->GetEventID() << ")" << std::endl;
    std::cout << "Total 2D hits in event: " << event->GetNHits() << std::endl;
    
    // Reconstruct 3D hits
    std::vector<Hit3D> hits3D = reconstructHits(event);
    
    // Print summary
    printSummary(hits3D, event->GetEventID());
    
    // Create output directory and filename
    TString outputDir = createOutputDirectory(inputFile);
    TString inputBaseName = gSystem->BaseName(inputFile);
    inputBaseName.ReplaceAll(".root", "");
    TString outputFile = Form("%s/%s_3D_event%d.root", outputDir.Data(), inputBaseName.Data(), eventNumber);
    
    std::cout << "\nSaving output to: " << outputFile << std::endl;
    
    // Create output file
    TFile* fOutput = new TFile(outputFile, "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
        file->Close();
        return;
    }
    
    // Create TH3 histograms for different hit types
    // Note: Y axis is vertical, X and Z are horizontal
    TH3F* h3D_Full = new TH3F("h3D_Full", 
                          Form("Event %d - Full 3D Matches (3 views);X [cm];Z [cm];Y [cm]", eventNumber),
                          24, 0, 24, 48, 0, 48, 8, 0, 8);

    TH3F* h3D_Partial = new TH3F("h3D_Partial", 
                             Form("Event %d - Partial 3D Matches (2 views);X [cm];Z [cm];Y [cm]", eventNumber),
                             24, 0, 24, 48, 0, 48, 8, 0, 8);

    TH3F* h3D_All = new TH3F("h3D_All", 
                         Form("Event %d - All 3D Matches;X [cm];Z [cm];Y [cm]", eventNumber),
                         24, 0, 24, 48, 0, 48, 8, 0, 8);

    // Create a separate histogram for transparency/alpha values
    TH3F* h3D_Alpha = new TH3F("h3D_Alpha", 
                           Form("Event %d - Alpha Values;X [cm];Z [cm];Y [cm]", eventNumber),
                           24, 0, 24, 48, 0, 48, 8, 0, 8);

    // Fill histograms with charge as weight and set transparency
    for (const auto& hit : hits3D) {
    // Get the bin centers - fill entire bin by using bin center coordinates
    int binX = h3D_All->GetXaxis()->FindBin(hit.x + 0.5);
    int binZ = h3D_All->GetYaxis()->FindBin(hit.z + 0.5);  // Z maps to Y axis in ROOT
    int binY = h3D_All->GetZaxis()->FindBin(hit.y + 0.5);  // Y maps to Z axis in ROOT
    
    // Set the bin content to the charge value (for color scaling)
    // Note: ROOT TH3 uses (X, Y, Z) ordering, but we map as (X, Z, Y)
    h3D_All->SetBinContent(binX, binZ, binY, hit.charge);
    
    // Set alpha value based on hit type
    double alphaValue;
    if (hit.nViews == 3) {
        alphaValue = 1.0;  // Full opacity for 3-view matches
        h3D_Full->SetBinContent(binX, binZ, binY, hit.charge);
    } else {
        alphaValue = 0.5;  // Semi-transparent for 2-view matches
        h3D_Partial->SetBinContent(binX, binZ, binY, hit.charge);
    }
    
    h3D_Alpha->SetBinContent(binX, binZ, binY, alphaValue);
    }

    // Set up color palette for charge visualization
    gStyle->SetPalette(kViridis);  // or use kViridis, kPlasma, etc.

    // Configure histogram display properties
    h3D_All->SetOption("BOX2Z");
    h3D_Full->SetOption("BOX2Z");
    h3D_Partial->SetOption("BOX2Z");

    // Set minimum to 0 for better color scaling
    h3D_All->SetMinimum(0);
    h3D_Full->SetMinimum(0);
    h3D_Partial->SetMinimum(0);

    // Write histograms
    h3D_Full->Write();
    h3D_Partial->Write();
    h3D_All->Write();
    h3D_Alpha->Write();
    
    // Also save a tree with the 3D hit information
    TTree* tree3D = new TTree("Hits3D", "Reconstructed 3D Hits");
    
    // Variables for the tree
    Double_t t_x, t_y, t_z, t_charge;
    Double_t t_timeXY, t_timeXZ, t_timeZY;
    Int_t t_nViews;
    Char_t t_matchType[50];
    
    tree3D->Branch("x", &t_x, "x/D");
    tree3D->Branch("y", &t_y, "y/D");
    tree3D->Branch("z", &t_z, "z/D");
    tree3D->Branch("charge", &t_charge, "charge/D");
    tree3D->Branch("timeXY", &t_timeXY, "timeXY/D");
    tree3D->Branch("timeXZ", &t_timeXZ, "timeXZ/D");
    tree3D->Branch("timeZY", &t_timeZY, "timeZY/D");
    tree3D->Branch("nViews", &t_nViews, "nViews/I");
    tree3D->Branch("matchType", t_matchType, "matchType/C");
    
    // Fill the tree
    for (const auto& hit : hits3D) {
        t_x = hit.x;
        t_y = hit.y;
        t_z = hit.z;
        t_charge = hit.charge;
        t_timeXY = hit.timeXY;
        t_timeXZ = hit.timeXZ;
        t_timeZY = hit.timeZY;
        t_nViews = hit.nViews;
        strcpy(t_matchType, hit.matchType.c_str());
        
        tree3D->Fill();
    }
    
    tree3D->Write();
    
    // Save summary information
    TNamed eventInfo("EventInfo", Form("Event %d: %zu 3D hits reconstructed", eventNumber, hits3D.size()));
    eventInfo.Write();
    
    fOutput->Close();
    file->Close();
    
    std::cout << "\n3D reconstruction complete!" << std::endl;
    std::cout << "Output saved to: " << outputFile << std::endl;
    std::cout << "\nTo view in ROOT:" << std::endl;
    std::cout << "  root " << outputFile << std::endl;
    std::cout << "  h3D_All->Draw(\"BOX2Z\")" << std::endl;
    std::cout << "  // For color-coded charge visualization:" << std::endl;
    std::cout << "  gStyle->SetPalette(kViridis)" << std::endl;
    std::cout << "  h3D_All->SetMinimum(0)" << std::endl;
    std::cout << "  h3D_All->Draw(\"BOX2Z\")" << std::endl;
    std::cout << "  // For transparency effects, you may need to combine with h3D_Alpha" << std::endl;
    std::cout << "  or" << std::endl;
    std::cout << "  h3D_Full->Draw(\"BOX2Z\")" << std::endl; exit(1);
    
    delete event;
}