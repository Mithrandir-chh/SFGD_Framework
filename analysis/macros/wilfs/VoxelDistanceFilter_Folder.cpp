// Calculates distances from high-confidence voxels to fitted trajectory line
#define THIS_NAME DistanceAnalysis_3D
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
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
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

// Include the framework headers
#include "../../src/tools/global_header.hh"

using namespace std;

// Structure to hold event statistics
struct EventStats {
    TString eventNumber;
    TString fileName;
    double mean;
    double stdev;
    int highConfidenceCount;
    int totalVoxels;
};

// Function to create directory if it doesn't exist
TString createDistanceOutputDirectory(const TString& inputFolder) {
    // Go one level up from input folder
    TString parentDir = gSystem->DirName(inputFolder);
    
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

// Function to extract event number from filename
TString extractEventNumber(const TString& fileName) {
    TString eventNumber = "";
    TString baseName = gSystem->BaseName(fileName);
    
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
    
    return eventNumber;
}

// Function to process a single file and return statistics
EventStats processSingleFile(const TString& inputFile, const TString& outputDir) {
    EventStats stats;
    stats.fileName = gSystem->BaseName(inputFile);
    stats.eventNumber = extractEventNumber(inputFile);
    stats.mean = 0;
    stats.stdev = 0;
    stats.highConfidenceCount = 0;
    stats.totalVoxels = 0;
    
    std::cout << "Processing file: " << stats.fileName << std::endl;
    
    // Open input file
    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << inputFile << std::endl;
        return stats;
    }
    
    // Load trajectory parameters
    TTree* trajTree = (TTree*)file->Get("Trajectory");
    if (!trajTree) {
        std::cerr << "Error: No trajectory tree found in file" << std::endl;
        file->Close();
        return stats;
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
    
    if (trajTree->GetEntries() == 0 || !valid) {
        std::cerr << "Error: Invalid trajectory data in " << stats.fileName << std::endl;
        file->Close();
        return stats;
    }
    
    trajTree->GetEntry(0);
    
    // Load h3D_Alpha histogram which contains confidence levels
    TH3F* h3D_Alpha = (TH3F*)file->Get("h3D_Alpha");
    if (!h3D_Alpha) {
        std::cerr << "Error: Cannot find h3D_Alpha histogram in " << stats.fileName << std::endl;
        file->Close();
        return stats;
    }
    
    // Create distance distribution histogram
    TH1F* h_distances = new TH1F("h_distances", 
                                "Distance from High-Confidence Voxels to Fitted Line;Distance [cm];Number of Voxels",
                                100, 0, 5);
    
    // Loop through all bins in h3D_Alpha
    for (int binX = 1; binX <= h3D_Alpha->GetNbinsX(); binX++) {
        for (int binY = 1; binY <= h3D_Alpha->GetNbinsY(); binY++) {
            for (int binZ = 1; binZ <= h3D_Alpha->GetNbinsZ(); binZ++) {
                
                double confidence = h3D_Alpha->GetBinContent(binX, binY, binZ);
                
                if (confidence > 0) {  // Non-empty voxel
                    stats.totalVoxels++;
                    
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
                        
                        stats.highConfidenceCount++;
                    }
                }
            }
        }
    }
    
    // Calculate statistics
    stats.mean = h_distances->GetMean();
    stats.stdev = h_distances->GetStdDev();
    
    // Create and save individual plot
    TCanvas* c1 = new TCanvas("c1", "Distance Distribution", 800, 600);
    h_distances->Draw();
    h_distances->SetLineColor(kBlue);
    h_distances->SetLineWidth(2);
    h_distances->SetFillColor(kBlue);
    h_distances->SetFillStyle(3001);
    
    // Add event number to plot title if available
    if (stats.eventNumber != "") {
        TString plotTitle = Form("Distance Distribution - Event %s", stats.eventNumber.Data());
        h_distances->SetTitle(plotTitle);
    }
    
    // Save individual PNG
    TString pngFile;
    if (stats.eventNumber != "") {
        pngFile = Form("%s/event_%s_distanceGood.png", outputDir.Data(), stats.eventNumber.Data());
    } else {
        pngFile = Form("%s/%s_distance_distribution.png", outputDir.Data(), stats.fileName.Data());
    }
    
    c1->SaveAs(pngFile);
    std::cout << "  Individual plot saved: " << pngFile << std::endl;
    
    // Cleanup
    delete c1;
    delete h_distances;
    file->Close();
    
    std::cout << "  Event " << stats.eventNumber << ": Mean=" << std::fixed << std::setprecision(2) 
              << stats.mean << " cm, StDev=" << stats.stdev << " cm" << std::endl;
    
    return stats;
}

// Main entry point function that matches THIS_NAME
void DistanceAnalysis_3D() {
    gROOT->SetBatch(kTRUE);
    
    // Get command line arguments from the global application
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Check arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_folder>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /path/to/root_files/" << std::endl;
        return;
    }
    
    TString inputFolder = argv[1];
    
    // Check if input folder exists
    if (gSystem->AccessPathName(inputFolder)) {
        std::cerr << "Error: Input folder " << inputFolder << " not found" << std::endl;
        return;
    }
    
    std::cout << "Processing all ROOT files from folder: " << inputFolder << std::endl;
    
    // Create output directory
    TString outputDir = createDistanceOutputDirectory(inputFolder);
    
    // Find all .root files in the input folder
    TSystemDirectory dir(inputFolder, inputFolder);
    TList* files = dir.GetListOfFiles();
    
    vector<TString> rootFiles;
    if (files) {
        TSystemFile* file;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            TString fileName = file->GetName();
            if (!file->IsDirectory() && fileName.EndsWith(".root")) {
                rootFiles.push_back(inputFolder + "/" + fileName);
            }
        }
    }
    
    if (rootFiles.empty()) {
        std::cerr << "Error: No .root files found in " << inputFolder << std::endl;
        return;
    }
    
    std::cout << "Found " << rootFiles.size() << " ROOT files to process" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Process all files
    vector<EventStats> allStats;
    vector<double> allMeans;
    vector<double> allStDevs;
    
    for (size_t i = 0; i < rootFiles.size(); i++) {
        std::cout << "[" << (i+1) << "/" << rootFiles.size() << "] ";
        EventStats stats = processSingleFile(rootFiles[i], outputDir);
        
        if (stats.highConfidenceCount > 0) {  // Only include valid results
            allStats.push_back(stats);
            allMeans.push_back(stats.mean);
            allStDevs.push_back(stats.stdev);
        }
        std::cout << "---" << std::endl;
    }
    
    if (allStats.empty()) {
        std::cerr << "Error: No valid events processed" << std::endl;
        return;
    }
    
    // Calculate summary statistics
    double sumMeans = 0, sumStDevs = 0;
    for (size_t i = 0; i < allMeans.size(); i++) {
        sumMeans += allMeans[i];
        sumStDevs += allStDevs[i];
    }
    
    double overallMeanOfMeans = sumMeans / allMeans.size();
    double overallMeanOfStDevs = sumStDevs / allStDevs.size();
    
    // Calculate standard deviation of means and StDevs
    double sumSquaredDiffMeans = 0, sumSquaredDiffStDevs = 0;
    for (size_t i = 0; i < allMeans.size(); i++) {
        sumSquaredDiffMeans += (allMeans[i] - overallMeanOfMeans) * (allMeans[i] - overallMeanOfMeans);
        sumSquaredDiffStDevs += (allStDevs[i] - overallMeanOfStDevs) * (allStDevs[i] - overallMeanOfStDevs);
    }
    
    double stdMeans = sqrt(sumSquaredDiffMeans / allMeans.size());
    double stdStDevs = sqrt(sumSquaredDiffStDevs / allStDevs.size());
    
    // Create summary histograms
    TH1F* h_means = new TH1F("h_means", "Distribution of Mean Distances;Mean Distance [cm];Number of Events", 
                            50, overallMeanOfMeans - 3*stdMeans, overallMeanOfMeans + 3*stdMeans);
    TH1F* h_stdev_dist = new TH1F("h_stdev_dist", "Distribution of StDev Values;StDev [cm];Number of Events", 
                                 50, overallMeanOfStDevs - 3*stdStDevs, overallMeanOfStDevs + 3*stdStDevs);
    
    // Fill summary histograms
    for (size_t i = 0; i < allMeans.size(); i++) {
        h_means->Fill(allMeans[i]);
        h_stdev_dist->Fill(allStDevs[i]);
    }
    
    // Create summary ROOT file
    TString summaryRootFile = outputDir + "/distance_analysis_summary.root";
    TFile* fSummary = new TFile(summaryRootFile, "RECREATE");
    
    // Create summary tree
    TTree* summaryTree = new TTree("EventSummary", "Distance Analysis Summary");
    TString eventNum, fileName;
    Double_t mean, stdev;
    Int_t highConfCount, totalVox;
    
    summaryTree->Branch("EventNumber", &eventNum);
    summaryTree->Branch("FileName", &fileName);
    summaryTree->Branch("Mean", &mean);
    summaryTree->Branch("StDev", &stdev);
    summaryTree->Branch("HighConfidenceCount", &highConfCount);
    summaryTree->Branch("TotalVoxels", &totalVox);
    
    // Fill summary tree
    for (size_t i = 0; i < allStats.size(); i++) {
        eventNum = allStats[i].eventNumber;
        fileName = allStats[i].fileName;
        mean = allStats[i].mean;
        stdev = allStats[i].stdev;
        highConfCount = allStats[i].highConfidenceCount;
        totalVox = allStats[i].totalVoxels;
        summaryTree->Fill();
    }
    
    // Write summary histograms and tree
    h_means->Write();
    h_stdev_dist->Write();
    summaryTree->Write();
    
    // Create summary info
    TNamed summaryInfo("SummaryStats", 
                      Form("Processed %d events, Mean of means: %.3f±%.3f cm, Mean of StDevs: %.3f±%.3f cm", 
                           (int)allStats.size(), overallMeanOfMeans, stdMeans, overallMeanOfStDevs, stdStDevs));
    summaryInfo.Write();
    
    fSummary->Close();
    
    // Create summary plots
    TCanvas* c_summary = new TCanvas("c_summary", "Summary Statistics", 1600, 600);
    c_summary->Divide(2, 1);
    
    c_summary->cd(1);
    h_means->Draw();
    h_means->SetLineColor(kRed);
    h_means->SetLineWidth(2);
    h_means->SetFillColor(kRed);
    h_means->SetFillStyle(3001);
    
    c_summary->cd(2);
    h_stdev_dist->Draw();
    h_stdev_dist->SetLineColor(kGreen);
    h_stdev_dist->SetLineWidth(2);
    h_stdev_dist->SetFillColor(kGreen);
    h_stdev_dist->SetFillStyle(3001);
    
    // Save summary PNG
    TString summaryPngFile = outputDir + "/distance_analysis_summary.png";
    c_summary->SaveAs(summaryPngFile);
    
    // Print final summary
    std::cout << "========================================" << std::endl;
    std::cout << "BATCH PROCESSING COMPLETE!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Processed events: " << allStats.size() << std::endl;
    std::cout << "Mean of mean distances: " << std::fixed << std::setprecision(3) 
              << overallMeanOfMeans << " ± " << stdMeans << " cm" << std::endl;
    std::cout << "Mean of StDev values: " << overallMeanOfStDevs << " ± " << stdStDevs << " cm" << std::endl;
    std::cout << std::endl;
    std::cout << "Output files saved:" << std::endl;
    std::cout << "  Summary ROOT file: " << summaryRootFile << std::endl;
    std::cout << "  Summary PNG plot: " << summaryPngFile << std::endl;
    std::cout << "  Individual event PNGs: " << outputDir << "/event_XXXXXX_distanceGood.png" << std::endl;
    
    // Cleanup
    delete c_summary;
    delete h_means;
    delete h_stdev_dist;
    
    exit(0);
}