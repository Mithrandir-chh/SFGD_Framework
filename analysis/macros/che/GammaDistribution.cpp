// Charge Analysis for Decay Electrons
// Analyzes all ROOT files in a folder and creates histogram of summed charges
#define THIS_NAME ChargeAnalysis_DecayElectrons
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <chrono>
#include <glob.h>
#include <dirent.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TNamed.h"

// Include the framework headers
#include "../../src/tools/global_header.hh"

using namespace std;

// Simple storage for data
vector<double> zValues_x23p5;
vector<int> hitCounts;

// Function to get all .root files in a directory
vector<string> getROOTFiles(const string& directory) {
    vector<string> files;
    
    DIR* dir = opendir(directory.c_str());
    if (dir == nullptr) {
        cout << "Error: Cannot open directory " << directory << endl;
        return files;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;
        if (filename.length() > 5 && filename.substr(filename.length() - 5) == ".root") {
            files.push_back(directory + "/" + filename);
        }
    }
    closedir(dir);
    
    sort(files.begin(), files.end());
    return files;
}

// Function to process a single ROOT file and return sum of decay electron charges
double processSingleFile(const string& filename) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cout << "Warning: Cannot open file " << filename << endl;
        return -999.0;
    }
    
    TTree* tree = (TTree*)file->Get("Hits3D");
    if (!tree) {
        cout << "Warning: Hits3D tree not found in " << filename << endl;
        file->Close();
        delete file;
        return -999.0;
    }
    
    Double_t confidence, charge, x, z;
    
    if (!tree->GetBranch("confidence") || !tree->GetBranch("charge")) {
        cout << "Warning: Required branches not found in " << filename << endl;
        file->Close();
        delete file;
        return -999.0;
    }
    
    tree->SetBranchAddress("confidence", &confidence);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("z", &z);  
    
    Long64_t nEntries = tree->GetEntries();
    double totalCharge = 0.0;
    int validHits = 0;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        if (confidence >= 0.5) {
            totalCharge += charge;
            validHits++;
        }
        
        if (x == 23.5) {
            zValues_x23p5.push_back(z);
        }
    }
    
    hitCounts.push_back(validHits);
    
    cout << "File: " << gSystem->BaseName(filename.c_str()) 
         << " - Hits: " << validHits 
         << ", Charge: " << fixed << setprecision(2) << totalCharge << endl;
    
    file->Close();
    delete file;
    
    return totalCharge;

}

// Simple z distribution plot
void plotZDistribution(const TString& outputDir) {
    if (zValues_x23p5.empty()) {
        cout << "No x=23.5 voxels found" << endl;
        return;
    }
    
    cout << "Creating z distribution for " << zValues_x23p5.size() << " voxels" << endl;
    
    double minZ = 0;
    double maxZ = 50;
    
    TH1D* h_z = new TH1D("h_z", "Gamma Z Distribution (x=23.5);Z;Count", 50, minZ, maxZ);
    
    for (double z : zValues_x23p5) {
        h_z->Fill(z);
    }
    
    TFile* file = new TFile(outputDir + "/z_distribution.root", "RECREATE");
    h_z->Write();
    file->Close();
    delete file;
    
    TCanvas* c = new TCanvas("c_z", "Z Distribution", 800, 600);
    h_z->SetFillColor(kRed-9);
    h_z->Draw("HIST");
    c->SaveAs(outputDir + "/z_distribution.png");
    
    delete c;
    delete h_z;
}

// NEW: Simple hits distribution plot
void plotHitsDistribution(const TString& outputDir) {
    if (hitCounts.empty()) {
        cout << "No hit data found" << endl;
        return;
    }
    
    cout << "Creating hits distribution for " << hitCounts.size() << " files" << endl;
    
    int minHits = 0;
    int maxHits = 50;
    
    TH1D* h_hits = new TH1D("h_hits", "Voxel per Gamma Event;Number of Voxel;Number of Events", 
                            50, minHits, maxHits);
    
    for (int hits : hitCounts) {
        h_hits->Fill(hits);
    }
    
    TFile* file = new TFile(outputDir + "/voxel_distribution.root", "RECREATE");
    h_hits->Write();
    file->Close();
    delete file;
    
    TCanvas* c = new TCanvas("c_hits", "Voxel Distribution", 800, 600);
    h_hits->SetFillColor(kBlue-9);
    h_hits->Draw("HIST");
    c->SaveAs(outputDir + "/voxel_distribution.png");
    
    cout << "Voxel stats - Mean: " << h_hits->GetMean() 
         << ", RMS: " << h_hits->GetRMS() << endl;
    
    delete c;
    delete h_hits;
}

TString createChargeAnalysisOutputDirectory(const TString& inputDir) {
    TString outputDir = inputDir + "/charge_analysis_results";
    
    if (gSystem->AccessPathName(outputDir)) {
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            return "./charge_analysis_results";
        }
    }
    
    return outputDir;
}

void ChargeAnalysis_DecayElectrons() {
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <directory_with_root_files>" << std::endl;
        return;
    }
    
    string inputDirectory = argv[1];
    
    cout << "=== Decay Electron Charge Analysis ===" << endl;
    cout << "Input directory: " << inputDirectory << endl;
    
    if (gSystem->AccessPathName(inputDirectory.c_str())) {
        cout << "Error: Directory does not exist: " << inputDirectory << endl;
        return;
    }
    
    vector<string> rootFiles = getROOTFiles(inputDirectory);
    
    if (rootFiles.empty()) {
        cout << "No ROOT files found in directory: " << inputDirectory << endl;
        return;
    }
    
    cout << "Found " << rootFiles.size() << " ROOT files" << endl;
    
    // Clear data
    zValues_x23p5.clear();
    hitCounts.clear();
    
    vector<double> chargeSums;
    
    for (const string& file : rootFiles) {
        double chargeSum = processSingleFile(file);
        if (chargeSum != -999.0) {
            chargeSums.push_back(chargeSum);
        }
    }
    
    if (chargeSums.empty()) {
        cout << "No valid files processed!" << endl;
        return;
    }
    
    cout << "Processed " << chargeSums.size() << " valid files" << endl;
    
    // Create charge histogram
    TH1D* h_charges = new TH1D("h_charges", 
                               "Gamma energy;Total energy;Events",
                               50, 0, 10);
    
    for (double charge : chargeSums) {
        h_charges->Fill(charge);
    }
    
    cout << "Charge stats - Mean: " << h_charges->GetMean() 
         << ", RMS: " << h_charges->GetRMS() << endl;
    
    TString outputDir = createChargeAnalysisOutputDirectory(inputDirectory);
    
    // Save charge histogram
    TFile* outputFile = new TFile(outputDir + "/charge_distribution.root", "RECREATE");
    h_charges->Write();
    outputFile->Close();
    delete outputFile;
    
    gROOT->SetBatch(kTRUE);
    TCanvas* canvas = new TCanvas("c1", "Charge Distribution", 800, 600);
    h_charges->SetFillColor(kAzure-9);
    h_charges->Draw("HIST");
    canvas->SaveAs(outputDir + "/charge_distribution.png");
    
    delete canvas;
    delete h_charges;
    
    // Create additional plots
    plotZDistribution(outputDir);
    plotHitsDistribution(outputDir);
    
    cout << "Analysis complete. Output in: " << outputDir << endl;
    
    exit(0);
}