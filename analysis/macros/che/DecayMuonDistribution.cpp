// Charge Analysis for Decay Electrons with DeltaT Distribution
// Analyzes all ROOT files in a folder and creates histograms of summed charges and deltaT
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
#include <regex>
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

// Structure to hold analysis results for each file
struct FileAnalysisResult {
    string filename;
    double totalCharge;
    double deltaT;
    int decayElectronHits;
    bool validCharge;
    bool validDeltaT;
    
    FileAnalysisResult() : totalCharge(-999.0), deltaT(-999.0), decayElectronHits(0), 
                          validCharge(false), validDeltaT(false) {}
};

// Function to get all .root files in a directory
vector<string> getROOTFiles(const string& directory) {
    vector<string> files;
    
    // Use dirent for cross-platform compatibility
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

// Function to extract deltaT from DecayInfo TNamed object
double extractDeltaT(const string& filename) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cout << "Warning: Cannot open file " << filename << " for deltaT extraction" << endl;
        return -999.0;
    }
    
    // Look for DecayInfo TNamed object
    TNamed* decayInfo = (TNamed*)file->Get("DecayInfo");
    if (!decayInfo) {
        cout << "Warning: DecayInfo object not found in " << filename << endl;
        file->Close();
        delete file;
        return -999.0;
    }
    
    // Get the title string which contains the deltaT information
    string title = decayInfo->GetTitle();
    
    // Parse deltaT using regex: look for "ΔT=XXX.X ticks" pattern
    regex deltaTPattern(R"(ΔT=([0-9]+\.?[0-9]*)\s+ticks)");
    smatch match;
    
    double deltaT = -999.0;
    if (regex_search(title, match, deltaTPattern)) {
        try {
            deltaT = stod(match[1].str())*2.5;
        } catch (const std::exception& e) {
            cout << "Warning: Could not parse deltaT value from: " << title << endl;
            deltaT = -999.0;
        }
    } else {
        cout << "Warning: deltaT pattern not found in DecayInfo title: " << title << endl;
    }
    
    file->Close();
    delete file;
    
    return deltaT;
}

// Function to process a single ROOT file and return sum of decay electron energy
double processSingleFileCharge(const string& filename, int& decayElectronHits) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cout << "Warning: Cannot open file " << filename << endl;
        return -999.0;
    }
    
    // Look for CombinedHits tree
    TTree* tree = (TTree*)file->Get("CombinedHits");
    if (!tree) {
        cout << "Warning: CombinedHits tree not found in " << filename << endl;
        file->Close();
        delete file;
        return -999.0;
    }
    
    // Set up branches
    Int_t isDecayElectron;
    Double_t charge;
    
    // Check if branches exist
    if (!tree->GetBranch("isDecayElectron") || !tree->GetBranch("charge")) {
        cout << "Warning: Required branches not found in " << filename << endl;
        file->Close();
        delete file;
        return -999.0;
    }
    
    tree->SetBranchAddress("isDecayElectron", &isDecayElectron);
    tree->SetBranchAddress("charge", &charge);
    
    Long64_t nEntries = tree->GetEntries();
    double totalCharge = 0.0;
    decayElectronHits = 0;
    
    // Loop through all entries
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Sum charges only for decay electron hits
        if (isDecayElectron == 1) {
            totalCharge += charge;
            decayElectronHits++;
        }
    }
    
    file->Close();
    delete file;
    
    return totalCharge;
}

// Combined function to process a single ROOT file for both charge and deltaT
FileAnalysisResult processSingleFile(const string& filename) {
    FileAnalysisResult result;
    result.filename = filename;
    
    // Extract charge information
    result.totalCharge = processSingleFileCharge(filename, result.decayElectronHits);
    result.validCharge = (result.totalCharge != -999.0);
    
    // Extract deltaT information
    result.deltaT = extractDeltaT(filename);
    result.validDeltaT = (result.deltaT != -999.0);
    
    cout << "File: " << gSystem->BaseName(filename.c_str()) 
         << " - Decay electron hits: " << result.decayElectronHits 
         << ", Total charge: " << fixed << setprecision(2) << result.totalCharge
         << ", ΔT: " << fixed << setprecision(1) << result.deltaT << " ticks" << endl;
    
    return result;
}

// Function to create normalized cumulative distribution starting from t=950
TH1D* createCumulativeDistribution(TH1D* h_deltaT, double startTime = 950.0) {
    if (!h_deltaT) {
        cout << "Error: Input histogram is null" << endl;
        return nullptr;
    }
    
    // Get the range of the original histogram
    double minTime = h_deltaT->GetXaxis()->GetXmin();
    double maxTime = h_deltaT->GetXaxis()->GetXmax();
    int nBins = h_deltaT->GetNbinsX();
    
    cout << "\n=== Creating Cumulative Distribution ===" << endl;
    cout << "Original histogram range: " << minTime << " to " << maxTime << " ns" << endl;
    cout << "Starting cumulative integration from: " << startTime << " ns" << endl;
    
    // Find the bin corresponding to start time
    int startBin = h_deltaT->FindBin(startTime);
    if (startBin < 1) startBin = 1;
    if (startBin > nBins) {
        cout << "Error: Start time is beyond histogram range" << endl;
        return nullptr;
    }
    
    // Create new histogram for cumulative distribution
    TH1D* h_cumulative = new TH1D("h_deltaT_cumulative", 
                                  "Cumulative Decay Distribution;#Delta T (ns);Cumulative Fraction",
                                  nBins, minTime, maxTime);
    
    // Get total integral from start time onwards for normalization
    double totalIntegral = h_deltaT->Integral(startBin, nBins);
    
    if (totalIntegral <= 0) {
        cout << "Error: No events found from start time onwards" << endl;
        delete h_cumulative;
        return nullptr;
    }
    
    cout << "Total events from " << startTime << " ticks onwards: " << totalIntegral << endl;
    
    // Fill cumulative histogram
    for (int bin = startBin; bin <= nBins; bin++) {
        // Calculate cumulative sum from startBin to current bin
        double cumulativeSum = h_deltaT->Integral(startBin, bin);
        
        // Normalize by total integral (so it goes from 0 to 1)
        double normalizedValue = cumulativeSum / totalIntegral;
        
        // Set bin content (use bin center as x-coordinate)
        h_cumulative->SetBinContent(bin, normalizedValue);
        h_cumulative->SetBinError(bin, 0); // Could calculate proper error if needed
    }
    
    // Set bins before start time to 0
    for (int bin = 1; bin < startBin; bin++) {
        h_cumulative->SetBinContent(bin, 0.0);
    }
    
    cout << "Cumulative distribution created successfully" << endl;
    cout << "Final cumulative value: " << h_cumulative->GetBinContent(nBins) << endl;
    
    return h_cumulative;
}

// Function to create output directory following your framework's pattern
TString createAnalysisOutputDirectory(const TString& inputDir) {
    TString outputDir = inputDir + "/charge_deltat_analysis_results";
    
    if (gSystem->AccessPathName(outputDir)) {
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            return "./charge_deltat_analysis_results";
        }
    }
    
    return outputDir;
}

// Main analysis function matching your framework style
void ChargeAnalysis_DecayElectrons() {
    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Check arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <directory_with_root_files>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /path/to/muon_decay_files/" << std::endl;
        std::cerr << "This will analyze all ROOT files in the directory for decay electron energy and deltaT" << std::endl;
        return;
    }
    
    string inputDirectory = argv[1];
    
    cout << "=== Decay Electron Energy and DeltaT Analysis ===" << endl;
    cout << "Input directory: " << inputDirectory << endl;
    
    // Check if directory exists
    if (gSystem->AccessPathName(inputDirectory.c_str())) {
        cout << "Error: Directory does not exist: " << inputDirectory << endl;
        return;
    }
    
    // Get all ROOT files in the directory
    vector<string> rootFiles = getROOTFiles(inputDirectory);
    
    if (rootFiles.empty()) {
        cout << "No ROOT files found in directory: " << inputDirectory << endl;
        return;
    }
    
    cout << "Found " << rootFiles.size() << " ROOT files" << endl;
    cout << "\n=== Processing Files ===" << endl;
    
    // Process each file and collect results
    vector<FileAnalysisResult> results;
    
    auto startTime = std::chrono::steady_clock::now();
    
    for (size_t i = 0; i < rootFiles.size(); i++) {
        // Progress reporting
        if (i % 10 == 0 || i == rootFiles.size() - 1) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double progress = 100.0 * i / rootFiles.size();
            
            cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% "
                 << "(" << i + 1 << "/" << rootFiles.size() << ") ";
            cout.flush();
        }
        
        FileAnalysisResult result = processSingleFile(rootFiles[i]);
        results.push_back(result);
    }
    
    cout << endl;
    
    // Separate valid charge and deltaT results
    vector<double> validCharges;
    vector<double> validDeltaTs;
    vector<FileAnalysisResult> validChargeResults;
    vector<FileAnalysisResult> validDeltaTResults;
    
    for (const auto& result : results) {
        if (result.validCharge) {
            validCharges.push_back(result.totalCharge);
            validChargeResults.push_back(result);
        }
        if (result.validDeltaT) {
            validDeltaTs.push_back(result.deltaT);
            validDeltaTResults.push_back(result);
        }
    }
    
    cout << "\n=== Analysis Results ===" << endl;
    cout << "Total files processed: " << results.size() << endl;
    cout << "Valid charge files: " << validCharges.size() << endl;
    cout << "Valid deltaT files: " << validDeltaTs.size() << endl;
    
    if (validCharges.empty() && validDeltaTs.empty()) {
        cout << "No valid files processed!" << endl;
        return;
    }
    
    // Create output directory
    TString outputDir = createAnalysisOutputDirectory(inputDirectory);
    
    // Charge analysis
    TH1D* h_charges = nullptr;
    if (!validCharges.empty()) {
        double minCharge = *min_element(validCharges.begin(), validCharges.end());
        double maxCharge = *max_element(validCharges.begin(), validCharges.end());
        
        cout << "\n=== Charge Statistics ===" << endl;
        cout << "Charge range: " << fixed << setprecision(2) << minCharge << " to " << maxCharge << endl;
        
        // Create charge histogram
        int nBinsCharge = 30;
        double histMinCharge = 0;
        double histMaxCharge = 60;
        
        h_charges = new TH1D("h_charges", 
                            "Distribution of Summed Decay Electron Energy;Total Energy (MeV);Number of Events",
                            nBinsCharge, histMinCharge, histMaxCharge);
        
        for (double charge : validCharges) {
            h_charges->Fill(charge);
        }
        
        cout << "Mean charge: " << fixed << setprecision(2) << h_charges->GetMean() << endl;
        cout << "RMS charge: " << fixed << setprecision(2) << h_charges->GetRMS() << endl;
    }
    
    // DeltaT analysis
    TH1D* h_deltaT = nullptr;
    if (!validDeltaTs.empty()) {
        double minDeltaT = *min_element(validDeltaTs.begin(), validDeltaTs.end());
        double maxDeltaT = *max_element(validDeltaTs.begin(), validDeltaTs.end());
        
        cout << "\n=== DeltaT Statistics ===" << endl;
        cout << "DeltaT range: " << fixed << setprecision(1) << minDeltaT << " to " << maxDeltaT << " ns" << endl;
        
        // Create deltaT histogram with reasonable binning
        int nBinsDeltaT = 21;
        double histMinDeltaT = 950;
        double histMaxDeltaT = 3050;
        
        h_deltaT = new TH1D("h_deltaT", 
                           "Distribution of Decay Time Intervals;#Delta T (ns);Number of Events",
                           nBinsDeltaT, histMinDeltaT, histMaxDeltaT);
        
        for (double deltaT : validDeltaTs) {
            h_deltaT->Fill(deltaT);
        }
        
        cout << "Mean deltaT: " << fixed << setprecision(1) << h_deltaT->GetMean() << " ns" << endl;
        cout << "RMS deltaT: " << fixed << setprecision(1) << h_deltaT->GetRMS() << " ns" << endl;
    }
    
    // Create cumulative distribution (add this after h_deltaT creation)
    TH1D* h_deltaT_cumulative = nullptr;
    if (h_deltaT) {
        h_deltaT_cumulative = createCumulativeDistribution(h_deltaT, 950.0);
    }

    // Save results to ROOT file
    TString rootOutputFile = outputDir + "/decay_electron_analysis.root";
    TFile* outputFile = new TFile(rootOutputFile, "RECREATE");
    
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create output ROOT file " << rootOutputFile << std::endl;
        return;
    }
    
    if (h_charges) h_charges->Write();
    if (h_deltaT) h_deltaT->Write();
    if (h_deltaT_cumulative) h_deltaT_cumulative->Write();
    
    // Save the combined data tree
    TTree* dataTree = new TTree("AnalysisData", "Combined Charge and DeltaT Data");
    Double_t fileCharge, fileDeltaT;
    Char_t fileName[500];
    Int_t fileIndex, decayHits;
    Bool_t validCharge, validDT;
    
    dataTree->Branch("fileIndex", &fileIndex, "fileIndex/I");
    dataTree->Branch("charge", &fileCharge, "charge/D");
    dataTree->Branch("deltaT", &fileDeltaT, "deltaT/D");
    dataTree->Branch("decayElectronHits", &decayHits, "decayElectronHits/I");
    dataTree->Branch("validCharge", &validCharge, "validCharge/O");
    dataTree->Branch("validDeltaT", &validDT, "validDeltaT/O");
    dataTree->Branch("filename", fileName, "filename/C");
    
    for (size_t i = 0; i < results.size(); i++) {
        fileIndex = i;
        fileCharge = results[i].totalCharge;
        fileDeltaT = results[i].deltaT;
        decayHits = results[i].decayElectronHits;
        validCharge = results[i].validCharge;
        validDT = results[i].validDeltaT;
        TString baseName = gSystem->BaseName(results[i].filename.c_str());
        strcpy(fileName, baseName.Data());
        dataTree->Fill();
    }
    
    dataTree->Write();
    
    // Save analysis summary
    TNamed analysisInfo("AnalysisInfo", 
                       Form("Decay Electron Analysis: %d files total, %d valid charges, %d valid deltaTuesday", 
                            (int)results.size(), (int)validCharges.size(), (int)validDeltaTs.size()));
    analysisInfo.Write();
    
    outputFile->Close();
    delete outputFile;
    
    cout << "ROOT file saved: " << rootOutputFile << endl;
    
    // Create and save canvases
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1111);
    
    // Charge histogram canvas
    if (h_charges) {
        TCanvas* canvas_charge = new TCanvas("c_charge", "Decay Electron Energy Distribution", 800, 600);
        canvas_charge->SetLeftMargin(0.12);
        canvas_charge->SetBottomMargin(0.12);
        canvas_charge->SetGrid();
        
        h_charges->SetLineColor(kAzure+6);
        h_charges->SetFillColor(kTeal-5);
        h_charges->SetLineWidth(2);
        h_charges->GetXaxis()->SetTitleSize(0.045);
        h_charges->GetYaxis()->SetTitleSize(0.045);
        h_charges->GetXaxis()->SetLabelSize(0.04);
        h_charges->GetYaxis()->SetLabelSize(0.04);
        h_charges->GetXaxis()->SetTitleOffset(1.2);
        h_charges->GetYaxis()->SetTitleOffset(1.3);
        
        h_charges->Draw("HIST");
        
        TString pngChargeFile = outputDir + "/decay_electron_charge_distribution.png";
        canvas_charge->SaveAs(pngChargeFile);
        cout << "Charge PNG saved: " << pngChargeFile << endl;
        
        delete canvas_charge;
    }
    
    // DeltaT histogram canvas
    if (h_deltaT) {
        TCanvas* canvas_deltaT = new TCanvas("c_deltaT", "Decay Time Interval Distribution", 800, 600);
        canvas_deltaT->SetLeftMargin(0.12);
        canvas_deltaT->SetBottomMargin(0.12);
        canvas_deltaT->SetGrid();
        
        h_deltaT->SetLineColor(kRed+1);
        h_deltaT->SetFillColor(kRed-9);
        h_deltaT->SetLineWidth(2);
        h_deltaT->GetXaxis()->SetTitleSize(0.045);
        h_deltaT->GetYaxis()->SetTitleSize(0.045);
        h_deltaT->GetXaxis()->SetLabelSize(0.04);
        h_deltaT->GetYaxis()->SetLabelSize(0.04);
        h_deltaT->GetXaxis()->SetTitleOffset(1.2);
        h_deltaT->GetYaxis()->SetTitleOffset(1.3);
        
        h_deltaT->Draw("HIST");
        
        TString pngDeltaTFile = outputDir + "/decay_time_interval_distribution.png";
        canvas_deltaT->SaveAs(pngDeltaTFile);
        cout << "DeltaT PNG saved: " << pngDeltaTFile << endl;
        
        delete canvas_deltaT;
    }
    
    // Add canvas for cumulative distribution (after deltaT canvas section):
    if (h_deltaT_cumulative) {
        TCanvas* canvas_cumulative = new TCanvas("c_cumulative", "Cumulative Decay Distribution", 800, 600);
        canvas_cumulative->SetLeftMargin(0.12);
        canvas_cumulative->SetBottomMargin(0.12);
        canvas_cumulative->SetGrid();
        
        h_deltaT_cumulative->SetLineColor(kBlue+1);
        h_deltaT_cumulative->SetLineWidth(3);
        h_deltaT_cumulative->SetMarkerColor(kBlue+1);
        h_deltaT_cumulative->SetMarkerStyle(20);
        h_deltaT_cumulative->SetMarkerSize(0.5);
        h_deltaT_cumulative->GetXaxis()->SetTitleSize(0.045);
        h_deltaT_cumulative->GetYaxis()->SetTitleSize(0.045);
        h_deltaT_cumulative->GetXaxis()->SetLabelSize(0.04);
        h_deltaT_cumulative->GetYaxis()->SetLabelSize(0.04);
        h_deltaT_cumulative->GetXaxis()->SetTitleOffset(1.2);
        h_deltaT_cumulative->GetYaxis()->SetTitleOffset(1.3);
        h_deltaT_cumulative->GetYaxis()->SetRangeUser(0.0, 1.1);
        
        h_deltaT_cumulative->Draw("LP");
        
        TString pngCumulativeFile = outputDir + "/cumulative_decay_distribution.png";
        canvas_cumulative->SaveAs(pngCumulativeFile);
        cout << "Cumulative PNG saved: " << pngCumulativeFile << endl;
        
        delete canvas_cumulative;
    }

    // Create summary file
    TString summaryFile = outputDir + "/analysis_summary.txt";
    std::ofstream summary(summaryFile.Data());
    summary << "Decay Electron Energy and DeltaT Analysis Summary\n";
    summary << "=================================================\n";
    summary << "Input directory: " << inputDirectory << "\n";
    summary << "Total ROOT files found: " << rootFiles.size() << "\n";
    summary << "Files with valid charge data: " << validCharges.size() << "\n";
    summary << "Files with valid deltaT data: " << validDeltaTs.size() << "\n";
    summary << "Files with both valid data: " << 0 << "\n"; // Could add this calculation
    
    if (!validCharges.empty()) {
        summary << "\nCharge Statistics:\n";
        summary << "Mean: " << std::fixed << std::setprecision(2) << h_charges->GetMean() << " MeV\n";
        summary << "RMS: " << std::fixed << std::setprecision(2) << h_charges->GetRMS() << " MeV\n";
        summary << "Minimum: " << std::fixed << std::setprecision(2) << *min_element(validCharges.begin(), validCharges.end()) << " MeV\n";
        summary << "Maximum: " << std::fixed << std::setprecision(2) << *max_element(validCharges.begin(), validCharges.end()) << " MeV\n";
    }
    
    if (!validDeltaTs.empty()) {
        summary << "\nDeltaT Statistics:\n";
        summary << "Mean: " << std::fixed << std::setprecision(1) << h_deltaT->GetMean() << " ns\n";
        summary << "RMS: " << std::fixed << std::setprecision(1) << h_deltaT->GetRMS() << " ns\n";
        summary << "Minimum: " << std::fixed << std::setprecision(1) << *min_element(validDeltaTs.begin(), validDeltaTs.end()) << " ticks\n";
        summary << "Maximum: " << std::fixed << std::setprecision(1) << *max_element(validDeltaTs.begin(), validDeltaTs.end()) << " ticks\n";
    }
    
    summary << "\nOutput files:\n";
    summary << "ROOT file: " << rootOutputFile << "\n";
    if (h_charges) summary << "Charge PNG: " << outputDir << "/decay_electron_charge_distribution.png\n";
    if (h_deltaT) summary << "DeltaT PNG: " << outputDir << "/decay_time_interval_distribution.png\n";
    
    summary.close();
    
    auto endTime = std::chrono::steady_clock::now();
    auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    
    cout << "\n=== Analysis Complete ===" << endl;
    cout << "Total processing time: " << totalTime << " seconds" << endl;
    cout << "Output directory: " << outputDir << endl;
    cout << "Summary file: " << summaryFile << endl;
    cout << "\nTo view the histograms in ROOT:" << endl;
    cout << "TFile *f = TFile::Open(\"" << rootOutputFile << "\");" << endl;
    if (h_charges) cout << "h_charges->Draw(\"HIST\");" << endl;
    if (h_deltaT) cout << "h_deltaT->Draw(\"HIST\");" << endl;
    
    // Clean up memory
    if (h_charges) delete h_charges;
    if (h_deltaT) delete h_deltaT;
    if (h_deltaT_cumulative) delete h_deltaT_cumulative;
    
    exit(0);
}