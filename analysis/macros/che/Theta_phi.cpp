// Analyzes trajectory data from multiple ROOT files and creates theta/phi distributions
#define THIS_NAME TrajectoryAnalysis
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPad.h>
#include <TGaxis.h>

// Include the framework headers
#include "../../src/tools/global_header.hh"

using namespace std;

class TrajectoryAnalyzer {
private:
    std::vector<std::string> rootFiles;
    TH1F* hTheta;
    TH1F* hPhi;
    TH2F* hThetaPhi;
    
    // Variables to read from tree
    Double_t theta, phi;
    
    // Data range tracking
    std::vector<Double_t> allTheta, allPhi;
    bool rangesSet;
    
public:
    TrajectoryAnalyzer() : rangesSet(false) {
        // Initialize with nullptr - will create after determining ranges
        hTheta = nullptr;
        hPhi = nullptr;
        hThetaPhi = nullptr;
    }
    
    ~TrajectoryAnalyzer() {
        if (hTheta) delete hTheta;
        if (hPhi) delete hPhi;
        if (hThetaPhi) delete hThetaPhi;
    }
    
    // Method to recursively find ROOT files by name pattern
    void findRootFiles(const std::string& directory, const std::string& namePattern = "") {
        try {
            for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
                if (entry.is_regular_file() && entry.path().extension() == ".root") {
                    std::string filename = entry.path().filename().string();
                    
                    // If no pattern specified, include all .root files
                    if (namePattern.empty()) {
                        rootFiles.push_back(entry.path().string());
                        std::cout << "Found ROOT file: " << entry.path().string() << std::endl;
                    }
                    // If pattern specified, check if filename contains the pattern
                    else if (filename.find(namePattern) != std::string::npos) {
                        rootFiles.push_back(entry.path().string());
                        std::cout << "Found matching ROOT file: " << entry.path().string() << std::endl;
                    }
                }
            }
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error accessing directory: " << e.what() << std::endl;
        }
        
        if (rootFiles.empty()) {
            if (namePattern.empty()) {
                std::cout << "No ROOT files found in directory: " << directory << std::endl;
            } else {
                std::cout << "No ROOT files matching pattern '" << namePattern 
                         << "' found in directory: " << directory << std::endl;
            }
        }
    }
    
    // Method to determine optimal bin ranges and create histograms
    void createAdaptiveHistograms() {
        if (allTheta.empty() || allPhi.empty()) {
            std::cerr << "No data collected for range determination!" << std::endl;
            return;
        }
        
        // Find min/max values
        auto thetaMinMax = std::minmax_element(allTheta.begin(), allTheta.end());
        auto phiMinMax = std::minmax_element(allPhi.begin(), allPhi.end());
        
        double thetaMin = *thetaMinMax.first;
        double thetaMax = *thetaMinMax.second;
        double phiMin = *phiMinMax.first;
        double phiMax = *phiMinMax.second;
        
        // Add 5% padding to ranges for better visualization
        double thetaRange = thetaMax - thetaMin;
        double phiRange = phiMax - phiMin;
        
        double thetaPadding = thetaRange * 0.05;
        double phiPadding = phiRange * 0.05;
        
        thetaMin -= thetaPadding;
        thetaMax += thetaPadding;
        phiMin -= phiPadding;
        phiMax += phiPadding;
        
        // For theta (typically 0 to 180 degrees), ensure minimum doesn't go below 0
        if (thetaMin < 0) {
            thetaMin = 0;
        }
        
        // For phi (typically -180 to 180 or 0 to 360), you might want similar constraints
        // Uncomment and adjust based on your phi convention:
        // if (phiMin < -180) phiMin = -180;
        // if (phiMax > 180) phiMax = 180;
        // OR for 0-360 convention:
        // if (phiMin < 0) phiMin = 0;
        // if (phiMax > 360) phiMax = 360;
        
        // Determine number of bins based on data size and range
        int nData = allTheta.size();
        int nBinsTheta = std::min(100, std::max(20, (int)sqrt(nData)));
        int nBinsPhi = std::min(100, std::max(20, (int)sqrt(nData)));
        
        // Calculate bin width
        double thetaBinWidth = (thetaMax - thetaMin) / nBinsTheta;
        double phiBinWidth = (phiMax - phiMin) / nBinsPhi;
        
        // Round to nice numbers for better bin edges, but respect physical constraints
        double thetaMinRounded = floor(thetaMin / thetaBinWidth) * thetaBinWidth;
        double thetaMaxRounded = ceil(thetaMax / thetaBinWidth) * thetaBinWidth;
        double phiMinRounded = floor(phiMin / phiBinWidth) * phiBinWidth;
        double phiMaxRounded = ceil(phiMax / phiBinWidth) * phiBinWidth;
        
        // Ensure theta doesn't go below 0 after rounding
        if (thetaMinRounded < 0) {
            thetaMinRounded = 0;
        }
        
        // Recalculate number of bins to maintain approximately the same bin width
        nBinsTheta = (int)round((thetaMaxRounded - thetaMinRounded) / thetaBinWidth);
        nBinsPhi = (int)round((phiMaxRounded - phiMinRounded) / phiBinWidth);
        
        // Ensure we have at least some minimum number of bins
        if (nBinsTheta < 10) nBinsTheta = 10;
        if (nBinsPhi < 10) nBinsPhi = 10;
        
        std::cout << "\n=== ADAPTIVE BINNING INFORMATION ===" << std::endl;
        std::cout << "Data points collected: " << nData << std::endl;
        std::cout << "Theta range: " << *thetaMinMax.first << " to " << *thetaMinMax.second << " deg" << std::endl;
        std::cout << "  -> Histogram range (aligned): " << thetaMinRounded << " to " << thetaMaxRounded << " deg" << std::endl;
        std::cout << "  -> Bins: " << nBinsTheta << " (width: " << (thetaMaxRounded - thetaMinRounded) / nBinsTheta << " deg)" << std::endl;
        std::cout << "Phi range: " << *phiMinMax.first << " to " << *phiMinMax.second << " deg" << std::endl;
        std::cout << "  -> Histogram range (aligned): " << phiMinRounded << " to " << phiMaxRounded << " deg" << std::endl;
        std::cout << "  -> Bins: " << nBinsPhi << " (width: " << (phiMaxRounded - phiMinRounded) / nBinsPhi << " deg)" << std::endl;
        std::cout << "=====================================" << std::endl;
        
        // Create histograms with adaptive ranges
        hTheta = new TH1F("hTheta", "Theta Distribution;Theta [deg];Counts", 
                         nBinsTheta, thetaMinRounded, thetaMaxRounded);
        hPhi = new TH1F("hPhi", "Phi Distribution;Phi [deg];Counts", 
                       nBinsPhi, phiMinRounded, phiMaxRounded);
        hThetaPhi = new TH2F("hThetaPhi", "Theta vs Phi Distribution;Phi [deg];Theta [deg];Counts", 
                           nBinsPhi, phiMinRounded, phiMaxRounded, 
                           nBinsTheta, thetaMinRounded, thetaMaxRounded);
        
        // Enable automatic error calculation (Sumw2)
        hTheta->Sumw2();
        hPhi->Sumw2();
        hThetaPhi->Sumw2();
        
        rangesSet = true;
    }
    
    // Method to process a single ROOT file
    bool processFile(const std::string& filename, bool collectOnly = false) {
        TFile* file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return false;
        }
        
        TTree* tree = dynamic_cast<TTree*>(file->Get("Trajectory"));
        if (!tree) {
            std::cerr << "Tree 'Trajectory' not found in file: " << filename << std::endl;
            file->Close();
            return false;
        }
        
        // Set branch addresses
        tree->SetBranchAddress("theta", &theta);
        tree->SetBranchAddress("phi", &phi);
        
        Long64_t nEntries = tree->GetEntries();
        
        if (collectOnly) {
            std::cout << "Collecting ranges from " << nEntries << " entries in " << filename << std::endl;
            // Collect all values for range determination
            for (Long64_t i = 0; i < nEntries; i++) {
                tree->GetEntry(i);
                allTheta.push_back(theta);
                allPhi.push_back(phi);
            }
        } else {
            std::cout << "Processing " << nEntries << " entries from " << filename << std::endl;
            
            if (!rangesSet) {
                std::cerr << "Error: Histograms not created yet!" << std::endl;
                file->Close();
                return false;
            }
            
            // Fill histograms
            for (Long64_t i = 0; i < nEntries; i++) {
                tree->GetEntry(i);
                
                // Fill histograms with weight 1.0 (automatic error calculation)
                hTheta->Fill(theta);
                hPhi->Fill(phi);
                hThetaPhi->Fill(phi, theta);
                
                // Progress indicator for large files
                if (i % 10000 == 0 && i > 0) {
                    std::cout << "Processed " << i << "/" << nEntries << " entries\r" << std::flush;
                }
            }
            
            // Add 10% systematic errors
            addSystematicErrors();
            std::cout << std::endl;
        }
        
        file->Close();
        return true;
    }
    
    // Method to process all found ROOT files
    void processAllFiles() {
        if (rootFiles.empty()) {
            std::cerr << "No ROOT files found to process!" << std::endl;
            return;
        }
        
        std::cout << "\n=== PHASE 1: COLLECTING DATA RANGES ===" << std::endl;
        // First pass: collect all data to determine ranges
        for (const auto& filename : rootFiles) {
            processFile(filename, true);  // collectOnly = true
        }
        
        // Create histograms with adaptive ranges
        createAdaptiveHistograms();
        
        std::cout << "\n=== PHASE 2: FILLING HISTOGRAMS ===" << std::endl;
        // Second pass: fill histograms
        int successCount = 0;
        for (const auto& filename : rootFiles) {
            if (processFile(filename, false)) {  // collectOnly = false
                successCount++;
            }
        }
        
        std::cout << "\nSuccessfully processed " << successCount << "/" << rootFiles.size() << " files" << std::endl;
    }
    
    // Method to set error bars to exactly 10% of bin content
    void addSystematicErrors() {
        if (!hTheta || !hPhi || !hThetaPhi) return;
        
        // Set error to exactly 10% of bin content for theta histogram
        for (int i = 1; i <= hTheta->GetNbinsX(); i++) {
            double content = hTheta->GetBinContent(i);
            double error = content * 0.10;  // 10% error
            hTheta->SetBinError(i, error);
        }
        
        // Set error to exactly 10% of bin content for phi histogram
        for (int i = 1; i <= hPhi->GetNbinsX(); i++) {
            double content = hPhi->GetBinContent(i);
            double error = content * 0.10;  // 10% error
            hPhi->SetBinError(i, error);
        }
        
        // Set error to exactly 10% of bin content for 2D histogram
        for (int i = 1; i <= hThetaPhi->GetNbinsX(); i++) {
            for (int j = 1; j <= hThetaPhi->GetNbinsY(); j++) {
                double content = hThetaPhi->GetBinContent(i, j);
                double error = content * 0.10;  // 10% error
                hThetaPhi->SetBinError(i, j, error);
            }
        }
        
        std::cout << "Set error bars to 10% of bin content for all histograms" << std::endl;
    }
    
    // Method to save histograms to a new ROOT file
    void saveToRootFile(const std::string& outputFile = "trajectory_analysis.root") {
        if (!rangesSet || !hTheta || !hPhi || !hThetaPhi) {
            std::cerr << "Error: Histograms not created yet!" << std::endl;
            return;
        }
        
        TFile* outFile = TFile::Open(outputFile.c_str(), "RECREATE");
        if (!outFile || outFile->IsZombie()) {
            std::cerr << "Error creating output file: " << outputFile << std::endl;
            return;
        }
        
        // Write histograms to file
        hTheta->Write();
        hPhi->Write();
        hThetaPhi->Write();
        
        outFile->Close();
        std::cout << "Histograms saved to " << outputFile << std::endl;
    }
    
    // Method to create and display plots
    void createPlots(const std::string& outputDir = "trajectory_analysis_output") {
        if (!rangesSet || !hTheta || !hPhi || !hThetaPhi) {
            std::cerr << "Error: Histograms not created yet!" << std::endl;
            return;
        }
        
        // Create output directory if it doesn't exist
        gSystem->mkdir(outputDir.c_str(), kTRUE);
        
        // Set ROOT style
        gStyle->SetOptStat(1111);
        gStyle->SetPadGridX(kTRUE);
        gStyle->SetPadGridY(kTRUE);
        
        // Create canvas for theta distribution with solid angle axis
        TCanvas* cTheta = new TCanvas("cTheta", "Theta Distribution", 800, 600);
        
        // Create the solid angle histogram
        TH1F* hOneMinusCosTheta = new TH1F("hOneMinusCosTheta", "", 
                                            hTheta->GetNbinsX(), 
                                            hTheta->GetXaxis()->GetXmin(), 
                                            hTheta->GetXaxis()->GetXmax());
        
        // Find the maximum bin content for normalization
        double maxCounts = hTheta->GetMaximum();
        
        // Fill the solid angle histogram with normalized values
        for (int i = 1; i <= hTheta->GetNbinsX(); i++) {
            double thetaDeg = hTheta->GetBinCenter(i);
            double thetaRad = thetaDeg * TMath::Pi() / 180.0;  // Convert to radians
            double oneMinusCos = 1.0 - cos(thetaRad);
            
            // Normalize so that max of solid angle corresponds to maxCounts
            double normalizedValue = (oneMinusCos) * maxCounts;
            hOneMinusCosTheta->SetBinContent(i, normalizedValue);
        }
        
        // Set up the main histogram
        hTheta->SetLineColor(kBlue);
        hTheta->SetFillColor(kBlue);
        hTheta->SetFillStyle(3005);
        hTheta->GetYaxis()->SetTitle("Counts");
        hTheta->GetYaxis()->SetTitleOffset(1.3);
        
        // Draw the main histogram
        hTheta->Draw("HIST");
        hTheta->Draw("E1 SAME");
        
        // Set up the solid angle curve
        hOneMinusCosTheta->SetLineColor(kRed);
        hOneMinusCosTheta->SetLineWidth(2);
        hOneMinusCosTheta->SetLineStyle(2);  // Dashed line
        hOneMinusCosTheta->Draw("L SAME");
        
        // Create a transparent pad for the second y-axis
        cTheta->cd();
        TPad* overlay = new TPad("overlay", "", 0, 0, 1, 1);
        overlay->SetFillStyle(4000);  // Transparent
        overlay->SetFrameFillStyle(0);
        overlay->Draw();
        overlay->cd();
        
        // Create a histogram just for the second axis
        TH1F* hAxis = (TH1F*)hTheta->Clone("hAxis");
        hAxis->Reset();
        hAxis->GetYaxis()->SetRangeUser(0, maxCounts);
        hAxis->GetYaxis()->SetTitle("Solid Angle");
        hAxis->GetYaxis()->SetTitleOffset(1.3);
        hAxis->GetYaxis()->SetLabelColor(kRed);
        hAxis->GetYaxis()->SetTitleColor(kRed);
        hAxis->GetXaxis()->SetLabelOffset(999);
        hAxis->GetXaxis()->SetTickLength(0);
        
        // Draw the second axis
        hAxis->Draw("AXIG Y+");
        
        // Add custom labels for the right y-axis showing solid angle values
        TGaxis* rightAxis = new TGaxis(hTheta->GetXaxis()->GetXmax(), 0,
                                       hTheta->GetXaxis()->GetXmax(), maxCounts,
                                       0, 1, 510, "+L");
        rightAxis->SetTitle("Solid Angle");
        rightAxis->SetTitleOffset(1.1);
        rightAxis->SetLabelColor(kRed);
        rightAxis->SetTitleColor(kRed);
        rightAxis->SetLineColor(kRed);
        rightAxis->Draw();
        
        // Add a legend
        cTheta->cd();
        TLegend* legend = new TLegend(0.55, 0.75, 0.85, 0.88);
        legend->AddEntry(hTheta, "Data", "f");
        legend->AddEntry(hOneMinusCosTheta, "Solid Angle", "l");
        legend->SetBorderSize(1);
        legend->Draw();
        
        cTheta->SaveAs((outputDir + "/theta_distribution.png").c_str());
        cTheta->SaveAs((outputDir + "/theta_distribution.pdf").c_str());  // Also save as PDF for better quality
        
        // Clean up
        delete hOneMinusCosTheta;
        delete overlay;
        delete rightAxis;
        
        // Create canvas for phi distribution (unchanged)
        TCanvas* cPhi = new TCanvas("cPhi", "Phi Distribution", 800, 600);
        hPhi->SetLineColor(kRed);
        hPhi->SetFillColor(kRed);
        hPhi->SetFillStyle(3005);
        hPhi->Draw("HIST");
        hPhi->Draw("E1 SAME");
        cPhi->SaveAs((outputDir + "/phi_distribution.png").c_str());
        
        // Create canvas for 2D theta-phi distribution
        TCanvas* cThetaPhi = new TCanvas("cThetaPhi", "Theta vs Phi Distribution", 900, 700);
        hThetaPhi->SetOption("COLZ");
        hThetaPhi->Draw("COLZ");
        cThetaPhi->SaveAs((outputDir + "/theta_phi_2d.png").c_str());
        
        std::cout << "Plots saved to " << outputDir << " directory" << std::endl;
        std::cout << "  - theta_distribution.png (histogram with error bands and solid angle)" << std::endl;
        std::cout << "  - theta_distribution.pdf (higher quality version)" << std::endl;
        std::cout << "  - phi_distribution.png (histogram with error bands)" << std::endl;
        std::cout << "  - theta_phi_2d.png" << std::endl;
    }
    
    // Method to print statistics
    void printStatistics() {
        if (!rangesSet || !hTheta || !hPhi) {
            std::cerr << "Error: Histograms not created yet!" << std::endl;
            return;
        }
        
        std::cout << "\n=== TRAJECTORY ANALYSIS STATISTICS ===" << std::endl;
        std::cout << "Total entries processed: " << hTheta->GetEntries() << std::endl;
        std::cout << "\nTheta Distribution:" << std::endl;
        std::cout << "  Range: " << hTheta->GetXaxis()->GetXmin() << " to " << hTheta->GetXaxis()->GetXmax() << " deg" << std::endl;
        std::cout << "  Bins: " << hTheta->GetNbinsX() << std::endl;
        std::cout << "  Mean: " << hTheta->GetMean() << " ± " << hTheta->GetMeanError() << " deg" << std::endl;
        std::cout << "  RMS:  " << hTheta->GetRMS() << " ± " << hTheta->GetRMSError() << " deg" << std::endl;
        std::cout << "\nPhi Distribution:" << std::endl;
        std::cout << "  Range: " << hPhi->GetXaxis()->GetXmin() << " to " << hPhi->GetXaxis()->GetXmax() << " deg" << std::endl;
        std::cout << "  Bins: " << hPhi->GetNbinsX() << std::endl;
        std::cout << "  Mean: " << hPhi->GetMean() << " ± " << hPhi->GetMeanError() << " deg" << std::endl;
        std::cout << "  RMS:  " << hPhi->GetRMS() << " ± " << hPhi->GetRMSError() << " deg" << std::endl;
        std::cout << "=======================================" << std::endl;
    }
};

// Function to create output directory if it doesn't exist
TString createOutputDirectory(const TString& inputDir) {
    // Create the subfolder path
    TString outputDir = inputDir + "/trajectory_analysis_output";
    
    // Check if directory exists, create if not
    if (gSystem->AccessPathName(outputDir)) {
        std::cout << "Creating output directory: " << outputDir << std::endl;
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            std::cerr << "Using current directory instead." << std::endl;
            return "./trajectory_analysis_output";
        }
    } else {
        std::cout << "Using existing output directory: " << outputDir << std::endl;
    }
    
    return outputDir;
}

// Main entry point function that matches THIS_NAME
void TrajectoryAnalysis() {
    // Get command line arguments from the global application
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Check arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_directory> [name_pattern] [output_filename.root]" << std::endl;
        std::cerr << "Examples:" << std::endl;
        std::cerr << "  " << argv[0] << " /path/to/root/files" << std::endl;
        std::cerr << "  " << argv[0] << " /path/to/root/files fitted" << std::endl;
        std::cerr << "  " << argv[0] << " /path/to/root/files trajectory results.root" << std::endl;
        std::cerr << "  " << argv[0] << " /path/to/root/files \"run_2024\" analysis.root" << std::endl;
        return;
    }
    
    TString inputDir = argv[1];
    TString namePattern = (argc > 2) ? argv[2] : "";
    TString outputFilename = (argc > 3) ? argv[3] : "trajectory_analysis.root";
    
    // If third argument looks like a .root file and no fourth argument, treat it as output filename
    if (argc == 3 && namePattern.EndsWith(".root")) {
        outputFilename = namePattern;
        namePattern = "";
    }
    
    std::cout << "Input directory: " << inputDir << std::endl;
    if (!namePattern.IsNull()) {
        std::cout << "Name pattern: " << namePattern << std::endl;
    }
    std::cout << "Output filename: " << outputFilename << std::endl;
    
    // Check if input directory exists
    if (gSystem->AccessPathName(inputDir)) {
        std::cerr << "Error: Input directory " << inputDir << " does not exist" << std::endl;
        return;
    }
    
    // Create analyzer instance
    TrajectoryAnalyzer analyzer;
    
    // Find ROOT files with optional name pattern
    if (namePattern.IsNull()) {
        std::cout << "Searching for all ROOT files in: " << inputDir << std::endl;
        analyzer.findRootFiles(inputDir.Data());
    } else {
        std::cout << "Searching for ROOT files containing '" << namePattern << "' in: " << inputDir << std::endl;
        analyzer.findRootFiles(inputDir.Data(), namePattern.Data());
    }
    
    // Process all files
    analyzer.processAllFiles();
    
    // Print statistics
    analyzer.printStatistics();
    
    // Create output directory and full output path
    TString outputDir = createOutputDirectory(inputDir);
    TString fullOutputPath = outputDir + "/" + outputFilename;
    
    // Save results and create plots
    TString plotsDir = outputDir + "/plots";
    analyzer.saveToRootFile(fullOutputPath.Data());
    analyzer.createPlots(plotsDir.Data());
    
    std::cout << "\nTrajectory analysis complete!" << std::endl;
    std::cout << "Results saved to: " << fullOutputPath << std::endl;
    std::cout << "Plots saved to: " << plotsDir << std::endl;
    std::cout << "\nTo view the results in ROOT:" << std::endl;
    std::cout << "  root " << fullOutputPath << std::endl;
    std::cout << "  hTheta->Draw()" << std::endl;
    std::cout << "  // or" << std::endl;
    std::cout << "  hPhi->Draw()" << std::endl;
    std::cout << "  // or for 2D plot:" << std::endl;
    std::cout << "  hThetaPhi->Draw(\"COLZ\")" << std::endl;
    exit(0);
}