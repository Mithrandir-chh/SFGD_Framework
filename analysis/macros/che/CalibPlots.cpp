// Gain Calibration Visualizer - Creates 2.5D gain maps and calibration files
#define THIS_NAME GainCalibrationVisualizer
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TPaveText.h"
#include "TText.h"

#include "../../src/tools/global_header.hh"

using namespace std;

// Structure to hold calibration data
struct CalibrationData {
    string view;
    int coord1, coord2;
    double gain;
    bool hasData;
    
    CalibrationData() : view(""), coord1(0), coord2(0), gain(1.0), hasData(false) {}
};

// Function to create output directory
string createOutputDirectory(const string& inputFile) {
    string inputDir = gSystem->DirName(inputFile.c_str());
    string outputDir = inputDir + "/GainCalibration";
    
    if (gSystem->AccessPathName(outputDir.c_str())) {
        cout << "Creating output directory: " << outputDir << endl;
        if (gSystem->mkdir(outputDir.c_str(), kTRUE) != 0) {
            cerr << "Warning: Could not create directory " << outputDir << endl;
            return "./GainCalibration";
        }
    }
    
    return outputDir;
}

void GainCalibrationVisualizer() {
    gROOT->SetBatch(kFALSE);  // Enable graphics for display
    
    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <ChannelGainAnalysis.root>" << endl;
        cerr << "Creates 2.5D gain calibration maps and calibration file" << endl;
        return;
    }
    
    string inputFile = argv[1];
    cout << "Input file: " << inputFile << endl;
    
    // Open input file
    TFile* fInput = TFile::Open(inputFile.c_str(), "READ");
    if (!fInput || fInput->IsZombie()) {
        cerr << "Error: Cannot open input file: " << inputFile << endl;
        return;
    }
    
    // Get the ChannelGainData tree
    TTree* tree = (TTree*)fInput->Get("ChannelGainData");
    if (!tree) {
        cerr << "Error: Cannot find ChannelGainData tree in input file" << endl;
        fInput->Close();
        return;
    }
    
    // Set up branch reading
    Char_t view_name[10];
    Int_t coord1, coord2;
    Double_t mean_gain, fit_mean;
    Bool_t fit_success;
    
    tree->SetBranchAddress("view", view_name);
    tree->SetBranchAddress("coord1", &coord1);
    tree->SetBranchAddress("coord2", &coord2);
    tree->SetBranchAddress("mean_gain", &mean_gain);
    tree->SetBranchAddress("fit_mean", &fit_mean);
    tree->SetBranchAddress("fit_success", &fit_success);
    
    Long64_t nEntries = tree->GetEntries();
    cout << "Number of channels in tree: " << nEntries << endl;
    
    // Store calibration data
    vector<CalibrationData> calibrationData;
    map<string, map<pair<int,int>, double>> viewGainMaps;
    
    // Read all channel data
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);
        
        string view = string(view_name);
        double gainValue = fit_success ? fit_mean : mean_gain;
        
        // Store for calibration file
        CalibrationData cal;
        cal.view = view;
        cal.coord1 = coord1;
        cal.coord2 = coord2;
        cal.gain = gainValue;
        cal.hasData = true;
        calibrationData.push_back(cal);
        
        // Store for histogram mapping
        viewGainMaps[view][make_pair(coord1, coord2)] = gainValue;
        
        if (entry % 1000 == 0) {
            cout << "Processing entry " << entry << "/" << nEntries << "\r" << flush;
        }
    }
    
    cout << "\nData loaded successfully!" << endl;
    
    // Create output directory and files
    string outputDir = createOutputDirectory(inputFile);
    string outputRoot = outputDir + "/GainCalibrationMaps.root";
    string outputCalib = outputDir + "/gain_calibration.txt";
    
    // Create output ROOT file
    TFile* fOutput = new TFile(outputRoot.c_str(), "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        cerr << "Error: Cannot create output file: " << outputRoot << endl;
        fInput->Close();
        return;
    }
    
    // Set up style
    gStyle->SetPalette(kViridis);
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.15);
    
    // Create 2.5D histograms for each view
    map<string, TH2F*> gainHistograms;
    
    // XY view: coord1=x (0-23), coord2=y (0-7)
    gainHistograms["XY"] = new TH2F("hGain_XY", "Gain Map - XY View;X Channel;Y Channel;Gain [PE/(2*TrackLength)]", 
                                    24, 0, 24, 8, 0, 8);
    
    // XZ view: coord1=x (0-23), coord2=z (0-47)
    gainHistograms["XZ"] = new TH2F("hGain_XZ", "Gain Map - XZ View;X Channel;Z Channel;Gain [PE/(2*TrackLength)]", 
                                    24, 0, 24, 48, 0, 48);
    
    // ZY view: coord1=z (0-47), coord2=y (0-7)
    gainHistograms["ZY"] = new TH2F("hGain_ZY", "Gain Map - ZY View;Z Channel;Y Channel;Gain [PE/(2*TrackLength)]", 
                                    48, 0, 48, 8, 0, 8);
    
    // Initialize all bins to -1.0 (dead channels)
    for (auto& histPair : gainHistograms) {
        TH2F* hist = histPair.second;
        for (int i = 1; i <= hist->GetNbinsX(); i++) {
            for (int j = 1; j <= hist->GetNbinsY(); j++) {
                hist->SetBinContent(i, j, -1.0);
                hist->SetMinimum(0.0);
                hist->SetMaximum(50.0);
            }
        }
    }
    
    // Fill histograms with actual gain data
    int channelsWithData = 0;
    for (const string& view : {"XY", "XZ", "ZY"}) {
        TH2F* hist = gainHistograms[view];
        
        for (const auto& channelPair : viewGainMaps[view]) {
            int c1 = channelPair.first.first;
            int c2 = channelPair.first.second;
            double gain = channelPair.second;
            
            // ROOT histogram bins are 1-indexed
            hist->SetBinContent(c1 + 1, c2 + 1, gain);
            channelsWithData++;
        }
        
        cout << "Filled " << viewGainMaps[view].size() << " channels for " << view << " view" << endl;
    }
    
    // Create canvas with all three views
    TCanvas* canvas = new TCanvas("canvas", "Gain Calibration Maps", 1800, 600);
    canvas->Divide(3, 1);
    
    // Draw XY view
    canvas->cd(1);
    gPad->SetRightMargin(0.15);
    gainHistograms["XY"]->Draw("COLZ");
    gainHistograms["XY"]->GetZaxis()->SetTitle("Gain");
    
    // Draw XZ view
    canvas->cd(2);
    gPad->SetRightMargin(0.15);
    gainHistograms["XZ"]->Draw("COLZ");
    gainHistograms["XZ"]->GetZaxis()->SetTitle("Gain");
    
    // Draw ZY view
    canvas->cd(3);
    gPad->SetRightMargin(0.15);
    gainHistograms["ZY"]->Draw("COLZ");
    gainHistograms["ZY"]->GetZaxis()->SetTitle("Gain");
    
    canvas->Update();
    
    // Save histograms to ROOT file
    fOutput->cd();
    for (auto& histPair : gainHistograms) {
        histPair.second->Write();
    }
    canvas->Write("GainMapsCanvas");
    
    // Save canvas as PNG
    string pngFile = outputDir + "/GainCalibrationMaps.png";
    canvas->SaveAs(pngFile.c_str());
    
    // Create calibration text file
    ofstream calibFile(outputCalib);
    if (!calibFile.is_open()) {
        cerr << "Error: Cannot create calibration file: " << outputCalib << endl;
    } else {
        calibFile << "# Gain Calibration File" << endl;
        calibFile << "# Generated from: " << inputFile << endl;
        calibFile << "# Format: view coord1 coord2 gain" << endl;
        calibFile << "# Views: XY (coord1=x, coord2=y), XZ (coord1=x, coord2=z), ZY (coord1=z, coord2=y)" << endl;
        calibFile << "# Gain = -1.0 stands for dead channels" << endl;
        calibFile << "" << endl;
        
        // Write all possible channels for each view
        int totalChannels = 0;
        int deadChannels = 0;
        
        // XY view: all combinations of x (0-23) and y (0-7)
        for (int x = 0; x < 24; x++) {
            for (int y = 0; y < 8; y++) {
                double gain = -1.0;  // Default for dead channels
                if (viewGainMaps["XY"].count(make_pair(x, y))) {
                    gain = viewGainMaps["XY"][make_pair(x, y)];
                } else {
                    deadChannels++;
                }
                calibFile << fixed << setprecision(6);
                calibFile << "XY " << x << " " << y << " " << gain << endl;
                totalChannels++;
            }
        }
        
        // XZ view: all combinations of x (0-23) and z (0-47)
        for (int x = 0; x < 24; x++) {
            for (int z = 0; z < 48; z++) {
                double gain = -1.0;  // Default for dead channels
                if (viewGainMaps["XZ"].count(make_pair(x, z))) {
                    gain = viewGainMaps["XZ"][make_pair(x, z)];
                } else {
                    deadChannels++;
                }
                calibFile << fixed << setprecision(6);
                calibFile << "XZ " << x << " " << z << " " << gain << endl;
                totalChannels++;
            }
        }
        
        // ZY view: all combinations of z (0-47) and y (0-7)
        for (int z = 0; z < 48; z++) {
            for (int y = 0; y < 8; y++) {
                double gain = -1.0;  // Default for dead channels
                if (viewGainMaps["ZY"].count(make_pair(z, y))) {
                    gain = viewGainMaps["ZY"][make_pair(z, y)];
                } else {
                    deadChannels++;
                }
                calibFile << fixed << setprecision(6);
                calibFile << "ZY " << z << " " << y << " " << gain << endl;
                totalChannels++;
            }
        }
        
        calibFile.close();
        
        cout << "\nCalibration file created: " << outputCalib << endl;
        cout << "Total channels: " << totalChannels << endl;
        cout << "Channels with gain data: " << channelsWithData << endl;
        cout << "Dead channels (gain = -1.0): " << deadChannels << endl;
    }
    
    // Print summary
    cout << "\n=== GAIN CALIBRATION SUMMARY ===" << endl;
    cout << "Input file: " << inputFile << endl;
    cout << "Output files created:" << endl;
    cout << "  ROOT file: " << outputRoot << endl;
    cout << "  PNG file: " << pngFile << endl;
    cout << "  Calibration file: " << outputCalib << endl;
    cout << "\nChannel statistics:" << endl;
    cout << "  XY channels with data: " << viewGainMaps["XY"].size() << " / " << (24*8) << endl;
    cout << "  XZ channels with data: " << viewGainMaps["XZ"].size() << " / " << (24*48) << endl;
    cout << "  ZY channels with data: " << viewGainMaps["ZY"].size() << " / " << (48*8) << endl;
    
    // Calculate gain statistics
    vector<double> allGains;
    for (const auto& cal : calibrationData) {
        if (cal.hasData && cal.gain > 0 && cal.gain < 200) {
            allGains.push_back(cal.gain);
        }
    }
    
    if (!allGains.empty()) {
        sort(allGains.begin(), allGains.end());
        double minGain = allGains.front();
        double maxGain = allGains.back();
        double medianGain = allGains[allGains.size()/2];
        double meanGain = 0;
        for (double g : allGains) meanGain += g;
        meanGain /= allGains.size();
        
        cout << "\nGain statistics:" << endl;
        cout << "  Range: [" << fixed << setprecision(3) << minGain << ", " << maxGain << "]" << endl;
        cout << "  Mean: " << meanGain << endl;
        cout << "  Median: " << medianGain << endl;
    }
    
    cout << "\nTo use the calibration file:" << endl;
    cout << "  - Read view, coord1, coord2, gain from each line" << endl;
    cout << "  - Apply: corrected_charge = raw_charge / gain" << endl;
    cout << "  - Channels with gain = -1.0 are dead (no correction applied)" << endl;
    
    // Cleanup
    fOutput->Close();
    fInput->Close();

    exit(0);
}