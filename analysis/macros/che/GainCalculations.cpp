// Channel Gain Analysis - Analyzes gain per channel from trajectory voxel data
#define THIS_NAME ChannelGainAnalysis
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TApplication.h"

#include "../../src/tools/global_header.hh"

using namespace std;

// Structure to hold channel identification
struct ChannelID {
    string view;  // "XY", "XZ", "ZY"
    int coord1;   // X for XY/XZ, Z for ZY
    int coord2;   // Y for XY/ZY, Z for XZ
    
    bool operator<(const ChannelID& other) const {
        if (view != other.view) return view < other.view;
        if (coord1 != other.coord1) return coord1 < other.coord1;
        return coord2 < other.coord2;
    }
    
    string getName() const {
        return view + "_" + to_string(coord1) + "_" + to_string(coord2);
    }
};

// Structure to hold voxel data for analysis
struct VoxelData {
    int x, y, z;
    double pathLength;
    double chargeXY, chargeXZ, chargeZY;
    bool onTrajectory;
    bool hasHit;
    bool usedInAnalysis;
};

// Structure to hold channel data for one entry
struct ChannelData {
    double totalCharge;
    double totalPathLength;
    int nVoxels;
    
    ChannelData() : totalCharge(0), totalPathLength(0), nVoxels(0) {}
};

// Structure to hold precomputed distance calculation info
struct ChannelDistanceInfo {
    ChannelID channelID;
    std::function<double(const VoxelData&)> distanceCalculator;
    string ruleDescription; // for debugging/validation
};

// Function to generate all possible channels and their distance calculators
map<ChannelID, std::function<double(const VoxelData&)>> generateChannelDistanceMap() {
    map<ChannelID, std::function<double(const VoxelData&)>> distanceMap;
    
    // Generate all possible channels based on your detector geometry
    // You'll need to know the coordinate ranges for each view
    
    // XY channels
    for (int x = 0; x < 24; x++) {        // Replace MAX_X with your actual range
        for (int y = 0; y < 8; y++) {    // Replace MAX_Y with your actual range
            ChannelID channel;
            channel.view = "XY";
            channel.coord1 = x;
            channel.coord2 = y;
            
            if (x % 2 == 1) {  // Odd X
                distanceMap[channel] = [](const VoxelData& v) -> double {
                    return abs(v.z);
                };
            } else {  // Even X
                distanceMap[channel] = [](const VoxelData& v) -> double {
                    return abs(48 - v.z - 1);
                };
            }
        }
    }
    
    // ZY channels
    for (int z = 0; z < 48; z++) {
        for (int y = 0; y < 8; y++) {
            ChannelID channel;
            channel.view = "ZY";
            channel.coord1 = z;
            channel.coord2 = y;
            
            if (z % 2 == 0) {
                distanceMap[channel] = [](const VoxelData& v) -> double {
                    return abs(24 - v.x - 1);
                };
            } 

            else {
                distanceMap[channel] = [](const VoxelData& v) -> double {
                    return abs(v.x);
                };
            }

        }
    }
    
    for (int x = 0; x < 24; x++) {
        for (int z = 0; z < 48; z++) {
            ChannelID channel;
            channel.view = "XZ";
            channel.coord1 = x;
            channel.coord2 = z;

            if (z % 2 == 0) {
                distanceMap[channel] = [](const VoxelData& v) -> double {
                        return abs(8 - v.y - 1); 
                };
            }

            else {
                distanceMap[channel] = [](const VoxelData& v) -> double {
                    return abs(v.y);
                };
            }
            
        }
    }
    return distanceMap;
}

// Structure to handle calibration data loading and application
struct CalibrationData {
    map<string, map<pair<int,int>, double>> gainMap;
    
    bool LoadCalibration(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cout << "Error: Cannot open calibration file " << filename << endl;
            return false;
        }
        
        string view;
        int coord1, coord2;
        double gain;
        
        while (file >> view >> coord1 >> coord2 >> gain) {
            // Skip invalid gains (marked as -1.000000)
            if (gain > 0) {
                gainMap[view][make_pair(coord1, coord2)] = gain;
            }
        }
        file.close();
        
        cout << "Loaded calibration for " << gainMap.size() << " views" << endl;
        
        // Print summary of loaded calibration data
        for (const auto& viewPair : gainMap) {
            cout << "  " << viewPair.first << ": " << viewPair.second.size() << " channels" << endl;
        }
        
        return true;
    }
    
    double GetGain(const string& view, int coord1, int coord2) {
        auto viewIt = gainMap.find(view);
        if (viewIt != gainMap.end()) {
            auto coordIt = viewIt->second.find(make_pair(coord1, coord2));
            if (coordIt != viewIt->second.end()) {
                return coordIt->second;
            }
        }
        return -1.0; // Return -1 if calibration not found
    }
};

// Global calibration data object
CalibrationData g_calibData;

// Global flag to control uncalibration
bool g_useUncalibration = false;

// Global distance map (initialize once at program start)
map<ChannelID, std::function<double(const VoxelData&)>> g_channelDistanceMap;

void initializeChannelDistanceMap() {
    cout << "Generating channel distance lookup table..." << endl;
    g_channelDistanceMap = generateChannelDistanceMap();
    cout << "Generated distance calculators for " << g_channelDistanceMap.size() << " channels" << endl;
}

// Simple lookup function
double calculateVoxelToChannelDistance(const VoxelData& voxel, const ChannelID& channel) {
    auto it = g_channelDistanceMap.find(channel);
    if (it != g_channelDistanceMap.end()) {
        return it->second(voxel);
    }
    
    cerr << "Warning: No distance calculator found for channel " << channel.getName() << endl;
    return 0.0;
}

double applyAttenuationCorrection(double rawGain, double distance, const string& view) {
    double alpha = 0.1399;
    double Ls = 63.06;
    double Ll = 4000;
    
    // Exponential attenuation correction
    double correctionFactor = (alpha * exp ((-distance)/Ls)
                                +(1.0 - alpha) * exp ((-distance)/Ll));
    
    return rawGain / correctionFactor;
}

// Function to get charge for a specific view
double getViewCharge(const VoxelData& voxel, const string& view) {
    if (view == "XY") return voxel.chargeXY;
    if (view == "XZ") return voxel.chargeXZ;
    if (view == "ZY") return voxel.chargeZY;
    return 0.0;
}

// Function to get uncalibrated (raw) charge by multiplying by gain
double getUncalibratedCharge(const VoxelData& voxel, const string& view) {
    double calibratedCharge = getViewCharge(voxel, view);
    if (calibratedCharge <= 0) return 0.0;
    
    // Determine coordinates based on view
    int coord1, coord2;
    if (view == "XY") {
        coord1 = voxel.x;
        coord2 = voxel.y;
    } else if (view == "XZ") {
        coord1 = voxel.x;
        coord2 = voxel.z;
    } else if (view == "ZY") {
        coord1 = voxel.z;
        coord2 = voxel.y;
    } else {
        cerr << "Warning: Unknown view " << view << endl;
        return calibratedCharge;
    }
    
    // Get gain from calibration data
    double gain = g_calibData.GetGain(view, coord1, coord2);
    if (gain > 0) {
        // Undo calibration: raw_charge = calibrated_charge * gain
        double rawCharge = calibratedCharge * gain;
        return rawCharge;
    } else {
        // No calibration available, assume input is already raw
        cout << "Warning: No calibration found for " << view << "(" << coord1 << "," << coord2 
             << "), using original charge" << endl;
        return calibratedCharge;
    }
}

// Function to create channel ID from voxel and view
ChannelID createChannelID(const VoxelData& voxel, const string& view) {
    ChannelID channel;
    channel.view = view;
    
    if (view == "XY") {
        channel.coord1 = voxel.x;
        channel.coord2 = voxel.y;
    } else if (view == "XZ") {
        channel.coord1 = voxel.x;
        channel.coord2 = voxel.z;
    } else if (view == "ZY") {
        channel.coord1 = voxel.z;
        channel.coord2 = voxel.y;
    }
    
    return channel;
}

// Function to create output directory
string createOutputDirectory(const string& inputFile) {
    string inputDir = gSystem->DirName(inputFile.c_str());
    string outputDir = inputDir + "/ChannelGainAnalysis";
    
    if (gSystem->AccessPathName(outputDir.c_str())) {
        cout << "Creating output directory: " << outputDir << endl;
        if (gSystem->mkdir(outputDir.c_str(), kTRUE) != 0) {
            cerr << "Warning: Could not create directory " << outputDir << endl;
            return "./ChannelGainAnalysis";
        }
    }
    
    return outputDir;
}

void ChannelGainAnalysis() {
    gROOT->SetBatch(kTRUE);
    
    initializeChannelDistanceMap();

    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    if (argc < 2 || argc > 3) {
        cerr << "Usage: " << argv[0] << " <Summary_AllFiles.root> [calibration_file.txt]" << endl;
        cerr << "Analyzes channel gain from consolidated trajectory voxel data" << endl;
        cerr << "Optional: calibration_file.txt to undo existing calibration in input data" << endl;
        cerr << "If no calibration file provided, assumes input data is already raw" << endl;
        return;
    }
    
    string inputFile = argv[1];
    string calibrationFile = "";
    
    if (argc == 3) {
        calibrationFile = argv[2];
        if (!calibrationFile.empty()) {
            g_useUncalibration = true;
        }
    }
    
    cout << "Input file: " << inputFile << endl;
    
    if (g_useUncalibration) {
        cout << "Calibration file: " << calibrationFile << endl;
        // Load calibration data to undo existing calibration
        cout << "Loading calibration data to undo existing calibration..." << endl;
        if (!g_calibData.LoadCalibration(calibrationFile)) {
            cerr << "Error: Failed to load calibration file: " << calibrationFile << endl;
            return;
        }
    } else {
        cout << "No calibration file provided - assuming input data is already raw" << endl;
    }

    // Open input file
    TFile* fInput = TFile::Open(inputFile.c_str(), "READ");
    if (!fInput || fInput->IsZombie()) {
        cerr << "Error: Cannot open input file: " << inputFile << endl;
        return;
    }
    
    // Get the ConsolidatedHits tree
    TTree* tree = (TTree*)fInput->Get("ConsolidatedHits");
    if (!tree) {
        cerr << "Error: Cannot find ConsolidatedHits tree in input file" << endl;
        fInput->Close();
        return;
    }
    
    // Set up branch reading
    vector<Int_t>* file_ids = nullptr;
    vector<Int_t>* voxel_xs = nullptr;
    vector<Int_t>* voxel_ys = nullptr;
    vector<Int_t>* voxel_zs = nullptr;
    vector<Double_t>* charges = nullptr;
    vector<Double_t>* confidences = nullptr;
    vector<Double_t>* chargeXYs = nullptr;
    vector<Double_t>* chargeXZs = nullptr;
    vector<Double_t>* chargeZYs = nullptr;
    vector<Double_t>* path_lengths = nullptr;
    vector<Bool_t>* has_hits = nullptr;
    vector<Bool_t>* on_trajectories = nullptr;
    vector<Bool_t>* used_in_analyses = nullptr;
    
    tree->SetBranchAddress("file_id", &file_ids);
    tree->SetBranchAddress("voxel_x", &voxel_xs);
    tree->SetBranchAddress("voxel_y", &voxel_ys);
    tree->SetBranchAddress("voxel_z", &voxel_zs);
    tree->SetBranchAddress("charge", &charges);
    tree->SetBranchAddress("confidence", &confidences);
    tree->SetBranchAddress("chargeXY", &chargeXYs);
    tree->SetBranchAddress("chargeXZ", &chargeXZs);
    tree->SetBranchAddress("chargeZY", &chargeZYs);
    tree->SetBranchAddress("path_length", &path_lengths);
    tree->SetBranchAddress("has_hit", &has_hits);
    tree->SetBranchAddress("on_trajectory", &on_trajectories);
    tree->SetBranchAddress("used_in_analysis", &used_in_analyses);
    
    Long64_t nEntries = tree->GetEntries();
    cout << "Number of entries (files): " << nEntries << endl;
    
    // Map to store all gain values per channel across all entries
    map<ChannelID, vector<double>> channelGains;
    
    // Map to store detailed gain data for logging
    map<ChannelID, vector<pair<int, double>>> channelGainDetails; // pair<entry_number, gain>
    
    // Process each entry
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);
        
        if (entry % 10 == 0 || entry == nEntries - 1) {
            cout << "Processing entry " << entry << "/" << nEntries 
                 << " (" << (100.0 * entry / nEntries) << "% done)\r" << flush;
        }
        
        // Convert vector data to VoxelData structures for this entry
        vector<VoxelData> voxels;
        size_t nVoxelsInEntry = voxel_xs->size();
        
        for (size_t i = 0; i < nVoxelsInEntry; i++) {
            VoxelData voxel;
            voxel.x = (*voxel_xs)[i];
            voxel.y = (*voxel_ys)[i];
            voxel.z = (*voxel_zs)[i];
            voxel.pathLength = (*path_lengths)[i];
            voxel.chargeXY = (*chargeXYs)[i];
            voxel.chargeXZ = (*chargeXZs)[i];
            voxel.chargeZY = (*chargeZYs)[i];
            voxel.onTrajectory = (*on_trajectories)[i];
            voxel.hasHit = (*has_hits)[i];
            voxel.usedInAnalysis = (*used_in_analyses)[i];
            
            voxels.push_back(voxel);
        }
        
        // Filter voxels: on_trajectory AND has_hit AND used_in_analysis
        vector<VoxelData> filteredVoxels;
        for (const auto& voxel : voxels) {
            if (voxel.onTrajectory && voxel.hasHit && voxel.usedInAnalysis && voxel.pathLength >= 0.7) {
                filteredVoxels.push_back(voxel);
            }
        }
        
        // Group voxels by channels for each view
        map<ChannelID, ChannelData> entryChannels;
        
        vector<string> views = {"XY", "XZ", "ZY"};
        
        for (const string& view : views) {
            // Group voxels by this view's channel
            map<ChannelID, vector<VoxelData>> viewChannels;
            
            for (const auto& voxel : filteredVoxels) {
                // Only include voxel if it has charge in this view
                double viewCharge = getViewCharge(voxel, view);
                if (viewCharge > 0) {
                    ChannelID channelID = createChannelID(voxel, view);
                    viewChannels[channelID].push_back(voxel);
                }
            }
            
            for (const auto& channelPair : viewChannels) {
                const ChannelID& channelID = channelPair.first;
                const vector<VoxelData>& channelVoxels = channelPair.second;
                
                ChannelData& channelData = entryChannels[channelID];
                
                for (const auto& voxel : channelVoxels) {
                    // Get charge for this view - conditionally apply uncalibration
                    double rawCharge;
                    if (g_useUncalibration) {
                        // Undo existing calibration to get raw charge
                        rawCharge = getUncalibratedCharge(voxel, view);
                    } else {
                        // Use charge directly (assume input is already raw)
                        rawCharge = getViewCharge(voxel, view);
                    }
                    
                    // Skip if we couldn't get valid charge
                    if (rawCharge <= 0) continue;
                    
                    // *** THIS IS WHERE WE CALCULATE DISTANCE ***
                    double distance = calculateVoxelToChannelDistance(voxel, channelID);
                    
                    // Apply attenuation correction to this voxel's charge
                    double correctedCharge = applyAttenuationCorrection(rawCharge, distance, view);
                    
                    // Add corrected charge to channel totals
                    channelData.totalCharge += correctedCharge;
                    channelData.totalPathLength += voxel.pathLength;
                    channelData.nVoxels++;
                }
            }
        }
        
        // Calculate gain for each channel in this entry
        for (const auto& channelPair : entryChannels) {
            const ChannelID& channelID = channelPair.first;
            const ChannelData& channelData = channelPair.second;
            
            if (channelData.totalPathLength > 0) {
                double gain = 3 * channelData.totalCharge / (2.0 * channelData.totalPathLength);
                channelGains[channelID].push_back(gain);
                channelGainDetails[channelID].push_back({entry, gain});
            }
        }
    }
    
    cout << "\nProcessing complete!" << endl;
    cout << "Total unique channels found: " << channelGains.size() << endl;
    
    // Create output directory and files
    string outputDir = createOutputDirectory(inputFile);
    string outputRoot = outputDir + "/ChannelGainAnalysis.root";
    string outputPNG = outputDir + "/PNG";
    
    gSystem->mkdir(outputPNG.c_str(), kTRUE);
    
    TFile* fOutput = new TFile(outputRoot.c_str(), "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        cerr << "Error: Cannot create output file: " << outputRoot << endl;
        fInput->Close();
        return;
    }
    
    // Set up style for histograms
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    
    // Calculate individual channel statistics for proper binning
    cout << "Calculating individual channel ranges..." << endl;
    
    // Create histograms and save gain details
    TTree* gainTree = new TTree("ChannelGainData", "Channel Gain Statistics");

    Char_t channel_name[100];
    Char_t view_name[10];
    Int_t coord1, coord2, n_entries;
    Double_t mean_gain, rms_gain, min_gain, max_gain, median_gain;
    Double_t fit_mean, fit_sigma, fit_chi2, fit_ndf;
    Bool_t fit_success;

    // Vectors to store all individual gain values for this channel
    vector<double> individual_gains;
    vector<int> file_entries;

    gainTree->Branch("channel_name", channel_name, "channel_name/C");
    gainTree->Branch("view", view_name, "view/C");
    gainTree->Branch("coord1", &coord1, "coord1/I");
    gainTree->Branch("coord2", &coord2, "coord2/I");
    gainTree->Branch("n_entries", &n_entries, "n_entries/I");
    gainTree->Branch("mean_gain", &mean_gain, "mean_gain/D");
    gainTree->Branch("rms_gain", &rms_gain, "rms_gain/D");
    gainTree->Branch("min_gain", &min_gain, "min_gain/D");
    gainTree->Branch("max_gain", &max_gain, "max_gain/D");
    gainTree->Branch("median_gain", &median_gain, "median_gain/D");
    gainTree->Branch("fit_mean", &fit_mean, "fit_mean/D");
    gainTree->Branch("fit_sigma", &fit_sigma, "fit_sigma/D");
    gainTree->Branch("fit_chi2", &fit_chi2, "fit_chi2/D");
    gainTree->Branch("fit_ndf", &fit_ndf, "fit_ndf/D");
    gainTree->Branch("fit_success", &fit_success, "fit_success/O");

    // Store individual measurements as vectors (optional - can be large)
    gainTree->Branch("individual_gains", &individual_gains);
    gainTree->Branch("file_entries", &file_entries);

    TCanvas* canvas = new TCanvas("canvas", "Channel Gain", 800, 600);

    int histCount = 0;

    for (const auto& channelPair : channelGains) {
        const ChannelID& channelID = channelPair.first;
        const vector<double>& gains = channelPair.second;
        
        if (histCount % 100 == 0) {
            cout << "Creating histogram " << histCount << "/" << channelGains.size() 
                << " for channel " << channelID.getName() << "\r" << flush;
        }
        
        // Filter valid gains (≤ 100.0)
        vector<double> validGains;
        vector<int> validFileEntries;
        
        for (size_t i = 0; i < gains.size(); i++) {
            if (gains[i] <= 200.0) {
                validGains.push_back(gains[i]);
                // Get corresponding file entry number
                validFileEntries.push_back(channelGainDetails[channelID][i].first);
            }
        }
        
        if (validGains.empty()) continue; // Skip channels with no valid gains
        
        // Calculate comprehensive statistics
        double sumGain = 0, sumGain2 = 0;
        for (double gain : validGains) {
            sumGain += gain;
            sumGain2 += gain * gain;
        }
        
        double meanGainVal = sumGain / validGains.size();
        double rmsGainVal = sqrt(sumGain2 / validGains.size() - meanGainVal * meanGainVal);
        
        // Calculate min, max, median
        vector<double> sortedGains = validGains;
        sort(sortedGains.begin(), sortedGains.end());
        
        double minGainVal = sortedGains.front();
        double maxGainVal = sortedGains.back();
        double medianGainVal;
        
        size_t n = sortedGains.size();
        if (n % 2 == 0) {
            medianGainVal = (sortedGains[n/2 - 1] + sortedGains[n/2]) / 2.0;
        } else {
            medianGainVal = sortedGains[n/2];
        }
        
        // Create histogram for fitting and visualization
        double histMin = 0.0;
        double histMax = 200.0;
        int nBins = 200;
        
        string histName = "h_gain_" + channelID.getName();
        string histTitle = "Gain Distribution - " + channelID.view + " Channel (" 
                        + to_string(channelID.coord1) + "," + to_string(channelID.coord2) + ")"
                        + ";Gain [PE/(2*TrackLength)];Entries";
        
        TH1F* hist = new TH1F(histName.c_str(), histTitle.c_str(), nBins, histMin, histMax);
        hist->SetFillColor(kTeal-5);       
        hist->SetLineColor(kTeal-5);         
        
        // Fill histogram
        for (double gain : validGains) {
            hist->Fill(gain);
        }
        
        // Initialize fit variables
        fit_mean = meanGainVal;
        fit_sigma = rmsGainVal;
        fit_chi2 = -1;
        fit_ndf = -1;
        fit_success = false;
        
        TF1* gaus = nullptr;
        
        // Try Gaussian fit if enough statistics
        if (validGains.size() >= 20 && hist->GetMaximum() > 5) {
            gaus = new TF1("gaus", "gaus", histMin, histMax);
            gaus->SetParameters(hist->GetMaximum(), meanGainVal, rmsGainVal);
            
            TFitResultPtr fitResult = hist->Fit(gaus, "QS");
            
            if (fitResult.Get() && fitResult->IsValid()) {
                fit_mean = gaus->GetParameter(1);
                fit_sigma = abs(gaus->GetParameter(2));
                fit_chi2 = fitResult->Chi2();
                fit_ndf = fitResult->Ndf();
                
                // Sanity check on fit results
                if (fit_mean > 0 && fit_mean < 100 && fit_sigma > 0 && fit_sigma < 50) {
                    fit_success = true;
                    gaus->SetLineColor(kRed);
                    gaus->SetLineWidth(2);
                } else {
                    fit_mean = meanGainVal;
                    fit_sigma = rmsGainVal;
                    fit_chi2 = -1;
                    fit_ndf = -1;
                    fit_success = false;
                }
            }
        }
        
        // Fill tree data - ONE ENTRY PER CHANNEL
        strcpy(channel_name, channelID.getName().c_str());
        strcpy(view_name, channelID.view.c_str());
        coord1 = channelID.coord1;
        coord2 = channelID.coord2;
        n_entries = validGains.size();
        mean_gain = meanGainVal;
        rms_gain = rmsGainVal;
        min_gain = minGainVal;
        max_gain = maxGainVal;
        median_gain = medianGainVal;
        
        // Store individual measurements (optional - comment out if tree gets too large)
        individual_gains = validGains;
        file_entries = validFileEntries;
        
        // Fill tree once per channel
        gainTree->Fill();
        
        // Clear vectors for next channel (important!)
        individual_gains.clear();
        file_entries.clear();
        
        // Save histogram and create PNG
        fOutput->cd();
        hist->Write();
        
        canvas->Clear();
        hist->Draw();
        
        if (gaus && fit_success) {
            gaus->Draw("same");
        }
        
        canvas->Update();
        
        string pngFile = outputPNG + "/gain_" + channelID.getName() + ".png";
        canvas->SaveAs(pngFile.c_str());
        
        // Clean up
        if (gaus) delete gaus;
        delete hist;
        histCount++;
    }
    
    cout << "\nCreated " << histCount << " histograms" << endl;
    
    // Save the gain tree
    fOutput->cd();
    gainTree->Write();
    
    // Create summary statistics with proper ranges
    // First find global range for summary histograms
    double globalMin = 1e9, globalMax = -1e9;
    vector<double> allMeans, allRMS;
    
    for (const auto& channelPair : channelGains) {
        const vector<double>& gains = channelPair.second;
        double sumGain = 0, sumGain2 = 0;
        int validGainCount = 0;
        
        for (double gain : gains) {
            if (gain <= 200.0) {
                sumGain += gain;
                sumGain2 += gain * gain;
                globalMin = min(globalMin, gain);
                globalMax = max(globalMax, gain);
                validGainCount++;
            }
        }
        
        if (validGainCount > 0) {
            double meanGain = sumGain / validGainCount;
            double rmsGain = sqrt(sumGain2 / validGainCount - meanGain * meanGain);
            allMeans.push_back(meanGain);
            allRMS.push_back(rmsGain);
        }
    }
    
    double meanMin = *min_element(allMeans.begin(), allMeans.end());
    double meanMax = *max_element(allMeans.begin(), allMeans.end());
    double rmsMax = *max_element(allRMS.begin(), allRMS.end());
    
    TH1F* hChannelCount = new TH1F("hChannelCount", "Number of Entries per Channel;Number of Entries;Channels", 50, 0, 50);
    TH1F* hMeanGain = new TH1F("hMeanGain", "Mean Gain Distribution;Mean Gain;Channels", 200, 0, 200);
    TH1F* hRMSGain = new TH1F("hRMSGain", "RMS Gain Distribution;RMS Gain;Channels", 200, 0, rmsMax * 1.1);

    map<string, int> viewCounts;
    

    
    fOutput->cd();
    hChannelCount->Write();
    hMeanGain->Write();
    hRMSGain->Write();
    
    // Print summary
    cout << "\n=== CHANNEL GAIN ANALYSIS SUMMARY ===" << endl;
    cout << "Total entries processed: " << nEntries << endl;
    cout << "Total unique channels: " << channelGains.size() << endl;
    cout << "Channels by view:" << endl;
    for (const auto& viewPair : viewCounts) {
        cout << "  " << viewPair.first << ": " << viewPair.second << " channels" << endl;
    }
    cout << "Global gain range: [" << globalMin << ", " << globalMax << "]" << endl;
    cout << "Mean gain range: [" << meanMin << ", " << meanMax << "]" << endl;
    cout << "\nOutput files:" << endl;
    cout << "  ROOT file: " << outputRoot << endl;
    cout << "  PNG files: " << outputPNG << "/" << endl;
    cout << "  Individual gain data saved in ChannelGainData tree" << endl;
    
    // Force write and close files properly
    fOutput->Write();
    fOutput->Close();
    
    fInput->Close();
    
    // Clean up heap-allocated objects
    delete canvas;
    
    exit(0);
}