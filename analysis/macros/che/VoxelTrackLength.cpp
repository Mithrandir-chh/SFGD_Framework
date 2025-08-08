// Trajectory Voxel Path Length Calculator
// Calculates the path length of particle trajectory through each detector voxel
#define THIS_NAME TrajectoryVoxelAnalyzer
#define OVERRIDE_OPTIONS

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <set>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TPolyLine3D.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "../../src/tools/global_header.hh"

using namespace std;

// Structure to hold voxel information
struct VoxelInfo {
    int x, y, z;
    double pathLength;
    bool hasHit;
    bool usedInAnalysis;
    double hitConfidence;
    double hitCharge;
    double hitChargeXY;
    double hitChargeXZ;
    double hitChargeZY;
    double entryT, exitT;
    bool onTrajectory;
    
    VoxelInfo() : x(0), y(0), z(0), pathLength(0), hasHit(false), 
                  usedInAnalysis(false), hitConfidence(0), hitCharge(0), 
                  hitChargeXY(0), hitChargeXZ(0), hitChargeZY(0), 
                  entryT(0), exitT(0), onTrajectory(false) {}
};

// Structure to hold trajectory parameters
struct TrajectoryParams {
    double x0, y0, z0;  // Starting point
    double vx, vy, vz;  // Direction vector
    bool valid;
    
    TrajectoryParams() : x0(0), y0(0), z0(0), vx(0), vy(0), vz(0), valid(false) {}
};

// Structure to hold hit information
struct HitInfo {
    double x, y, z;
    double confidence;
    double charge;
    double chargeXY;
    double chargeXZ;
    double chargeZY;
    
    HitInfo() : x(0), y(0), z(0), confidence(0), charge(0), chargeXY(0), chargeXZ(0),
    chargeZY(0) {}
};

// Structure to hold summary statistics for each file
struct FileSummary {
    string filename;
    int totalVoxels;
    int voxelsWithHits;
    int hitsOnTrajectory;
    int hitsOffTrajectory;
    double totalPathLength;
    double avgPathLength;
    bool valid;
    
    FileSummary() : totalVoxels(0), voxelsWithHits(0), hitsOnTrajectory(0), 
                   hitsOffTrajectory(0), totalPathLength(0), avgPathLength(0), valid(false) {}
};

// Function to get list of ROOT files in a directory
vector<string> getRootFiles(const string& directory) {
    vector<string> rootFiles;
    DIR* dir = opendir(directory.c_str());
    
    if (dir == nullptr) {
        cerr << "Error: Cannot open directory " << directory << endl;
        return rootFiles;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;
        // Check if file has .root extension
        if (filename.length() > 5 && 
            filename.substr(filename.length() - 5) == ".root") {
            rootFiles.push_back(directory + "/" + filename);
        }
    }
    
    closedir(dir);
    sort(rootFiles.begin(), rootFiles.end());
    
    // cout << "Found " << rootFiles.size() << " ROOT files in " << directory << endl;
    return rootFiles;
}

// Function to create output directory path
string getOutputPath(const string& inputPath) {
    // Find the last directory separator
    size_t lastSlash = inputPath.find_last_of("/\\");
    
    if (lastSlash == string::npos) {
        // No directory separator found, file is in current directory
        return "./";
    }
    
    // Extract directory path (same level as input file/folder)
    string inputDir = inputPath.substr(0, lastSlash);
    return inputDir + "/";
}

// Function to extract filename without path
string getFilename(const string& fullPath) {
    size_t lastSlash = fullPath.find_last_of("/\\");
    if (lastSlash == string::npos) {
        return fullPath;
    }
    return fullPath.substr(lastSlash + 1);
}

// Function to create directory if it doesn't exist
bool createDirectory(const string& path) {
    struct stat st = {0};
    if (stat(path.c_str(), &st) == -1) {
        return mkdir(path.c_str(), 0755) == 0;
    }
    return true; // Directory already exists
}

bool isHighConfidence(const HitInfo& hit, double threshold = 0.5) {
    return hit.confidence >= threshold;
}

// Load trajectory parameters from ROOT file
TrajectoryParams loadTrajectory(TFile* file) {
    TrajectoryParams traj;
    
    TTree* trajTree = (TTree*)file->Get("Trajectory");
    if (!trajTree) {
        // cout << "Error: Cannot find Trajectory tree" << endl;
        return traj;
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
    
    if (trajTree->GetEntries() > 0) {
        trajTree->GetEntry(0);
        traj.x0 = x0;
        traj.y0 = y0;
        traj.z0 = z0;
        traj.vx = vx;
        traj.vy = vy;
        traj.vz = vz;
        traj.valid = (valid == 1);
    }
    
    // cout << "Loaded trajectory: " << endl;
    // cout << "  Origin: (" << traj.x0 << ", " << traj.y0 << ", " << traj.z0 << ")" << endl;
    // cout << "  Direction: (" << traj.vx << ", " << traj.vy << ", " << traj.vz << ")" << endl;
    // cout << "  Valid: " << (traj.valid ? "Yes" : "No") << endl;
    
    return traj;
}

// Load hits from ROOT file
vector<HitInfo> loadAllHits(TFile* file) {  // Remove confidenceThreshold parameter
    vector<HitInfo> hits;
    
    TTree* hitsTree = (TTree*)file->Get("Hits3D");
    if (!hitsTree) {
        // cout << "Error: Cannot find Hits3D tree" << endl;
        return hits;
    }
    
    Double_t x, y, z, confidence, charge, chargeXY, chargeXZ, chargeZY;
    
    hitsTree->SetBranchAddress("x", &x);
    hitsTree->SetBranchAddress("y", &y);
    hitsTree->SetBranchAddress("z", &z);
    hitsTree->SetBranchAddress("confidence", &confidence);
    hitsTree->SetBranchAddress("charge", &charge);
    hitsTree->SetBranchAddress("chargeXY", &chargeXY);
    hitsTree->SetBranchAddress("chargeXZ", &chargeXZ);
    hitsTree->SetBranchAddress("chargeZY", &chargeZY);
    
    Long64_t nEntries = hitsTree->GetEntries();
    // cout << "Loading " << nEntries << " hits (all confidence levels)..." << endl;
    
    int highConfidenceCount = 0;
    double confidenceThreshold = 0.5;  // Keep threshold as local variable
    
    for (Long64_t i = 0; i < nEntries; i++) {
        hitsTree->GetEntry(i);
        
        HitInfo hit;
        hit.x = x;
        hit.y = y;
        hit.z = z;
        hit.confidence = confidence;
        hit.charge = charge;
        hit.chargeXY = chargeXY;
        hit.chargeXZ = chargeXZ;
        hit.chargeZY = chargeZY;
        hits.push_back(hit);  // Add ALL hits
        
        if (confidence >= confidenceThreshold) {
            highConfidenceCount++;
        }
    }
    
    // cout << "Loaded " << hits.size() << " total hits" << endl;
    // cout << "High-confidence hits (>= " << confidenceThreshold << "): " 
    //      << highConfidenceCount << endl;
    // cout << "Low-confidence hits: " << (hits.size() - highConfidenceCount) << endl;
    
    return hits;
}

// Calculate intersection of trajectory with detector boundaries
pair<double, double> findTrajectoryRange(const TrajectoryParams& traj) {
    double tMin = -1e9, tMax = 1e9;
    
    // Check intersection with x boundaries
    if (abs(traj.vx) > 1e-9) {
        double t1 = (0.0 - traj.x0) / traj.vx;
        double t2 = (24.0 - traj.x0) / traj.vx;
        tMin = max(tMin, min(t1, t2));
        tMax = min(tMax, max(t1, t2));
    } else {
        if (traj.x0 < 0.0 || traj.x0 > 24.0) {
            return make_pair(0, -1); // Invalid range
        }
    }
    
    // Check intersection with y boundaries
    if (abs(traj.vy) > 1e-9) {
        double t1 = (0.0 - traj.y0) / traj.vy;
        double t2 = (8.0 - traj.y0) / traj.vy;
        tMin = max(tMin, min(t1, t2));
        tMax = min(tMax, max(t1, t2));
    } else {
        if (traj.y0 < 0.0 || traj.y0 > 8.0) {
            return make_pair(0, -1); // Invalid range
        }
    }
    
    // Check intersection with z boundaries
    if (abs(traj.vz) > 1e-9) {
        double t1 = (0.0 - traj.z0) / traj.vz;
        double t2 = (48.0 - traj.z0) / traj.vz;
        tMin = max(tMin, min(t1, t2));
        tMax = min(tMax, max(t1, t2));
    } else {
        if (traj.z0 < 0.0 || traj.z0 > 48.0) {
            return make_pair(0, -1); // Invalid range
        }
    }
    
    // cout << "Trajectory parameter range: t = [" << tMin << ", " << tMax << "]" << endl;
    
    return make_pair(tMin, tMax);
}

// Get 3D position along trajectory
void getTrajectoryPoint(const TrajectoryParams& traj, double t, double& x, double& y, double& z) {
    x = traj.x0 + t * traj.vx;
    y = traj.y0 + t * traj.vy;
    z = traj.z0 + t * traj.vz;
}

// Calculate intersection of trajectory with voxel boundaries
vector<double> findVoxelIntersections(const TrajectoryParams& traj, int vx, int vy, int vz,
                                     double tMin, double tMax) {
    vector<double> intersections;
    
    // Voxel boundaries
    double xMin = vx, xMax = vx + 1;
    double yMin = vy, yMax = vy + 1;
    double zMin = vz, zMax = vz + 1;
    
    // Find intersections with x-planes
    if (abs(traj.vx) > 1e-9) {
        double t1 = (xMin - traj.x0) / traj.vx;
        double t2 = (xMax - traj.x0) / traj.vx;
        if (t1 >= tMin && t1 <= tMax) intersections.push_back(t1);
        if (t2 >= tMin && t2 <= tMax) intersections.push_back(t2);
    }
    
    // Find intersections with y-planes
    if (abs(traj.vy) > 1e-9) {
        double t1 = (yMin - traj.y0) / traj.vy;
        double t2 = (yMax - traj.y0) / traj.vy;
        if (t1 >= tMin && t1 <= tMax) intersections.push_back(t1);
        if (t2 >= tMin && t2 <= tMax) intersections.push_back(t2);
    }
    
    // Find intersections with z-planes
    if (abs(traj.vz) > 1e-9) {
        double t1 = (zMin - traj.z0) / traj.vz;
        double t2 = (zMax - traj.z0) / traj.vz;
        if (t1 >= tMin && t1 <= tMax) intersections.push_back(t1);
        if (t2 >= tMin && t2 <= tMax) intersections.push_back(t2);
    }
    
    // Sort intersections
    sort(intersections.begin(), intersections.end());
    
    // Remove duplicates
    intersections.erase(unique(intersections.begin(), intersections.end(),
                              [](double a, double b) { return abs(a - b) < 1e-9; }),
                       intersections.end());
    
    return intersections;
}

// Check if a point is inside a voxel
bool isInsideVoxel(double x, double y, double z, int vx, int vy, int vz) {
    return (x >= vx && x <= vx + 1 && 
            y >= vy && y <= vy + 1 && 
            z >= vz && z <= vz + 1);
}

// Calculate path length through a voxel
double calculateVoxelPathLength(const TrajectoryParams& traj, int vx, int vy, int vz,
                               double tMin, double tMax) {
    // Find all intersections with voxel boundaries
    vector<double> intersections = findVoxelIntersections(traj, vx, vy, vz, tMin, tMax);
    
    if (intersections.empty()) {
        // Check if entire trajectory segment is inside voxel
        double x1, y1, z1, x2, y2, z2;
        getTrajectoryPoint(traj, tMin, x1, y1, z1);
        getTrajectoryPoint(traj, tMax, x2, y2, z2);
        
        if (isInsideVoxel(x1, y1, z1, vx, vy, vz) && 
            isInsideVoxel(x2, y2, z2, vx, vy, vz)) {
            // Entire segment is inside voxel
            double dx = traj.vx * (tMax - tMin);
            double dy = traj.vy * (tMax - tMin);
            double dz = traj.vz * (tMax - tMin);
            return sqrt(dx*dx + dy*dy + dz*dz);
        }
        return 0.0;
    }
    
    // Add trajectory endpoints if they're inside the voxel
    double x, y, z;
    getTrajectoryPoint(traj, tMin, x, y, z);
    if (isInsideVoxel(x, y, z, vx, vy, vz)) {
        intersections.insert(intersections.begin(), tMin);
    }
    
    getTrajectoryPoint(traj, tMax, x, y, z);
    if (isInsideVoxel(x, y, z, vx, vy, vz)) {
        intersections.push_back(tMax);
    }
    
    // Sort and remove duplicates again
    sort(intersections.begin(), intersections.end());
    intersections.erase(unique(intersections.begin(), intersections.end(),
                              [](double a, double b) { return abs(a - b) < 1e-9; }),
                       intersections.end());
    
    if (intersections.size() < 2) return 0.0;
    
    // Calculate path length between valid intersection pairs
    double totalLength = 0.0;
    for (size_t i = 0; i < intersections.size() - 1; i++) {
        double t1 = intersections[i];
        double t2 = intersections[i + 1];
        double tMid = (t1 + t2) / 2.0;
        
        // Check if midpoint is inside voxel
        getTrajectoryPoint(traj, tMid, x, y, z);
        if (isInsideVoxel(x, y, z, vx, vy, vz)) {
            double dx = traj.vx * (t2 - t1);
            double dy = traj.vy * (t2 - t1);
            double dz = traj.vz * (t2 - t1);
            totalLength += sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    
    return totalLength;
}

// Comprehensive voxel traversal using plane intersections
vector<VoxelInfo> traverseVoxels(const TrajectoryParams& traj, const vector<HitInfo>& hits) {
    vector<VoxelInfo> voxels;
    
    // Find trajectory range within detector
    auto range = findTrajectoryRange(traj);
    double tMin = range.first;
    double tMax = range.second;
    
    if (tMax <= tMin) {
        // cout << "No valid trajectory range found in detector" << endl;
        return voxels;
    }
    
    // Create map for ALL hits (not just high confidence)
    map<tuple<int, int, int>, HitInfo> hitMap;
    int highConfidenceHits = 0;
    
    for (const auto& hit : hits) {
        int hx = (int)floor(hit.x);
        int hy = (int)floor(hit.y);
        int hz = (int)floor(hit.z);
        hitMap[make_tuple(hx, hy, hz)] = hit;
        
        if (isHighConfidence(hit)) {
            highConfidenceHits++;
        }
    }
    
    // cout << "Total hits in map: " << hitMap.size() << endl;
    // cout << "High-confidence hits for analysis: " << highConfidenceHits << endl;
    
    // Use fine sampling to find all intersected voxels (ONLY high-confidence hits)
    set<tuple<int, int, int>> intersectedVoxels;
    int nSamples = 100000;
    double dt = (tMax - tMin) / nSamples;
    
    for (int i = 0; i <= nSamples; i++) {
        double t = tMin + i * dt;
        double x, y, z;
        getTrajectoryPoint(traj, t, x, y, z);
        
        if (x >= 0 && x <= 24 && y >= 0 && y <= 8 && z >= 0 && z <= 48) {
            int vx = max(0, min(23, (int)floor(x)));
            int vy = max(0, min(7, (int)floor(y)));
            int vz = max(0, min(47, (int)floor(z)));
            
            intersectedVoxels.insert(make_tuple(vx, vy, vz));
        }
    }
    
    // cout << "Found " << intersectedVoxels.size() << " intersected voxels" << endl;
    
    // Calculate path length for each intersected voxel
    for (const auto& voxel : intersectedVoxels) {
        int vx = get<0>(voxel);
        int vy = get<1>(voxel);
        int vz = get<2>(voxel);
        
        VoxelInfo info;
        info.x = vx;
        info.y = vy;
        info.z = vz;
        info.pathLength = calculateVoxelPathLength(traj, vx, vy, vz, tMin, tMax);
        info.onTrajectory = true;
        
        // Check if this voxel has ANY hit (regardless of confidence)
        auto hitIt = hitMap.find(make_tuple(vx, vy, vz));
        if (hitIt != hitMap.end()) {
            info.hasHit = true;
            info.hitConfidence = hitIt->second.confidence;
            info.hitCharge = hitIt->second.charge;
            info.hitChargeXY = hitIt->second.chargeXY;
            info.hitChargeXZ = hitIt->second.chargeXZ;
            info.hitChargeZY = hitIt->second.chargeZY;
            info.usedInAnalysis = isHighConfidence(hitIt->second);  // NEW
        }
        
        if (info.pathLength > 1e-9) {
            voxels.push_back(info);
        }
    }
    
    // Add ALL hits that are NOT on the trajectory (including low confidence)
    for (const auto& hit : hits) {
        int hx = max(0, min(23, (int)floor(hit.x)));
        int hy = max(0, min(7, (int)floor(hit.y)));
        int hz = max(0, min(47, (int)floor(hit.z)));
        
        // Check if this voxel is already in our list
        bool alreadyIncluded = false;
        for (const auto& existingVoxel : voxels) {
            if (existingVoxel.x == hx && existingVoxel.y == hy && existingVoxel.z == hz) {
                alreadyIncluded = true;
                break;
            }
        }
        
        if (!alreadyIncluded) {
            VoxelInfo info;
            info.x = hx;
            info.y = hy;
            info.z = hz;
            info.pathLength = 0.0;
            info.hasHit = true;
            info.hitConfidence = hit.confidence;
            info.hitCharge = hit.charge;
            info.hitChargeXY = hit.chargeXY;
            info.hitChargeXZ = hit.chargeXZ;
            info.hitChargeZY = hit.chargeZY;
            info.onTrajectory = false;
            info.usedInAnalysis = isHighConfidence(hit);  // NEW
            
            voxels.push_back(info);
        }
    }
    
    // cout << "Total voxels in analysis: " << voxels.size() << " (including all hits)" << endl;
    
    return voxels;
}

// Single file analysis function (extracted from main)
bool analyzeFile(const string& inputFile, const string& outputDir) {
    // cout << "\n" << string(60, '=') << endl;
    // cout << "Processing file: " << inputFile << endl;
    // cout << string(60, '=') << endl;
    
    // Open input file
    TFile* file = TFile::Open(inputFile.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << inputFile << endl;
        return false;
    }
    
    // Load trajectory parameters
    TrajectoryParams traj = loadTrajectory(file);
    if (!traj.valid) {
        cerr << "Error: Invalid trajectory found in " << inputFile << endl;
        file->Close();
        return false;
    }
    
    // Load ALL hits
    vector<HitInfo> hits = loadAllHits(file);  // Remove confidence parameter
    
    // Perform voxel traversal
    // cout << "\nPerforming voxel traversal..." << endl;
    vector<VoxelInfo> voxels = traverseVoxels(traj, hits);
    
    // Create output file path
    string filename = getFilename(inputFile);
    string outputFilename = filename;
    size_t dotPos = outputFilename.find_last_of(".");
    if (dotPos != string::npos) {
        outputFilename = outputFilename.substr(0, dotPos) + "_VoxelPathLengths.root";
    } else {
        outputFilename += "_VoxelPathLengths.root";
    }
    
    string outputFile = outputDir + "/" + outputFilename;
    
    TFile* fOut = new TFile(outputFile.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        cerr << "Error: Cannot create output file " << outputFile << endl;
        file->Close();
        return false;
    }
    
    // Create histograms
    TH1F* hPathLength = new TH1F("hPathLength", "Path Length Distribution;Path Length [cm];Entries", 
                                100, 0, 2.0);
    TH1F* hPathLengthWithHits = new TH1F("hPathLengthWithHits", 
                                        "Path Length Distribution (Voxels with Hits);Path Length [cm];Entries",
                                        100, 0, 2.0);
    TH1F* hPathLengthNoHits = new TH1F("hPathLengthNoHits", 
                                      "Path Length Distribution (Voxels without Hits);Path Length [cm];Entries",
                                      100, 0, 2.0);
    
    // Separate histograms for hits on/off trajectory
    TH1F* hPathLengthHitsOnTrajectory = new TH1F("hPathLengthHitsOnTrajectory", 
                                                "Path Length Distribution (Hits on Trajectory);Path Length [cm];Entries",
                                                100, 0, 2.0);
    TH1F* hPathLengthHitsOffTrajectory = new TH1F("hPathLengthHitsOffTrajectory", 
                                                 "Path Length Distribution (Hits off Trajectory);Path Length [cm];Entries",
                                                 100, 0, 2.0);
    TH1F* hPathLengthAllHits = new TH1F("hPathLengthAllHits", 
                                       "Path Length Distribution (All Hits - Including Off-Track);Path Length [cm];Entries",
                                       100, 0, 2.0);
    
    TH3F* hVoxelMap = new TH3F("hVoxelMap", "Voxel Path Lengths with Trajectory;Z;Y;X", 
                               48, 0, 48, 8, 0, 8, 24, 0, 24);
    
    // Create trajectory visualization
    auto range = findTrajectoryRange(traj);
    double tMin = range.first;
    double tMax = range.second;
    
    TPolyLine3D* trajectoryLine = nullptr;
    if (tMax > tMin) {
        int nTrajPoints = 1000;
        trajectoryLine = new TPolyLine3D(nTrajPoints);
        
        for (int i = 0; i < nTrajPoints; i++) {
            double t = tMin + (tMax - tMin) * i / (nTrajPoints - 1);
            double x, y, z;
            getTrajectoryPoint(traj, t, x, y, z);
            
            // Only add points within detector bounds
            if (x >= 0 && x <= 24 && y >= 0 && y <= 8 && z >= 0 && z <= 48) {
                trajectoryLine->SetPoint(i, z, y, x);  // Note: ROOT uses Z,Y,X order for TH3
            }
        }
        
        trajectoryLine->SetLineColor(kRed);
        trajectoryLine->SetLineWidth(3);
        trajectoryLine->SetLineStyle(1);
    }
    
    // Create output tree
    TTree* voxelTree = new TTree("VoxelPathLengths", "Voxel Path Length Analysis");
    
    Int_t voxel_x, voxel_y, voxel_z;
    Double_t path_length, hit_confidence, hit_charge, hit_chargeXY, hit_chargeXZ, hit_chargeZY;
    Bool_t has_hit, on_trajectory, used_in_analysis;
    
    voxelTree->Branch("voxel_x", &voxel_x, "voxel_x/I");
    voxelTree->Branch("voxel_y", &voxel_y, "voxel_y/I");
    voxelTree->Branch("voxel_z", &voxel_z, "voxel_z/I");
    voxelTree->Branch("path_length", &path_length, "path_length/D");
    voxelTree->Branch("has_hit", &has_hit, "has_hit/O");
    voxelTree->Branch("hit_confidence", &hit_confidence, "hit_confidence/D");
    voxelTree->Branch("hit_charge", &hit_charge, "hit_charge/D");
    voxelTree->Branch("hit_chargeXY", &hit_chargeXY, "hit_chargeXY/D");    // NEW
    voxelTree->Branch("hit_chargeXZ", &hit_chargeXZ, "hit_chargeXZ/D");    // NEW
    voxelTree->Branch("hit_chargeZY", &hit_chargeZY, "hit_chargeZY/D");  
    voxelTree->Branch("used_in_analysis", &used_in_analysis, "used_in_analysis/O");
    voxelTree->Branch("on_trajectory", &on_trajectory, "on_trajectory/O");
    
    // Fill histograms and tree
    double totalPathLength = 0;
    int voxelsWithHits = 0;
    int voxelsWithoutHits = 0;
    int hitsOnTrajectory = 0;
    int hitsOffTrajectory = 0;
    
    for (const auto& voxel : voxels) {
        // Always fill the tree (ALL hits)
        voxel_x = voxel.x;
        voxel_y = voxel.y;
        voxel_z = voxel.z;
        path_length = voxel.pathLength;
        has_hit = voxel.hasHit;
        hit_confidence = voxel.hitConfidence;
        hit_charge = voxel.hitCharge;
        hit_chargeXY = voxel.hitChargeXY;
        hit_chargeXZ = voxel.hitChargeXZ;
        hit_chargeZY = voxel.hitChargeZY;
        on_trajectory = voxel.onTrajectory;
        used_in_analysis = voxel.usedInAnalysis;  // NEW
        
        voxelTree->Fill();  // ALL hits go into tree
        
        // Only fill histograms with high-confidence hits OR trajectory voxels
        bool fillHistograms = !voxel.hasHit || voxel.usedInAnalysis;
        
        if (fillHistograms) {
            hPathLength->Fill(path_length);
            
            if (voxel.hasHit && voxel.usedInAnalysis) {  // Only high-confidence hits
                voxelsWithHits++;
                hPathLengthWithHits->Fill(path_length);
                hPathLengthAllHits->Fill(path_length);
                
                if (voxel.onTrajectory) {
                    hitsOnTrajectory++;
                    hPathLengthHitsOnTrajectory->Fill(path_length);
                } else {
                    hitsOffTrajectory++;
                    hPathLengthHitsOffTrajectory->Fill(path_length);
                }
            } else {
                hPathLengthNoHits->Fill(path_length);
                voxelsWithoutHits++;
            }
            
            // Fill 3D histogram
            if (voxel.pathLength > 0) {
                hVoxelMap->SetBinContent(voxel.z + 1, voxel.y + 1, voxel.x + 1, path_length);
            }
        }
        
        totalPathLength += path_length;
    }
    
    // Write results
    hPathLength->Write();
    hPathLengthWithHits->Write();
    hPathLengthNoHits->Write();
    hPathLengthHitsOnTrajectory->Write();    
    hPathLengthHitsOffTrajectory->Write();    
    hPathLengthAllHits->Write();              
    hVoxelMap->Write();
    voxelTree->Write();
    
    // Write trajectory
    if (trajectoryLine) {
        trajectoryLine->Write("trajectory");
    }
    
    // Print summary
    // cout << "\n=== Voxel Path Length Analysis Summary ===" << endl;
    // cout << "Total voxels traversed by trajectory: " << (voxels.size() - hitsOffTrajectory) << endl;
    // cout << "Total voxels in analysis: " << voxels.size() << endl;
    // cout << "Voxels with hits: " << voxelsWithHits << endl;
    // cout << "  - Hits on trajectory: " << hitsOnTrajectory << endl;
    // cout << "  - Hits off trajectory: " << hitsOffTrajectory << endl;
    // cout << "Voxels without hits: " << voxelsWithoutHits << endl;
    // cout << "Total path length: " << totalPathLength << " cm" << endl;
    // cout << "Average path length per voxel: " << (voxels.size() > 0 ? totalPathLength / voxels.size() : 0) << " cm" << endl;
    
    if (hitsOnTrajectory > 0) {
        double avgPathLengthOnTrajectory = 0;
        for (const auto& voxel : voxels) {
            if (voxel.hasHit && voxel.onTrajectory) avgPathLengthOnTrajectory += voxel.pathLength;
        }
        avgPathLengthOnTrajectory /= hitsOnTrajectory;
        // cout << "Average path length (hits on trajectory): " << avgPathLengthOnTrajectory << " cm" << endl;
    }
    
    // cout << "\nOutput saved to: " << outputFile << endl;
    
    fOut->Close();
    file->Close();
    
    return true;
}

void createSummaryFile(const vector<string>& inputFiles, const string& outputDir, 
                      const vector<FileSummary>& summaries) {
    
    // cout << "\n" << string(60, '=') << endl;
    // cout << "Creating consolidated summary ROOT file..." << endl;
    // cout << string(60, '=') << endl;
    
    string summaryFile = outputDir + "/A_Summary_AllFiles.root";
    TFile* fSummary = new TFile(summaryFile.c_str(), "RECREATE");
    
    if (!fSummary || fSummary->IsZombie()) {
        cerr << "Error: Cannot create summary file " << summaryFile << endl;
        return;
    }
    
    fSummary->SetCompressionLevel(1);
    
    // === NEW: Consolidated hits tree ===
    TTree* consolidatedHits = new TTree("ConsolidatedHits", "All Hits from All Files");
    
    // Vectors to store ALL hits data
    vector<Int_t> file_ids;
    vector<Int_t> voxel_xs, voxel_ys, voxel_zs;
    vector<Double_t> charges, confidences;
    vector<Double_t> chargeXYs, chargeXZs, chargeZYs;
    vector<Double_t> path_lengths;
    vector<Bool_t> has_hits, on_trajectories, used_in_analyses;
    
    // Create branches for vectors
    consolidatedHits->Branch("file_id", &file_ids);
    consolidatedHits->Branch("voxel_x", &voxel_xs);
    consolidatedHits->Branch("voxel_y", &voxel_ys);
    consolidatedHits->Branch("voxel_z", &voxel_zs);
    consolidatedHits->Branch("charge", &charges);
    consolidatedHits->Branch("confidence", &confidences);
    consolidatedHits->Branch("chargeXY", &chargeXYs);
    consolidatedHits->Branch("chargeXZ", &chargeXZs);
    consolidatedHits->Branch("chargeZY", &chargeZYs);
    consolidatedHits->Branch("path_length", &path_lengths);
    consolidatedHits->Branch("has_hit", &has_hits);
    consolidatedHits->Branch("on_trajectory", &on_trajectories);
    consolidatedHits->Branch("used_in_analysis", &used_in_analyses);
    
    // === NEW: File mapping tree ===
    TTree* fileMapping = new TTree("FileMapping", "File ID to Filename Mapping");
    
    Int_t file_id_map;
    Char_t file_name_map[256];
    Int_t total_hits_in_file;
    
    fileMapping->Branch("file_id", &file_id_map, "file_id/I");
    fileMapping->Branch("file_name", file_name_map, "file_name/C");
    fileMapping->Branch("total_hits_in_file", &total_hits_in_file, "total_hits_in_file/I");
    
    // [Keep existing histogram creation code...]
    TH1F* hPathLengthAll = new TH1F("hPathLengthAll", 
                                   "Combined Path Length Distribution (All Files);Path Length [cm];Entries", 
                                   100, 0, 2.0);
    TH1F* hPathLengthWithHitsAll = new TH1F("hPathLengthWithHitsAll", 
                                           "Combined Path Length Distribution - All Hits (All Files);Path Length [cm];Entries",
                                           100, 0, 2.0);
    TH1F* hPathLengthHitsOnTrajectoryAll = new TH1F("hPathLengthHitsOnTrajectoryAll", 
                                                   "Combined Path Length Distribution - Hits on Trajectory (All Files);Path Length [cm];Entries",
                                                   100, 0, 2.0);
    TH1F* hPathLengthHitsOffTrajectoryAll = new TH1F("hPathLengthHitsOffTrajectoryAll", 
                                                    "Combined Path Length Distribution - Hits off Trajectory (All Files);Path Length [cm];Entries",
                                                    100, 0, 2.0);
    
    // 3D combined histogram for all trajectories and hits
    TH3F* hVoxelMapAll = new TH3F("hVoxelMapAll", 
                                 "Combined Voxel Map - All Trajectories;Z;Y;X", 
                                 48, 0, 48, 8, 0, 8, 24, 0, 24);
                        
    TH3F* hVoxelMapFiltered = new TH3F("hVoxelMapFiltered", 
                                      "Filtered Hits: On-Track + High Confidence + Path Length >= 0.7cm;Z;Y;X", 
                                      48, 0, 48, 8, 0, 8, 24, 0, 24);
    
    // Summary statistics histograms
    int nFiles = 0;
    for (const auto& summary : summaries) {
        if (summary.valid) nFiles++;
    }
    
    TH1F* hTotalVoxelsPerFile = new TH1F("hTotalVoxelsPerFile", 
                                        "Total Voxels per File;File Number;Total Voxels", 
                                        nFiles, 0, nFiles);
    TH1F* hHitsOnTrajectoryPerFile = new TH1F("hHitsOnTrajectoryPerFile", 
                                             "Hits on Trajectory per File;File Number;Hits on Trajectory", 
                                             nFiles, 0, nFiles);
    TH1F* hHitsOffTrajectoryPerFile = new TH1F("hHitsOffTrajectoryPerFile", 
                                              "Hits off Trajectory per File;File Number;Hits off Trajectory", 
                                              nFiles, 0, nFiles);
    TH1F* hTotalPathLengthPerFile = new TH1F("hTotalPathLengthPerFile", 
                                            "Total Path Length per File;File Number;Total Path Length [cm]", 
                                            nFiles, 0, nFiles);
    TH1F* hAvgPathLengthPerFile = new TH1F("hAvgPathLengthPerFile", 
                                          "Average Path Length per File;File Number;Average Path Length [cm]", 
                                          nFiles, 0, nFiles);
    
    // Collection to store all trajectory lines
    vector<TPolyLine3D*> allTrajectories;
    
    // Process each successful file
    int fileIndex = 0;
    int validFiles = 0;
    int totalHitsAllFiles = 0;

    for (size_t i = 0; i < inputFiles.size(); i++) {
        const FileSummary& summary = summaries[i];
        
        if (!summary.valid) continue;
        
        // Get corresponding output file
        string inputFilename = getFilename(inputFiles[i]);
        string outputFilename = inputFilename;
        size_t dotPos = outputFilename.find_last_of(".");
        if (dotPos != string::npos) {
            outputFilename = outputFilename.substr(0, dotPos) + "_VoxelPathLengths.root";
        } else {
            outputFilename += "_VoxelPathLengths.root";
        }
        
        string outputFile = outputDir + "/" + outputFilename;
        
        // Open the individual analysis file
        TFile* fIndividual = TFile::Open(outputFile.c_str(), "READ");
        if (!fIndividual || fIndividual->IsZombie()) {
            cout << "Warning: Cannot open " << outputFile << " for summary" << endl;
            continue;
        }
        
        // === CLEAR vectors for this file ===
        file_ids.clear();
        voxel_xs.clear(); voxel_ys.clear(); voxel_zs.clear();
        charges.clear(); confidences.clear();
        chargeXYs.clear(); chargeXZs.clear(); chargeZYs.clear();
        path_lengths.clear();
        has_hits.clear(); on_trajectories.clear(); used_in_analyses.clear();
        
        // === Read all voxel data from THIS file ===
        TTree* voxelTree = (TTree*)fIndividual->Get("VoxelPathLengths");
        if (voxelTree) {
            Int_t voxel_x, voxel_y, voxel_z;
            Double_t path_length, hit_confidence, hit_charge;
            Double_t hit_chargeXY, hit_chargeXZ, hit_chargeZY;
            Bool_t has_hit, on_trajectory, used_in_analysis;
            
            // Set branch addresses
            voxelTree->SetBranchAddress("voxel_x", &voxel_x);
            voxelTree->SetBranchAddress("voxel_y", &voxel_y);
            voxelTree->SetBranchAddress("voxel_z", &voxel_z);
            voxelTree->SetBranchAddress("path_length", &path_length);
            voxelTree->SetBranchAddress("hit_confidence", &hit_confidence);
            voxelTree->SetBranchAddress("hit_charge", &hit_charge);
            voxelTree->SetBranchAddress("hit_chargeXY", &hit_chargeXY);
            voxelTree->SetBranchAddress("hit_chargeXZ", &hit_chargeXZ);
            voxelTree->SetBranchAddress("hit_chargeZY", &hit_chargeZY);
            voxelTree->SetBranchAddress("has_hit", &has_hit);
            voxelTree->SetBranchAddress("on_trajectory", &on_trajectory);
            voxelTree->SetBranchAddress("used_in_analysis", &used_in_analysis);
            
            Long64_t nEntries = voxelTree->GetEntries();
            int filteredCount = 0;
            // cout << "Reading " << nEntries << " voxels from " << inputFilename << endl;
            
            // Read all entries for THIS file and add to vectors
            for (Long64_t entry = 0; entry < nEntries; entry++) {
                voxelTree->GetEntry(entry);
                
                // Add to vectors for this file
                file_ids.push_back(fileIndex);
                voxel_xs.push_back(voxel_x);
                voxel_ys.push_back(voxel_y);
                voxel_zs.push_back(voxel_z);
                charges.push_back(hit_charge);
                confidences.push_back(hit_confidence);
                chargeXYs.push_back(hit_chargeXY);
                chargeXZs.push_back(hit_chargeXZ);
                chargeZYs.push_back(hit_chargeZY);
                path_lengths.push_back(path_length);
                has_hits.push_back(has_hit);
                on_trajectories.push_back(on_trajectory);
                used_in_analyses.push_back(used_in_analysis);

                if (on_trajectory && has_hit && used_in_analysis && path_length >= 0.7) {
                    // Fill the filtered histogram (note: ROOT uses z,y,x convention)
                    hVoxelMapFiltered->SetBinContent(voxel_z + 1, voxel_y + 1, voxel_x + 1, path_length);
                    filteredCount++;
                }
            }
            
            totalHitsAllFiles += nEntries;
            
            // Fill file mapping
            file_id_map = fileIndex;
            strncpy(file_name_map, inputFilename.c_str(), 255);
            file_name_map[255] = '\0';
            total_hits_in_file = nEntries;
            fileMapping->Fill();
            
            // === FILL consolidated tree for THIS file ===
            consolidatedHits->Fill();
            
            // cout << "Added entry " << fileIndex << " to ConsolidatedHits with " 
            //     << nEntries << " voxels" << endl;
        }
        
        // Add histograms to combined histograms
        TH1F* hPath = (TH1F*)fIndividual->Get("hPathLength");
        TH1F* hPathWithHits = (TH1F*)fIndividual->Get("hPathLengthWithHits");
        TH1F* hPathOnTraj = (TH1F*)fIndividual->Get("hPathLengthHitsOnTrajectory");
        TH1F* hPathOffTraj = (TH1F*)fIndividual->Get("hPathLengthHitsOffTrajectory");
        TH3F* hVoxelMap = (TH3F*)fIndividual->Get("hVoxelMap");
        
        if (hPath) hPathLengthAll->Add(hPath);
        if (hPathWithHits) hPathLengthWithHitsAll->Add(hPathWithHits);
        if (hPathOnTraj) hPathLengthHitsOnTrajectoryAll->Add(hPathOnTraj);
        if (hPathOffTraj) hPathLengthHitsOffTrajectoryAll->Add(hPathOffTraj);
        if (hVoxelMap) hVoxelMapAll->Add(hVoxelMap);
        
        // Get trajectory line and add to collection
        TPolyLine3D* traj = (TPolyLine3D*)fIndividual->Get("trajectory");
        if (traj) {
            TPolyLine3D* trajCopy = new TPolyLine3D(*traj);
            trajCopy->SetLineColor(kRed + (validFiles % 8));  // Different colors for different files
            trajCopy->SetLineWidth(2);
            trajCopy->SetLineStyle(1);
            allTrajectories.push_back(trajCopy);
        }
        
        // Fill summary histograms
        hTotalVoxelsPerFile->SetBinContent(fileIndex + 1, summary.totalVoxels);
        hTotalVoxelsPerFile->GetXaxis()->SetBinLabel(fileIndex + 1, inputFilename.c_str());
        
        hHitsOnTrajectoryPerFile->SetBinContent(fileIndex + 1, summary.hitsOnTrajectory);
        hHitsOnTrajectoryPerFile->GetXaxis()->SetBinLabel(fileIndex + 1, inputFilename.c_str());
        
        hHitsOffTrajectoryPerFile->SetBinContent(fileIndex + 1, summary.hitsOffTrajectory);
        hHitsOffTrajectoryPerFile->GetXaxis()->SetBinLabel(fileIndex + 1, inputFilename.c_str());
        
        hTotalPathLengthPerFile->SetBinContent(fileIndex + 1, summary.totalPathLength);
        hTotalPathLengthPerFile->GetXaxis()->SetBinLabel(fileIndex + 1, inputFilename.c_str());
        
        hAvgPathLengthPerFile->SetBinContent(fileIndex + 1, summary.avgPathLength);
        hAvgPathLengthPerFile->GetXaxis()->SetBinLabel(fileIndex + 1, inputFilename.c_str());
        
        fIndividual->Close();
        delete fIndividual;  // Explicitly delete to free memory
        fileIndex++;
        validFiles++;
    }
    
    // === Fill the consolidated tree (single entry with all vectors) ===
    consolidatedHits->Fill();
    
    // cout << "\nConsolidated " << totalHitsAllFiles << " total voxels from " 
    //      << fileIndex << " files" << endl;

    // Create summary tree
    TTree* summaryTree = new TTree("SummaryStats", "Summary Statistics for All Files");
    
    Int_t file_number;
    Char_t file_name[256];  // Use fixed-size char array instead of TString
    Int_t total_voxels, voxels_with_hits, hits_on_trajectory, hits_off_trajectory;
    Double_t total_path_length, avg_path_length;
    
    summaryTree->Branch("file_number", &file_number, "file_number/I");
    summaryTree->Branch("file_name", file_name, "file_name/C");
    summaryTree->Branch("total_voxels", &total_voxels, "total_voxels/I");
    summaryTree->Branch("voxels_with_hits", &voxels_with_hits, "voxels_with_hits/I");
    summaryTree->Branch("hits_on_trajectory", &hits_on_trajectory, "hits_on_trajectory/I");
    summaryTree->Branch("hits_off_trajectory", &hits_off_trajectory, "hits_off_trajectory/I");
    summaryTree->Branch("total_path_length", &total_path_length, "total_path_length/D");
    summaryTree->Branch("avg_path_length", &avg_path_length, "avg_path_length/D");
    
    // Fill summary tree
    int validIndex = 0;
    for (size_t i = 0; i < summaries.size(); i++) {
        if (!summaries[i].valid) continue;
        
        file_number = validIndex;
        strncpy(file_name, getFilename(inputFiles[i]).c_str(), 255);
        file_name[255] = '\0';  // Ensure null termination
        total_voxels = summaries[i].totalVoxels;
        voxels_with_hits = summaries[i].voxelsWithHits;
        hits_on_trajectory = summaries[i].hitsOnTrajectory;
        hits_off_trajectory = summaries[i].hitsOffTrajectory;
        total_path_length = summaries[i].totalPathLength;
        avg_path_length = summaries[i].avgPathLength;
        
        summaryTree->Fill();
        validIndex++;
    }

    // Change to the summary file's directory before writing
    fSummary->cd();
    
    // Write all histograms first
    hPathLengthAll->Write();
    hPathLengthWithHitsAll->Write();
    hPathLengthHitsOnTrajectoryAll->Write();
    hPathLengthHitsOffTrajectoryAll->Write();
    hVoxelMapAll->Write();
    hVoxelMapFiltered->Write();
    
    hTotalVoxelsPerFile->Write();
    hHitsOnTrajectoryPerFile->Write();
    hHitsOffTrajectoryPerFile->Write();
    hTotalPathLengthPerFile->Write();
    hAvgPathLengthPerFile->Write();
    
    // Write TTrees
    consolidatedHits->Write();
    fileMapping->Write();
    summaryTree->Write();
    
    // Write all trajectory lines with proper naming
    for (size_t i = 0; i < allTrajectories.size(); i++) {
        string trajName = Form("trajectory_file_%d", (int)i);
        allTrajectories[i]->Write(trajName.c_str());
    }
    
    // Force write of all file data to disk
    fSummary->Write();
    fSummary->Flush();
    
    // Calculate and print summary statistics
    // cout << "\n=== SUMMARY STATISTICS FOR ALL FILES ===" << endl;
    // cout << "Number of valid files processed: " << validFiles << endl;
    
    if (validFiles > 0) {
        int totalVoxelsAll = 0, totalHitsOnTraj = 0, totalHitsOffTraj = 0;
        double totalPathLengthAll = 0;
        
        for (const auto& summary : summaries) {
            if (summary.valid) {
                totalVoxelsAll += summary.totalVoxels;
                totalHitsOnTraj += summary.hitsOnTrajectory;
                totalHitsOffTraj += summary.hitsOffTrajectory;
                totalPathLengthAll += summary.totalPathLength;
            }
        }
        
        // cout << "Combined statistics:" << endl;
        // cout << "  Total voxels: " << totalVoxelsAll << endl;
        // cout << "  Total hits on trajectory: " << totalHitsOnTraj << endl;
        // cout << "  Total hits off trajectory: " << totalHitsOffTraj << endl;
        // cout << "  Total path length: " << totalPathLengthAll << " cm" << endl;
        // cout << "  Average path length across all files: " << (totalVoxelsAll > 0 ? totalPathLengthAll / totalVoxelsAll : 0) << " cm" << endl;
        
        if (totalHitsOnTraj + totalHitsOffTraj > 0) {
            // cout << "  Hit efficiency (on-track): " << (100.0 * totalHitsOnTraj / (totalHitsOnTraj + totalHitsOffTraj)) << "%" << endl;
        }
    }
    
    // Clean up dynamically allocated trajectories
    for (auto* traj : allTrajectories) {
        delete traj;
    }
    allTrajectories.clear();
    
    // Close the file
    fSummary->Close();
    delete fSummary;
    
    // cout << "\nSummary file saved to: " << summaryFile << endl;
    // cout << "\nTo view summary results in ROOT:" << endl;
    // cout << "  root " << summaryFile << endl;
    // cout << "  // Combined distributions:" << endl;
    // cout << "  hPathLengthAll->Draw()  // All path lengths combined" << endl;
    // cout << "  hPathLengthWithHitsAll->Draw()  // All hits combined" << endl;
    // cout << "  hPathLengthHitsOnTrajectoryAll->Draw()  // All on-track hits" << endl;
    // cout << "  hPathLengthHitsOffTrajectoryAll->Draw()  // All off-track hits" << endl;
    // cout << "  " << endl;
    // cout << "  // Per-file statistics:" << endl;
    // cout << "  hTotalVoxelsPerFile->Draw()  // Voxels per file" << endl;
    // cout << "  hHitsOnTrajectoryPerFile->Draw()  // On-track hits per file" << endl;
    // cout << "  hHitsOffTrajectoryPerFile->Draw()  // Off-track hits per file" << endl;
    // cout << "  " << endl;
    // cout << "  // 3D visualization with all trajectories:" << endl;
    // cout << "  hVoxelMapAll->Draw(\"BOX\")  // All voxels combined" << endl;
    // cout << "  TPolyLine3D *traj0 = (TPolyLine3D*)gDirectory->Get(\"trajectory_file_0\")" << endl;
    // cout << "  if(traj0) traj0->Draw()  // Draw first trajectory" << endl;
    // cout << "  TPolyLine3D *traj1 = (TPolyLine3D*)gDirectory->Get(\"trajectory_file_1\")" << endl;
    // cout << "  if(traj1) traj1->Draw()  // Draw second trajectory (if exists)" << endl;
    // cout << "  // ... and so on for each file" << endl;
    // cout << "  " << endl;
    // cout << "  // Summary tree analysis:" << endl;
    // cout << "  SummaryStats->Draw(\"hits_on_trajectory:file_number\")  // Hits vs file" << endl;
    // cout << "  SummaryStats->Scan(\"*\")  // Print all summary data" << endl;
    // cout << "\nTo analyze consolidated data in ROOT:" << endl;
    // cout << "  root " << summaryFile << endl;
    // cout << "  // Access all data:" << endl;
    // cout << "  ConsolidatedHits->Draw(\"charge\")  // All charges from all files" << endl;
    // cout << "  ConsolidatedHits->Draw(\"chargeXY\")  // All XY projections" << endl;
    // cout << "  ConsolidatedHits->Draw(\"confidence\", \"file_id==0\")  // Only first file" << endl;
    // cout << "  ConsolidatedHits->Draw(\"charge:file_id\")  // Charge vs file" << endl;
    // cout << "  ConsolidatedHits->Draw(\"path_length\", \"has_hit && used_in_analysis\")  // High-confidence hits only" << endl;
    // cout << "  " << endl;
    // cout << "  // File mapping:" << endl;
    // cout << "  FileMapping->Scan(\"*\")  // See file ID to name mapping" << endl;
}

// Main analysis function
void TrajectoryVoxelAnalyzer() {
    gROOT->SetBatch(kTRUE);
    
    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_3D_file.root_OR_input_directory>" << endl;
        cerr << "  If a single .root file is provided, only that file will be processed." << endl;
        cerr << "  If a directory is provided, all .root files in that directory will be processed." << endl;
        return;
    }
    
    string inputPath = argv[1];
    // cout << "Input path: " << inputPath << endl;
    
    // Check if input is a file or directory
    struct stat pathStat;
    if (stat(inputPath.c_str(), &pathStat) != 0) {
        cerr << "Error: Cannot access " << inputPath << endl;
        return;
    }
    
    vector<string> inputFiles;
    string outputDir;
    
    if (S_ISDIR(pathStat.st_mode)) {
        // Input is a directory
        // cout << "Processing directory: " << inputPath << endl;
        inputFiles = getRootFiles(inputPath);
        
        if (inputFiles.empty()) {
            cerr << "No ROOT files found in directory " << inputPath << endl;
            return;
        }
        
        // Create output directory one level up from input directory
        outputDir = getOutputPath(inputPath);
        
        // Extract the name of the input directory to create a results folder
        size_t lastSlash = inputPath.find_last_of("/\\");
        string inputDirName;
        if (lastSlash == string::npos) {
            inputDirName = inputPath;
        } else {
            inputDirName = inputPath.substr(lastSlash + 1);
        }
        
        outputDir += inputDirName + "_VoxelAnalysis";
        
        // Create output directory
        if (!createDirectory(outputDir)) {
            cerr << "Error: Cannot create output directory " << outputDir << endl;
            return;
        }
        
        // cout << "Output directory: " << outputDir << endl;
        
    } else if (S_ISREG(pathStat.st_mode)) {
        // Input is a single file
        // cout << "Processing single file: " << inputPath << endl;
        
        // Check if it's a ROOT file
        if (inputPath.length() < 5 || inputPath.substr(inputPath.length() - 5) != ".root") {
            cerr << "Error: Input file must have .root extension" << endl;
            return;
        }
        
        inputFiles.push_back(inputPath);
        string fileDir = getOutputPath(inputPath);
        outputDir = getOutputPath(fileDir) + "VoxelAnalysis";
        
        // Create output directory
        if (!createDirectory(outputDir)) {
            cerr << "Error: Cannot create output directory " << outputDir << endl;
            return;
        }
        
        // cout << "Output directory: " << outputDir << endl;
        
    } else {
        cerr << "Error: Input path is neither a regular file nor a directory" << endl;
        return;
    }
    
    // Process all files
    int successCount = 0;
    int totalFiles = inputFiles.size();
    vector<FileSummary> summaries(totalFiles);
    
    // cout << "\n" << string(80, '=') << endl;
    // cout << "Starting batch analysis of " << totalFiles << " files" << endl;
    // cout << string(80, '=') << endl;
    
    for (int i = 0; i < totalFiles; i++) {
        // cout << "\nProcessing file " << (i + 1) << "/" << totalFiles << ": " << getFilename(inputFiles[i]) << endl;
        
        if (analyzeFile(inputFiles[i], outputDir)) {
            successCount++;
            // cout << "✓ Successfully processed: " << getFilename(inputFiles[i]) << endl;
            
            // Extract summary statistics from the processed file
            string inputFilename = getFilename(inputFiles[i]);
            string outputFilename = inputFilename;
            size_t dotPos = outputFilename.find_last_of(".");
            if (dotPos != string::npos) {
                outputFilename = outputFilename.substr(0, dotPos) + "_VoxelPathLengths.root";
            } else {
                outputFilename += "_VoxelPathLengths.root";
            }
            
            string outputFile = outputDir + "/" + outputFilename;
            
            // Read summary statistics from the output file
            TFile* fTemp = TFile::Open(outputFile.c_str(), "READ");
            if (fTemp && !fTemp->IsZombie()) {
                TTree* voxelTree = (TTree*)fTemp->Get("VoxelPathLengths");
                if (voxelTree) {
                    summaries[i].filename = inputFilename;
                    summaries[i].totalVoxels = voxelTree->GetEntries();
                    summaries[i].voxelsWithHits = voxelTree->GetEntries("has_hit");
                    summaries[i].hitsOnTrajectory = voxelTree->GetEntries("has_hit && on_trajectory");
                    summaries[i].hitsOffTrajectory = voxelTree->GetEntries("has_hit && !on_trajectory");
                    
                    // Calculate total path length
                    TH1F* tempHist = new TH1F("temp", "temp", 1000, 0, 1000);
                    voxelTree->Draw("path_length>>temp", "", "goff");
                    summaries[i].totalPathLength = tempHist->GetMean() * tempHist->GetEntries();
                    summaries[i].avgPathLength = tempHist->GetMean();
                    summaries[i].valid = true;
                    delete tempHist;
                }
                fTemp->Close();
            }
        } else {
            cout << "✗ Failed to process: " << getFilename(inputFiles[i]) << endl;
            summaries[i].valid = false;
        }
    }
    
    // Create summary file if we have multiple successful files
    if (successCount > 1) {
        createSummaryFile(inputFiles, outputDir, summaries);
    } else if (successCount == 1) {
        cout << "\nNote: Only one file was successfully processed. Summary file not created." << endl;
        cout << "Summary files are created when multiple files are processed." << endl;
    }
    
    // Final summary
    cout << "\n" << string(80, '=') << endl;
    cout << "BATCH ANALYSIS COMPLETE" << endl;
    cout << string(80, '=') << endl;
    cout << "Total files processed: " << successCount << "/" << totalFiles << endl;
    cout << "Success rate: " << (totalFiles > 0 ? (100.0 * successCount / totalFiles) : 0) << "%" << endl;
    cout << "Output directory: " << outputDir << endl;
    
    if (successCount > 0) {
        cout << "\nTo view individual results in ROOT:" << endl;
        cout << "  cd " << outputDir << endl;
        cout << "  root <filename>_VoxelPathLengths.root" << endl;
        cout << "  // Then use any of the following commands:" << endl;
        cout << "  hPathLength->Draw()  // All voxels on track + All voxels with detected hits" << endl;
        cout << "  hPathLengthWithHits->Draw()  // All 3D voxels with detected hits (including off-track)" << endl;
        cout << "  hPathLengthHitsOnTrajectory->Draw()   // Only voxels with detected hits and on-track" << endl;
        cout << "  hPathLengthHitsOffTrajectory->Draw()  // Only voxels with detected hits but off-track" << endl;
        cout << "  hVoxelMap->Draw(\"BOX\")  // 3D visualization" << endl;
        cout << "  VoxelPathLengths->Scan(\"*\", \"has_hit && !on_trajectory\") // List off-track hits" << endl;
        
        if (successCount > 1) {
            cout << "\nTo view SUMMARY results for all files:" << endl;
            cout << "  root " << outputDir << "/Summary_AllFiles.root" << endl;
            cout << "  // Combined analysis commands:" << endl;
            cout << "  hPathLengthAll->Draw()  // All files combined path lengths" << endl;
            cout << "  hPathLengthWithHitsAll->Draw()  // All hit voxels combined (including off-track)" << endl;
            cout << "  hPathLengthHitsOnTrajectoryAll->Draw()  // All on-track hits combined" << endl;
            cout << "  hVoxelMapAll->Draw(\"BOX\")  // All trajectories and hits combined" << endl;
            cout << "  trajectory_file_0->Draw(); //Trajectory for the first event root" << endl;
            cout << "  trajectory_file_1->Draw(SAME); //Trajectory for the second event root" << endl;
            cout << "  hTotalVoxelsPerFile->Draw()  // Statistics per file" << endl;
            cout << "  SummaryStats->Scan(\"*\")  // Print summary table" << endl;
        }
    }
    
    exit(0);
}