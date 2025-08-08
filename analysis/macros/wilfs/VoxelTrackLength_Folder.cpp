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
    int x, y, z;           // Voxel indices
    double pathLength;     // Path length through this voxel
    bool hasHit;           // Whether this voxel contains a high-confidence hit
    double hitConfidence;  // Confidence of hit in this voxel (0 if no hit)
    double hitCharge;      // Charge of hit in this voxel (0 if no hit)
    double entryT, exitT;  // Parameter values for entry/exit points
    bool onTrajectory;     // Whether this voxel is traversed by the trajectory
    
    VoxelInfo() : x(0), y(0), z(0), pathLength(0), hasHit(false), 
                  hitConfidence(0), hitCharge(0), entryT(0), exitT(0), onTrajectory(false) {}
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
    
    HitInfo() : x(0), y(0), z(0), confidence(0), charge(0) {}
};

// Load trajectory parameters from ROOT file
TrajectoryParams loadTrajectory(TFile* file) {
    TrajectoryParams traj;
    
    TTree* trajTree = (TTree*)file->Get("Trajectory");
    if (!trajTree) {
        cout << "Error: Cannot find Trajectory tree" << endl;
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
    
    cout << "Loaded trajectory: " << endl;
    cout << "  Origin: (" << traj.x0 << ", " << traj.y0 << ", " << traj.z0 << ")" << endl;
    cout << "  Direction: (" << traj.vx << ", " << traj.vy << ", " << traj.vz << ")" << endl;
    cout << "  Valid: " << (traj.valid ? "Yes" : "No") << endl;
    
    return traj;
}

// Load high-confidence hits from ROOT file
vector<HitInfo> loadHighConfidenceHits(TFile* file, double confidenceThreshold = 0.5) {
    vector<HitInfo> hits;
    
    TTree* hitsTree = (TTree*)file->Get("Hits3D");
    if (!hitsTree) {
        cout << "Error: Cannot find Hits3D tree" << endl;
        return hits;
    }
    
    Double_t x, y, z, confidence, charge;
    
    hitsTree->SetBranchAddress("x", &x);
    hitsTree->SetBranchAddress("y", &y);
    hitsTree->SetBranchAddress("z", &z);
    hitsTree->SetBranchAddress("confidence", &confidence);
    hitsTree->SetBranchAddress("charge", &charge);
    
    Long64_t nEntries = hitsTree->GetEntries();
    cout << "Processing " << nEntries << " hits..." << endl;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        hitsTree->GetEntry(i);
        
        if (confidence >= confidenceThreshold) {
            HitInfo hit;
            hit.x = x;
            hit.y = y;
            hit.z = z;
            hit.confidence = confidence;
            hit.charge = charge;
            hits.push_back(hit);
        }
    }
    
    cout << "Found " << hits.size() << " high-confidence hits (threshold >= " 
         << confidenceThreshold << ")" << endl;
    
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
    
    cout << "Trajectory parameter range: t = [" << tMin << ", " << tMax << "]" << endl;
    
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
        cout << "No valid trajectory range found in detector" << endl;
        return voxels;
    }
    
    // Create map for quick hit lookup
    map<tuple<int, int, int>, HitInfo> hitMap;
    for (const auto& hit : hits) {
        int hx = (int)floor(hit.x);
        int hy = (int)floor(hit.y);
        int hz = (int)floor(hit.z);
        hitMap[make_tuple(hx, hy, hz)] = hit;
    }
    
    // Use fine sampling to find all intersected voxels
    set<tuple<int, int, int>> intersectedVoxels;
    int nSamples = 100000;
    double dt = (tMax - tMin) / nSamples;
    
    for (int i = 0; i <= nSamples; i++) {
        double t = tMin + i * dt;
        double x, y, z;
        getTrajectoryPoint(traj, t, x, y, z);
        
        // Check bounds
        if (x >= 0 && x <= 24 && y >= 0 && y <= 8 && z >= 0 && z <= 48) {
            int vx = (int)floor(x);
            int vy = (int)floor(y);
            int vz = (int)floor(z);
            
            // Clamp to valid voxel indices
            vx = max(0, min(23, vx));
            vy = max(0, min(7, vy));
            vz = max(0, min(47, vz));
            
            intersectedVoxels.insert(make_tuple(vx, vy, vz));
        }
    }
    
    cout << "Found " << intersectedVoxels.size() << " intersected voxels" << endl;
    
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
        
        // Check if this voxel has a hit
        auto hitIt = hitMap.find(make_tuple(vx, vy, vz));
        if (hitIt != hitMap.end()) {
            info.hasHit = true;
            info.hitConfidence = hitIt->second.confidence;
            info.hitCharge = hitIt->second.charge;
        }
        
        if (info.pathLength > 1e-9) {  // Only include voxels with non-zero path length
            voxels.push_back(info);
        }
    }
    
    // Now add all hits that are NOT on the trajectory (path length = 0)
    for (const auto& hit : hits) {
        int hx = (int)floor(hit.x);
        int hy = (int)floor(hit.y);
        int hz = (int)floor(hit.z);
        
        // Clamp to valid voxel indices
        hx = max(0, min(23, hx));
        hy = max(0, min(7, hy));
        hz = max(0, min(47, hz));
        
        tuple<int, int, int> voxelKey = make_tuple(hx, hy, hz);
        
        // Check if this voxel is already in our list (i.e., on trajectory)
        bool alreadyIncluded = false;
        for (const auto& existingVoxel : voxels) {
            if (existingVoxel.x == hx && existingVoxel.y == hy && existingVoxel.z == hz) {
                alreadyIncluded = true;
                break;
            }
        }
        
        // If not already included, add it with zero path length
        if (!alreadyIncluded) {
            VoxelInfo info;
            info.x = hx;
            info.y = hy;
            info.z = hz;
            info.pathLength = 0.0;  // Zero path length for hits not on trajectory
            info.hasHit = true;
            info.hitConfidence = hit.confidence;
            info.hitCharge = hit.charge;
            info.onTrajectory = false;
            
            voxels.push_back(info);
        }
    }
    
    cout << "Total voxels in analysis: " << voxels.size() << " (including hits not on trajectory)" << endl;
    
    return voxels;
}

// Main analysis function
void TrajectoryVoxelAnalyzer() {
    gROOT->SetBatch(kTRUE);
    
    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_3D_file.root>" << endl;
        return;
    }
    
    TString inputFile = argv[1];
    cout << "Processing file: " << inputFile << endl;
    
    // Open input file
    TFile* file = TFile::Open(inputFile, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << inputFile << endl;
        return;
    }
    
    // Load trajectory parameters
    TrajectoryParams traj = loadTrajectory(file);
    if (!traj.valid) {
        cerr << "Error: Invalid trajectory found" << endl;
        file->Close();
        return;
    }
    
    // Load high-confidence hits
    vector<HitInfo> hits = loadHighConfidenceHits(file, 0.5);
    
    // Perform voxel traversal
    cout << "\nPerforming voxel traversal..." << endl;
    vector<VoxelInfo> voxels = traverseVoxels(traj, hits);
    
    // Create output file
    TString outputFile = inputFile;
    outputFile.ReplaceAll(".root", "_VoxelPathLengths.root");
    TFile* fOut = new TFile(outputFile, "RECREATE");
    
    // Create histograms
    TH1F* hPathLength = new TH1F("hPathLength", "Path Length Distribution;Path Length [cm];Entries", 
                                100, 0, 2.0);
    TH1F* hPathLengthWithHits = new TH1F("hPathLengthWithHits", 
                                        "Path Length Distribution (Voxels with Hits);Path Length [cm];Entries",
                                        100, 0, 2.0);
    TH1F* hPathLengthNoHits = new TH1F("hPathLengthNoHits", 
                                      "Path Length Distribution (Voxels without Hits);Path Length [cm];Entries",
                                      100, 0, 2.0);
    
    // NEW: Separate histograms for hits on/off trajectory
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
    Double_t path_length, hit_confidence, hit_charge;
    Bool_t has_hit, on_trajectory;
    
    voxelTree->Branch("voxel_x", &voxel_x, "voxel_x/I");
    voxelTree->Branch("voxel_y", &voxel_y, "voxel_y/I");
    voxelTree->Branch("voxel_z", &voxel_z, "voxel_z/I");
    voxelTree->Branch("path_length", &path_length, "path_length/D");
    voxelTree->Branch("has_hit", &has_hit, "has_hit/O");
    voxelTree->Branch("hit_confidence", &hit_confidence, "hit_confidence/D");
    voxelTree->Branch("hit_charge", &hit_charge, "hit_charge/D");
    voxelTree->Branch("on_trajectory", &on_trajectory, "on_trajectory/O");  // NEW branch
    
    // Fill histograms and tree
    double totalPathLength = 0;
    int voxelsWithHits = 0;
    int voxelsWithoutHits = 0;
    int hitsOnTrajectory = 0;
    int hitsOffTrajectory = 0;
    
    for (const auto& voxel : voxels) {
        voxel_x = voxel.x;
        voxel_y = voxel.y;
        voxel_z = voxel.z;
        path_length = voxel.pathLength;
        has_hit = voxel.hasHit;
        hit_confidence = voxel.hitConfidence;
        hit_charge = voxel.hitCharge;
        on_trajectory = voxel.onTrajectory;
        
        voxelTree->Fill();
        
        hPathLength->Fill(path_length);
        if (voxel.hasHit) {
            voxelsWithHits++;
            hPathLengthWithHits->Fill(path_length);  // This now includes zero path length hits
            hPathLengthAllHits->Fill(path_length);   // NEW: All hits histogram
            
            if (voxel.onTrajectory) {
                hitsOnTrajectory++;
                hPathLengthHitsOnTrajectory->Fill(path_length);
            } else {
                hitsOffTrajectory++;
                hPathLengthHitsOffTrajectory->Fill(path_length);  // These will be at zero
            }
        } else {
            hPathLengthNoHits->Fill(path_length);
            voxelsWithoutHits++;
        }
        
        // Fill 3D histogram (note: ROOT uses z,y,x convention)
        // Only fill if path length > 0 to avoid confusion in 3D visualization
        if (voxel.pathLength > 0) {
            hVoxelMap->SetBinContent(voxel.z + 1, voxel.y + 1, voxel.x + 1, path_length);
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
    cout << "\n=== Voxel Path Length Analysis Summary ===" << endl;
    cout << "Total voxels traversed by trajectory: " << (voxels.size() - hitsOffTrajectory) << endl;
    cout << "Total voxels in analysis: " << voxels.size() << endl;
    cout << "Voxels with hits: " << voxelsWithHits << endl;
    cout << "  - Hits on trajectory: " << hitsOnTrajectory << endl;
    cout << "  - Hits off trajectory: " << hitsOffTrajectory << endl;
    cout << "Voxels without hits: " << voxelsWithoutHits << endl;
    cout << "Total path length: " << totalPathLength << " cm" << endl;
    cout << "Average path length per voxel: " << (voxels.size() > 0 ? totalPathLength / voxels.size() : 0) << " cm" << endl;
    
    if (hitsOnTrajectory > 0) {
        double avgPathLengthOnTrajectory = 0;
        for (const auto& voxel : voxels) {
            if (voxel.hasHit && voxel.onTrajectory) avgPathLengthOnTrajectory += voxel.pathLength;
        }
        avgPathLengthOnTrajectory /= hitsOnTrajectory;
        cout << "Average path length (hits on trajectory): " << avgPathLengthOnTrajectory << " cm" << endl;
    }
    
    cout << "\nOutput saved to: " << outputFile << endl;
    cout << "\nTo view results in ROOT:" << endl;
    cout << "  root " << outputFile << endl;
    cout << "  hPathLength->Draw()  // Overall path length distribution" << endl;
    cout << "  hPathLengthWithHits->Draw()  // Path lengths for ALL voxels with hits (including off-track)" << endl;
    cout << "  hPathLengthAllHits->Draw()   // Same as above - all hits including zero path length" << endl;
    cout << "  hPathLengthHitsOnTrajectory->Draw()   // Only hits on trajectory" << endl;
    cout << "  hPathLengthHitsOffTrajectory->Draw()  // Only hits off trajectory (will show peak at zero)" << endl;
    cout << "  hVoxelMap->Draw(\"BOX\")  // 3D visualization of path lengths" << endl;
    cout << "  " << endl;
    cout << "  // To compare on-track vs off-track hits:" << endl;
    cout << "  hPathLengthHitsOnTrajectory->SetLineColor(2)" << endl;
    cout << "  hPathLengthHitsOffTrajectory->SetLineColor(4)" << endl;
    cout << "  hPathLengthHitsOnTrajectory->Draw()" << endl;
    cout << "  hPathLengthHitsOffTrajectory->Draw(\"same\")" << endl;
    cout << "  " << endl;
    cout << "  // To view voxels with trajectory overlay:" << endl;
    cout << "  hVoxelMap->Draw(\"BOX\")" << endl;
    cout << "  TPolyLine3D *traj = (TPolyLine3D*)gDirectory->Get(\"trajectory\")" << endl;
    cout << "  traj->SetLineColor(2)" << endl;
    cout << "  traj->SetLineWidth(3)" << endl;
    cout << "  traj->Draw()" << endl;
    cout << "  " << endl;
    cout << "  // Tree analysis examples:" << endl;
    cout << "  VoxelPathLengths->Scan(\"voxel_x:voxel_y:voxel_z:path_length:hit_confidence\", \"has_hit && !on_trajectory\") // List off-track hits" << endl;
    
    fOut->Close();
    file->Close();
    
    exit(0);
}