// Modified EventDisplays_3D_Side for batch processing with muon decay detection
// Saves reconstructed hits to TH3 histograms and detects muon decay candidates
#define THIS_NAME EventDisplays_3D_Side_Batch
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
#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPolyLine3D.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TNamed.h"

// Include the framework headers
#include "../../src/tools/global_header.hh"
#include "CalibrationManager.hh"
#include "Hit3DTypes.hh"
#include "Hit3DReconstructor.hh"

using namespace std;

CalibrationManager calibManager;

// New structure to store processed event information
struct ProcessedEvent {
    int eventNumber;
    int spillTag;
    double minTime, maxTime;  // Earliest and latest hit times
    std::vector<Hit3D> hits3D;
    TrajectoryParams trajectory;
    int nHits2D, nHits3D;
    Event* originalEvent;  // Keep reference for later processing
};

// Structure for muon decay candidates
struct MuonDecayCandidate {
    ProcessedEvent event1;  // Muon track
    ProcessedEvent event2;  // Decay electron
    double timeSeparation;  // Time between events (ticks)
    std::vector<std::pair<int, int>> spatialMatches;  // Indices of matching voxels
    double spatialOverlapScore;  // Quality of spatial overlap
};

// Function to create directories
TString createOutputDirectory(const TString& inputFile, const std::string& subdir = "") {
    TString inputDir = gSystem->DirName(inputFile);
    TString outputDir;
    
    if (subdir == "decay") {
        outputDir = inputDir + "/3D_Display_root_MuonDecay";
    } else {
        outputDir = inputDir;
    }
    
    if (gSystem->AccessPathName(outputDir)) {
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            return (subdir == "decay") ? "./3D_Display_root_MuonDecay" : inputDir;
        }
    }
    
    return outputDir;
}

double getCalibratedCharge(Hit* hit, CalibrationManager& calibManager) {
    if (!hit) return 0.0;
    
    std::string viewName;
    int coord1, coord2;

    // Determine view and coordinates
    if (hit->GetView() == 0) {  // XY view
        viewName = "XY";
        coord1 = hit->GetX();
        coord2 = hit->GetY();
    }
    else if (hit->GetView() == 1) {  // XZ view
        viewName = "XZ";
        coord1 = hit->GetX();
        coord2 = hit->GetZ();
    }
    else if (hit->GetView() == 2) {  // ZY view
        viewName = "ZY";
        coord1 = hit->GetZ();
        coord2 = hit->GetY();
    }
    else {
        return hit->GetPE();
    }

    double gain = calibManager.GetGain(viewName, coord1, coord2);
    if (gain > 0) {
        return hit->GetPE() / gain;
    } else {
        return hit->GetPE();
    }
}

// Define the concrete reconstructor class
class ConcreteHit3DReconstructor : public Hit3DReconstructorCRTP<ConcreteHit3DReconstructor> {
public:
    ConcreteHit3DReconstructor(CalibrationManager* calibManager) : calibManager_(calibManager) {}
    
    // This is the method that CRTP calls - framework access here
    std::vector<Hit3D> DoReconstructHits(void* eventPtr) {
        Event* event = static_cast<Event*>(eventPtr);
        
        // Organize hits by view (your existing code from main file)
        ViewHits viewXY, viewXZ, viewZY;
        
        TClonesArray* hits = event->GetHits();
        if (!hits) return std::vector<Hit3D>();
        
        // Your existing hit organization logic here...
        for (int i = 0; i < hits->GetEntries(); i++) {
            Hit* hit = (Hit*)hits->At(i);
            if (!hit) continue;

            double calibratedCharge = getCalibratedCharge(hit, *calibManager_);

            if (hit->GetView() == 0) {  // XY view
            viewXY.indices.push_back(i);
            viewXY.positions.push_back({hit->GetX(), hit->GetY()});
            viewXY.charges.push_back(calibratedCharge);
            viewXY.times.push_back(hit->GetTfromSpill());
            }
            else if (hit->GetView() == 1) {  // XZ view
                viewXZ.indices.push_back(i);
                viewXZ.positions.push_back({hit->GetX(), hit->GetZ()});
                viewXZ.charges.push_back(calibratedCharge);
                viewXZ.times.push_back(hit->GetTfromSpill());
            }
            else if (hit->GetView() == 2) {  // ZY view
                viewZY.indices.push_back(i);
                viewZY.positions.push_back({hit->GetZ(), hit->GetY()});
                viewZY.charges.push_back(calibratedCharge);
                viewZY.times.push_back(hit->GetTfromSpill());
            }
        }
        
        // Use inherited static methods for the algorithms
        std::vector<bool> usedXY(viewXY.indices.size(), false);
        std::vector<bool> usedXZ(viewXZ.indices.size(), false);
        std::vector<bool> usedZY(viewZY.indices.size(), false);
        
        std::vector<Hit3D> hits3D;
        
        // Call static methods from base class
        auto threeViewHits = FindThreeViewMatches(viewXY, viewXZ, viewZY, usedXY, usedXZ, usedZY);
        hits3D.insert(hits3D.end(), threeViewHits.begin(), threeViewHits.end());
        
        auto twoViewHits = FindTwoViewMatches(viewXY, viewXZ, viewZY, usedXY, usedXZ, usedZY);
        hits3D.insert(hits3D.end(), twoViewHits.begin(), twoViewHits.end());
        
        hits3D = ResolveAdjacentHits(hits3D, viewXY, viewXZ, viewZY);
        
        return hits3D;
    }
    
private:
    CalibrationManager* calibManager_;
};

// Function to reconstruct 3D hits from 2D projections
std::vector<Hit3D> reconstructHits(Event* event) {
    ConcreteHit3DReconstructor reconstructor(&calibManager);
    return reconstructor.ReconstructHits(event);  // CRTP magic happens here!
}

// Function to find entry point based on earliest time - now includes resolved partial hits
Hit3D findEntryPoint(const std::vector<Hit3D>& allHits) {
    if (allHits.empty()) {
        Hit3D dummy;
        dummy.x = dummy.y = dummy.z = 0;
        return dummy;
    }
    
    // Filter hits to include both full hits and resolved partial hits (now 3-view)
    std::vector<Hit3D> candidateHits;
    for (const auto& hit : allHits) {
        // Include original full hits and resolved partial hits that now have 3 views
        if ((hit.matchType == "Full" && hit.confidence >= 1.0)) {
            candidateHits.push_back(hit);
        }
    }
    
    if (candidateHits.empty()) {
        // Fallback to any 3-view hit if no full or resolved hits available
        for (const auto& hit : allHits) {
            if (hit.nViews == 3) {
                candidateHits.push_back(hit);
            }
        }
    }
    
    if (candidateHits.empty()) {
        // Last resort: return first hit
        Hit3D dummy = {}; // Use aggregate initialization to zero-initialize all members
        dummy.x = 0.0;
        dummy.y = 0.0;
        dummy.z = 0.0;
        dummy.timeXY = -1.0; // Use -1 to indicate invalid/missing data
        dummy.timeXZ = -1.0;
        dummy.timeZY = -1.0;
        dummy.chargeXY = -1.0;
        dummy.chargeXZ = -1.0;
        dummy.chargeZY = -1.0;
        dummy.matchType = "Dummy"; // Give it a clear name
        dummy.nViews = 0;
        dummy.confidence = 0.0;
        dummy.redistributed = false;
        return dummy;
    }
    
    // std::cout << "Entry point candidates: " << candidateHits.size() 
    //           << " (Full: " << std::count_if(candidateHits.begin(), candidateHits.end(), 
    //                                        [](const Hit3D& h) { return h.matchType == "Full"; })
    //           << ", Resolved: " << std::count_if(candidateHits.begin(), candidateHits.end(), 
    //                                            [](const Hit3D& h) { return h.matchType.find("Resolved") != std::string::npos; })
    //           << ")" << std::endl;
    
    // Primary selection: Find hits closest to x=0 (detector entrance)
    double minXDist = std::numeric_limits<double>::max();
    std::vector<size_t> spatialCandidates;
    
    for (size_t i = 0; i < candidateHits.size(); i++) {
        double xDist = std::abs(candidateHits[i].x - 0.0);
        if (xDist < minXDist - 0.01) {  // New minimum (with small tolerance)
            minXDist = xDist;
            spatialCandidates.clear();
            spatialCandidates.push_back(i);
        } else if (std::abs(xDist - minXDist) < 0.01) {  // Same minimum (within tolerance)
            spatialCandidates.push_back(i);
        }
    }
    
    // If only one candidate closest to x=0, return it
    if (spatialCandidates.size() == 1) {
        // std::cout << "Single hit closest to x=0 found at distance " << minXDist 
        //           << " (type: " << candidateHits[spatialCandidates[0]].matchType << ")" << std::endl;
        return candidateHits[spatialCandidates[0]];
    }
    
    // Multiple candidates with same x distance - use time as tie-breaker
    // std::cout << "Multiple hits at x-distance " << minXDist << ", using earliest time as tie-breaker..." << std::endl;
    
    double minTime = std::numeric_limits<double>::max();
    std::vector<size_t> timeBasedCandidates;
    
    for (size_t idx : spatialCandidates) {
        const Hit3D& hit = candidateHits[idx];
        
        // Collect all valid timestamps for this hit
        std::vector<double> hitTimes;
        if (hit.timeXY >= 0) hitTimes.push_back(hit.timeXY);
        if (hit.timeXZ >= 0) hitTimes.push_back(hit.timeXZ);
        if (hit.timeZY >= 0) hitTimes.push_back(hit.timeZY);
        
        if (hitTimes.empty()) continue;
        
        // Use the minimum time from this hit's available timestamps
        double hitMinTime = *std::min_element(hitTimes.begin(), hitTimes.end());
        
        if (hitMinTime < minTime - 0.01) {  // New minimum (with small tolerance)
            minTime = hitMinTime;
            timeBasedCandidates.clear();
            timeBasedCandidates.push_back(idx);
        } else if (std::abs(hitMinTime - minTime) < 0.01) {  // Same minimum (within tolerance)
            timeBasedCandidates.push_back(idx);
        }
    }
    
    // If only one candidate with earliest time, return it
    if (timeBasedCandidates.size() == 1) {
        // std::cout << "Selected hit with earliest time " << minTime 
        //           << " (type: " << candidateHits[timeBasedCandidates[0]].matchType << ")" << std::endl;
        return candidateHits[timeBasedCandidates[0]];
    }
    
    // Multiple hits with same x distance and same earliest time - average their positions
    double avgX = 0, avgY = 0, avgZ = 0;
    for (size_t idx : timeBasedCandidates) {
        avgX += candidateHits[idx].x;
        avgY += candidateHits[idx].y;
        avgZ += candidateHits[idx].z;
    }
    
    Hit3D entryPoint = candidateHits[timeBasedCandidates[0]];  // Use first hit as template
    entryPoint.x = avgX / timeBasedCandidates.size();
    entryPoint.y = avgY / timeBasedCandidates.size();
    entryPoint.z = avgZ / timeBasedCandidates.size();
    
    // std::cout << "Averaged entry point: (" << entryPoint.x << ", " << entryPoint.y << ", " << entryPoint.z 
    //           << ") from " << timeBasedCandidates.size() << " equivalent hits (x-dist=" << minXDist 
    //           << ", time=" << minTime << ")" << std::endl;
    
    return entryPoint;
}

// Function to perform PCA on hit collection
PCATrajectoryResult performPCA(const std::vector<Hit3D>& hits, double confidenceThreshold = 0.5) {
    PCATrajectoryResult result;
    result.valid = false;
    
    // Filter for high-confidence hits
    std::vector<Hit3D> highConfidenceHits;
    double totalWeight = 0;
    for (const auto& hit : hits) {
        if (hit.confidence >= confidenceThreshold) {
            highConfidenceHits.push_back(hit);
            totalWeight += hit.confidence * hit.charge;
        }
    }
    
    result.nHighConfidenceHits = highConfidenceHits.size();
    
    if (highConfidenceHits.size() < 3) {
        // std::cout << "Not enough high-confidence hits for PCA (need at least 3, have " 
        //           << highConfidenceHits.size() << ")" << std::endl;
        return result;
    }
    
    // Calculate weighted centroid
    double cx = 0, cy = 0, cz = 0;
    for (const auto& hit : highConfidenceHits) {
        double weight = hit.confidence * hit.charge;
        cx += hit.x * weight;
        cy += hit.y * weight;
        cz += hit.z * weight;
    }
    cx /= totalWeight;
    cy /= totalWeight;
    cz /= totalWeight;
    
    result.centroid_x = cx;
    result.centroid_y = cy;
    result.centroid_z = cz;
    
    // Build covariance matrix with weighted hits
    TMatrixD covMatrix(3, 3);
    covMatrix.Zero();
    
    for (const auto& hit : highConfidenceHits) {
        double weight = hit.confidence * hit.charge / totalWeight;
        double dx = hit.x - cx;
        double dy = hit.y - cy;
        double dz = hit.z - cz;
        
        covMatrix(0,0) += weight * dx * dx;
        covMatrix(0,1) += weight * dx * dy;
        covMatrix(0,2) += weight * dx * dz;
        covMatrix(1,0) += weight * dy * dx;
        covMatrix(1,1) += weight * dy * dy;
        covMatrix(1,2) += weight * dy * dz;
        covMatrix(2,0) += weight * dz * dx;
        covMatrix(2,1) += weight * dz * dy;
        covMatrix(2,2) += weight * dz * dz;
    }
    
    // Perform eigenvalue decomposition
    TVectorD eigenValues;
    TMatrixD eigenVectors = covMatrix.EigenVectors(eigenValues);
    
    // Sort eigenvalues and eigenvectors in descending order
    std::vector<std::pair<double, int>> eigenPairs;
    for (int i = 0; i < 3; i++) {
        eigenPairs.push_back({eigenValues[i], i});
    }
    std::sort(eigenPairs.begin(), eigenPairs.end(), 
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Store sorted eigenvalues and eigenvectors
    for (int i = 0; i < 3; i++) {
        result.eigenvalues[i] = eigenPairs[i].first;
        int idx = eigenPairs[i].second;
        
        if (i == 0) {  // First principal component
            result.pc1_x = eigenVectors(0, idx);
            result.pc1_y = eigenVectors(1, idx);
            result.pc1_z = eigenVectors(2, idx);
        } else if (i == 1) {  // Second principal component
            result.pc2_x = eigenVectors(0, idx);
            result.pc2_y = eigenVectors(1, idx);
            result.pc2_z = eigenVectors(2, idx);
        } else {  // Third principal component
            result.pc3_x = eigenVectors(0, idx);
            result.pc3_y = eigenVectors(1, idx);
            result.pc3_z = eigenVectors(2, idx);
        }
    }
    
    // Calculate explained variance
    result.totalVariance = result.eigenvalues[0] + result.eigenvalues[1] + result.eigenvalues[2];
    if (result.totalVariance > 0) {
        for (int i = 0; i < 3; i++) {
            result.explainedVariance[i] = result.eigenvalues[i] / result.totalVariance;
            result.cumulativeVariance[i] = (i == 0) ? result.explainedVariance[i] : 
                                          result.cumulativeVariance[i-1] + result.explainedVariance[i];
        }
    }
    
    // Calculate track quality metrics
    result.linearity = (result.eigenvalues[0] > 0) ? 
                      (result.eigenvalues[0] - result.eigenvalues[1]) / result.eigenvalues[0] : 0;
    
    result.planarity = (result.eigenvalues[0] > 0) ? 
                      (result.eigenvalues[1] - result.eigenvalues[2]) / result.eigenvalues[0] : 0;
    
    result.trackWidth = (result.eigenvalues[0] > 0) ? 
                       sqrt(result.eigenvalues[1] / result.eigenvalues[0]) : 0;
    
    result.trackThickness = (result.eigenvalues[0] > 0) ? 
                           sqrt(result.eigenvalues[2] / result.eigenvalues[0]) : 0;
    
    // Find entry point along the principal axis
    Hit3D entryPoint = findEntryPoint(highConfidenceHits);
    
    // Project entry point onto principal axis
    double t_entry = (entryPoint.x - cx) * result.pc1_x + 
                     (entryPoint.y - cy) * result.pc1_y + 
                     (entryPoint.z - cz) * result.pc1_z;
    
    // Set trajectory parameters
    result.x0 = cx + t_entry * result.pc1_x;
    result.y0 = cy + t_entry * result.pc1_y;
    result.z0 = cz + t_entry * result.pc1_z;
    
    // Direction vector is the first principal component
    // Ensure it points in positive x direction
    if (result.pc1_x < 0) {
        result.vx = -result.pc1_x;
        result.vy = -result.pc1_y;
        result.vz = -result.pc1_z;
    } else {
        result.vx = result.pc1_x;
        result.vy = result.pc1_y;
        result.vz = result.pc1_z;
    }
    
    // Calculate angles (matching original convention)
    // Theta: angle from +X axis (beam direction)
    result.theta = acos(result.vx) * 180.0 / TMath::Pi();
    
    // Phi: angle from -Z axis in XZ plane, rotated by 90 degrees
    // This matches the original calculation
    double phi_from_z = atan2(result.vz, result.vx);
    result.phi = (phi_from_z + TMath::Pi()/2.0) * 180.0 / TMath::Pi();
    
    // Calculate residuals and outlier detection
    result.rmsResidual = 0;
    result.maxResidual = 0;
    result.nOutliers = 0;
    result.outlierThreshold = 2.5 * sqrt(result.eigenvalues[1] + result.eigenvalues[2]);  // 2 sigma in transverse plane
    
    double sumWeightedDist2 = 0;
    double sumLongitudinal = 0, sumTransverse = 0, sumNormal = 0;
    int nValidHits = 0;
    
    for (const auto& hit : highConfidenceHits) {
        // Project hit onto principal axes
        double dx = hit.x - cx;
        double dy = hit.y - cy;
        double dz = hit.z - cz;
        
        double t_long = dx * result.pc1_x + dy * result.pc1_y + dz * result.pc1_z;
        double t_trans = dx * result.pc2_x + dy * result.pc2_y + dz * result.pc2_z;
        double t_norm = dx * result.pc3_x + dy * result.pc3_y + dz * result.pc3_z;
        
        // Calculate residual distance from principal axis
        double residual = sqrt(t_trans * t_trans + t_norm * t_norm);
        double weight = hit.confidence * hit.charge / totalWeight;
        
        sumWeightedDist2 += weight * residual * residual;
        result.maxResidual = std::max(result.maxResidual, residual);
        
        if (residual > result.outlierThreshold) {
            result.nOutliers++;
        }
        
        sumLongitudinal += std::abs(t_long);
        sumTransverse += std::abs(t_trans);
        sumNormal += std::abs(t_norm);
        nValidHits++;
    }
    
    result.rmsResidual = sqrt(sumWeightedDist2);
    result.meanLongitudinalDistance = sumLongitudinal / nValidHits;
    result.meanTransverseDistance = sumTransverse / nValidHits;
    result.meanNormalDistance = sumNormal / nValidHits;
    
    // Calculate weighted distance from centroid to entry point
    double centroidDist = sqrt((result.x0 - cx) * (result.x0 - cx) +
                              (result.y0 - cy) * (result.y0 - cy) +
                              (result.z0 - cz) * (result.z0 - cz));
    result.weightedCentroidDistance = centroidDist;
    
    // Calculate confidence score for single track
    // High linearity, low planarity, low RMS residual indicate good single track
    double linearityScore = result.linearity;  // 0 to 1, higher is better
    double planarityPenalty = result.planarity;  // 0 to 1, lower is better
    double residualAccpetancy = 0.5;
    double residualPenalty = exp(-0.5 * pow((result.rmsResidual)/0.70, 2));  // approaches 1 for small residuals
    double outlierRatio = (double)result.nOutliers / (double)result.nHighConfidenceHits;
    double outlierPenalty = exp(-outlierRatio / 0.3);;  // exponential penalty for outliers
    
    result.confidenceScore = linearityScore * residualPenalty * (1.0 - 0.5 * planarityPenalty);
    
    result.valid = true;
    
    // Print PCA results
    // std::cout << "\n=== PCA Analysis Results ===" << std::endl;
    // std::cout << "Number of high-confidence hits: " << result.nHighConfidenceHits << std::endl;
    // std::cout << "Centroid: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;
    // std::cout << "Eigenvalues: " << result.eigenvalues[0] << ", " 
    //           << result.eigenvalues[1] << ", " << result.eigenvalues[2] << std::endl;
    // std::cout << "Explained variance: " << result.explainedVariance[0] * 100 << "%, "
    //           << result.explainedVariance[1] * 100 << "%, "
    //           << result.explainedVariance[2] * 100 << "%" << std::endl;
    // std::cout << "Linearity: " << result.linearity << std::endl;
    // std::cout << "Planarity: " << result.planarity << std::endl;
    // std::cout << "Residual Penalty: " << residualPenalty << std::endl;
    // std::cout << "Outlier Penalty: " << outlierPenalty << std::endl;
    // std::cout << "Planarity Penalty: " << planarityPenalty << std::endl;
    // std::cout << "Track width: " << result.trackWidth << " cm" << std::endl;
    // std::cout << "Track thickness: " << result.trackThickness << " cm" << std::endl;
    // std::cout << "RMS residual: " << result.rmsResidual << " cm" << std::endl;
    // std::cout << "Max residual: " << result.maxResidual << " cm" << std::endl;
    // std::cout << "Number of outliers: " << result.nOutliers << std::endl;
    // std::cout << "Single-track confidence: " << result.confidenceScore << std::endl;
    
    return result;
}

// Simplified function - check if PCA is viable
bool canPerformPCA(const std::vector<Hit3D>& hits, double confidenceThreshold = 0.5) {
    int highConfidenceCount = 0;
    for (const auto& hit : hits) {
        if (hit.confidence >= confidenceThreshold) {
            highConfidenceCount++;
        }
    }
    return highConfidenceCount >= 3;
}

// Extract spill tag from event hits
int getSpillTag(Event* event) {
    TClonesArray* hits = event->GetHits();
    if (!hits || hits->GetEntries() == 0) return -1;
    
    Hit* firstHit = (Hit*)hits->At(0);
    if (!firstHit) return -1;
    
    return firstHit->GetSpillTag();  // Get spill tag from first hit
}

// Get time range of event (min and max hit times)
std::pair<double, double> getEventTimeRange(const std::vector<Hit3D>& hits3D) {
    if (hits3D.empty()) return {-1, -1};
    
    double minTime = std::numeric_limits<double>::max();
    double maxTime = std::numeric_limits<double>::lowest();
    
    for (const auto& hit : hits3D) {
        // Check all available times
        if (hit.timeXY >= 0) {
            minTime = std::min(minTime, hit.timeXY);
            maxTime = std::max(maxTime, hit.timeXY);
        }
        if (hit.timeXZ >= 0) {
            minTime = std::min(minTime, hit.timeXZ);
            maxTime = std::max(maxTime, hit.timeXZ);
        }
        if (hit.timeZY >= 0) {
            minTime = std::min(minTime, hit.timeZY);
            maxTime = std::max(maxTime, hit.timeZY);
        }
    }
    
    return {minTime, maxTime};
}

// Check temporal matching for muon decay
bool checkTemporalMatching(const ProcessedEvent& event1, const ProcessedEvent& event2) {
    double timeSeparation = event2.minTime - event1.maxTime;
    return (timeSeparation >= 400.0 && timeSeparation <= 1200.0);
}

// Check spatial matching for muon decay
std::vector<std::pair<int, int>> checkSpatialMatching(const ProcessedEvent& event1, const ProcessedEvent& event2, double confidenceThreshold = 0.5) {
    std::vector<std::pair<int, int>> matches;
    
    // Get high-confidence hits from both events
    std::vector<int> highConf1, highConf2;
    for (size_t i = 0; i < event1.hits3D.size(); i++) {
        if (event1.hits3D[i].confidence >= confidenceThreshold) {
            highConf1.push_back(i);
        }
    }
    for (size_t i = 0; i < event2.hits3D.size(); i++) {
        if (event2.hits3D[i].confidence >= confidenceThreshold) {
            highConf2.push_back(i);
        }
    }
    
    // Check for spatial overlap
    for (int idx1 : highConf1) {
        for (int idx2 : highConf2) {
            const Hit3D& hit1 = event1.hits3D[idx1];
            const Hit3D& hit2 = event2.hits3D[idx2];
            
            // Check if coordinates match
            if (std::abs(hit1.x - hit2.x) < 1.0 && 
                std::abs(hit1.y - hit2.y) < 1.0 && 
                std::abs(hit1.z - hit2.z) < 1.0) {
                matches.push_back({idx1, idx2});
            }
        }
    }
    
    return matches;
}

// Simplified trajectory fitting
TrajectoryParams fitTrajectory(const std::vector<Hit3D>& allHits) {
    TrajectoryParams params;
    params.valid = false;
    
    // Simple PCA implementation
    if (allHits.size() < 3) return params;
    
    // Calculate centroid
    double cx = 0, cy = 0, cz = 0;
    for (const auto& hit : allHits) {
        cx += hit.x;
        cy += hit.y;
        cz += hit.z;
    }
    cx /= allHits.size();
    cy /= allHits.size();
    cz /= allHits.size();
    
    // Simple direction calculation (can be enhanced with full PCA)
    params.x0 = cx;
    params.y0 = cy;
    params.z0 = cz;
    params.vx = 1.0;  // Simplified
    params.vy = 0.0;
    params.vz = 0.0;
    params.theta = 0.0;
    params.phi = 0.0;
    params.valid = true;
    
    return params;
}

// Process single event and store in memory
ProcessedEvent* ProcessEventForDecayDetection(Event* event, int eventNumber) {
    // Filter 1: Skip events with < 5 hits
    if (event->GetNHits() < 5) {
        return nullptr;
    }
    
    // Reconstruct 3D hits
    std::vector<Hit3D> hits3D = reconstructHits(event);
    
    // Filter 2: Check if we can perform PCA
    if (!canPerformPCA(hits3D)) {
        return nullptr;
    }
    
    // Fit trajectory
    TrajectoryParams trajectory = fitTrajectory(hits3D);
    
    if (!trajectory.valid) {
        return nullptr;
    }
    
    // Get spill tag and time range
    int spillTag = getSpillTag(event);
    auto timeRange = getEventTimeRange(hits3D);
    
    // Create processed event
    ProcessedEvent* processedEvent = new ProcessedEvent();
    processedEvent->eventNumber = eventNumber;
    processedEvent->spillTag = spillTag;
    processedEvent->minTime = timeRange.first;
    processedEvent->maxTime = timeRange.second;
    processedEvent->hits3D = hits3D;
    processedEvent->trajectory = trajectory;
    processedEvent->nHits2D = event->GetNHits();
    processedEvent->nHits3D = hits3D.size();
    processedEvent->originalEvent = event;  // Keep reference
    
    return processedEvent;
}

// Save individual event to file
void SaveIndividualEvent(const ProcessedEvent& processedEvent, const TString& inputFile) {
    TString outputDir = createOutputDirectory(inputFile, "");
    TString inputBaseName = gSystem->BaseName(inputFile);
    inputBaseName.ReplaceAll(".root", "");
    
    TString outputFile = Form("%s/%s_3D_event%d_PCA.root", 
                             outputDir.Data(), inputBaseName.Data(), processedEvent.eventNumber);
    
    TFile* fOutput = new TFile(outputFile, "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
        return;
    }
    
    // Create TH3 histograms
    TH3F* h3D_All = new TH3F("h3D_All", 
                         Form("Event %d - All 3D Matches;Z [cm];Y [cm];X [cm]", processedEvent.eventNumber),
                         48, 0, 48, 8, 0, 8, 24, 0, 24);
    
    // Fill histograms
    for (const auto& hit : processedEvent.hits3D) {
        int binZ = h3D_All->GetXaxis()->FindBin(hit.z);
        int binY = h3D_All->GetYaxis()->FindBin(hit.y);
        int binX = h3D_All->GetZaxis()->FindBin(hit.x);
        h3D_All->SetBinContent(binZ, binY, binX, hit.charge);
    }
    
    h3D_All->Write();
    
    // Save trajectory tree
    TTree* trajTree = new TTree("Trajectory", "Basic Trajectory Parameters");
    
    Double_t tr_x0 = processedEvent.trajectory.x0;
    Double_t tr_y0 = processedEvent.trajectory.y0;
    Double_t tr_z0 = processedEvent.trajectory.z0;
    Double_t tr_vx = processedEvent.trajectory.vx;
    Double_t tr_vy = processedEvent.trajectory.vy;
    Double_t tr_vz = processedEvent.trajectory.vz;
    Int_t tr_spillTag = processedEvent.spillTag;
    Double_t tr_minTime = processedEvent.minTime;
    Double_t tr_maxTime = processedEvent.maxTime;
    
    trajTree->Branch("x0", &tr_x0, "x0/D");
    trajTree->Branch("y0", &tr_y0, "y0/D");
    trajTree->Branch("z0", &tr_z0, "z0/D");
    trajTree->Branch("vx", &tr_vx, "vx/D");
    trajTree->Branch("vy", &tr_vy, "vy/D");
    trajTree->Branch("vz", &tr_vz, "vz/D");
    trajTree->Branch("spillTag", &tr_spillTag, "spillTag/I");
    trajTree->Branch("minTime", &tr_minTime, "minTime/D");
    trajTree->Branch("maxTime", &tr_maxTime, "maxTime/D");
    
    trajTree->Fill();
    trajTree->Write();
    
    fOutput->Close();
    delete fOutput;
}

// Save muon decay candidate to merged file
void SaveMuonDecayCandidate(const MuonDecayCandidate& candidate, const TString& inputFile, int candidateNumber) {
    TString outputDir = createOutputDirectory(inputFile, "decay");
    TString inputBaseName = gSystem->BaseName(inputFile);
    inputBaseName.ReplaceAll(".root", "");
    
    TString outputFile = Form("%s/%s_MuonDecay_candidate%d_events%d_%d.root", 
                             outputDir.Data(), inputBaseName.Data(), candidateNumber,
                             candidate.event1.eventNumber, candidate.event2.eventNumber);
    
    TFile* fOutput = new TFile(outputFile, "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Cannot create decay candidate file " << outputFile << std::endl;
        return;
    }
    
    // Create combined TH3 histograms
    TH3F* h3D_Combined = new TH3F("h3D_Combined", 
                              Form("Muon Decay Candidate %d - Events %d + %d;Z [cm];Y [cm];X [cm]", 
                                   candidateNumber, candidate.event1.eventNumber, candidate.event2.eventNumber),
                              48, 0, 48, 8, 0, 8, 24, 0, 24);
    
    TH3F* h3D_Muon = new TH3F("h3D_Muon", 
                          Form("Muon Track - Event %d;Z [cm];Y [cm];X [cm]", candidate.event1.eventNumber),
                          48, 0, 48, 8, 0, 8, 24, 0, 24);
    
    TH3F* h3D_Decay = new TH3F("h3D_Decay", 
                           Form("Decay Electron - Event %d;Z [cm];Y [cm];X [cm]", candidate.event2.eventNumber),
                           48, 0, 48, 8, 0, 8, 24, 0, 24);
    
    // Fill histograms
    // Muon track (Event 1)
    for (const auto& hit : candidate.event1.hits3D) {
        int binZ = h3D_Combined->GetXaxis()->FindBin(hit.z);
        int binY = h3D_Combined->GetYaxis()->FindBin(hit.y);
        int binX = h3D_Combined->GetZaxis()->FindBin(hit.x);
        
        h3D_Combined->SetBinContent(binZ, binY, binX, hit.charge);
        h3D_Muon->SetBinContent(binZ, binY, binX, hit.charge);
    }
    
    // Decay electron (Event 2) - add to combined, don't overwrite
    for (const auto& hit : candidate.event2.hits3D) {
        int binZ = h3D_Combined->GetXaxis()->FindBin(hit.z);
        int binY = h3D_Combined->GetYaxis()->FindBin(hit.y);
        int binX = h3D_Combined->GetZaxis()->FindBin(hit.x);
        
        double existingCharge = h3D_Combined->GetBinContent(binZ, binY, binX);
        h3D_Combined->SetBinContent(binZ, binY, binX, existingCharge + hit.charge);
        h3D_Decay->SetBinContent(binZ, binY, binX, hit.charge);
    }
    
    h3D_Combined->Write();
    h3D_Muon->Write();
    h3D_Decay->Write();
    
    // Save combined hits tree with event labeling
    TTree* combinedTree = new TTree("CombinedHits", "Combined Hits from Muon Decay Events");
    
    Double_t t_x, t_y, t_z, t_charge, t_time, t_confidence;
    Int_t t_eventNumber, t_isDecayElectron;
    Char_t t_matchType[100];
    
    combinedTree->Branch("x", &t_x, "x/D");
    combinedTree->Branch("y", &t_y, "y/D");
    combinedTree->Branch("z", &t_z, "z/D");
    combinedTree->Branch("charge", &t_charge, "charge/D");
    combinedTree->Branch("time", &t_time, "time/D");
    combinedTree->Branch("confidence", &t_confidence, "confidence/D");
    combinedTree->Branch("eventNumber", &t_eventNumber, "eventNumber/I");
    combinedTree->Branch("isDecayElectron", &t_isDecayElectron, "isDecayElectron/I");
    combinedTree->Branch("matchType", t_matchType, "matchType/C");
    
    // Add muon hits (Event 1)
    for (const auto& hit : candidate.event1.hits3D) {
        t_x = hit.x;
        t_y = hit.y;
        t_z = hit.z;
        t_charge = hit.charge;
        t_time = std::min({hit.timeXY >= 0 ? hit.timeXY : 1e9, 
                          hit.timeXZ >= 0 ? hit.timeXZ : 1e9, 
                          hit.timeZY >= 0 ? hit.timeZY : 1e9});
        t_confidence = hit.confidence;
        t_eventNumber = candidate.event1.eventNumber;
        t_isDecayElectron = 0;  // Muon track
        strcpy(t_matchType, hit.matchType.c_str());
        
        combinedTree->Fill();
    }
    
    // Add decay electron hits (Event 2)
    for (const auto& hit : candidate.event2.hits3D) {
        t_x = hit.x;
        t_y = hit.y;
        t_z = hit.z;
        t_charge = hit.charge;
        t_time = std::min({hit.timeXY >= 0 ? hit.timeXY : 1e9, 
                          hit.timeXZ >= 0 ? hit.timeXZ : 1e9, 
                          hit.timeZY >= 0 ? hit.timeZY : 1e9});
        t_confidence = hit.confidence;
        t_eventNumber = candidate.event2.eventNumber;
        t_isDecayElectron = 1;  // Decay electron
        strcpy(t_matchType, hit.matchType.c_str());
        
        combinedTree->Fill();
    }
    
    combinedTree->Write();
    
    // Save decay information
    TTree* decayTree = new TTree("DecayInfo", "Muon Decay Information");
    
    Int_t muonEvent = candidate.event1.eventNumber;
    Int_t decayEvent = candidate.event2.eventNumber;
    Int_t spillTag = candidate.event1.spillTag;
    Double_t timeSep = candidate.timeSeparation;
    Int_t nSpatialMatches = candidate.spatialMatches.size();
    Double_t overlapScore = candidate.spatialOverlapScore;
    Double_t muonMinTime = candidate.event1.minTime;
    Double_t muonMaxTime = candidate.event1.maxTime;
    Double_t decayMinTime = candidate.event2.minTime;
    Double_t decayMaxTime = candidate.event2.maxTime;
    
    decayTree->Branch("muonEventNumber", &muonEvent, "muonEventNumber/I");
    decayTree->Branch("decayEventNumber", &decayEvent, "decayEventNumber/I");
    decayTree->Branch("spillTag", &spillTag, "spillTag/I");
    decayTree->Branch("timeSeparation", &timeSep, "timeSeparation/D");
    decayTree->Branch("nSpatialMatches", &nSpatialMatches, "nSpatialMatches/I");
    decayTree->Branch("spatialOverlapScore", &overlapScore, "spatialOverlapScore/D");
    decayTree->Branch("muonMinTime", &muonMinTime, "muonMinTime/D");
    decayTree->Branch("muonMaxTime", &muonMaxTime, "muonMaxTime/D");
    decayTree->Branch("decayMinTime", &decayMinTime, "decayMinTime/D");
    decayTree->Branch("decayMaxTime", &decayMaxTime, "decayMaxTime/D");
    
    decayTree->Fill();
    decayTree->Write();
    
    // Save summary information
    TNamed decayInfo("DecayInfo", Form("Muon Decay Candidate: Events %d->%d, ΔT=%.1f ticks, %d spatial matches", 
                                      candidate.event1.eventNumber, candidate.event2.eventNumber,
                                      candidate.timeSeparation, (int)candidate.spatialMatches.size()));
    decayInfo.Write();
    
    fOutput->Close();
    delete fOutput;
}

// Main muon decay detection function
std::vector<MuonDecayCandidate> detectMuonDecays(const std::vector<ProcessedEvent*>& processedEvents) {
    std::vector<MuonDecayCandidate> candidates;
    
    // Group events by spill
    std::map<int, std::vector<ProcessedEvent*>> spillGroups;
    for (auto* event : processedEvents) {
        if (event && event->spillTag >= 0) {
            spillGroups[event->spillTag].push_back(event);
        }
    }
    
    std::cout << "\n=== Muon Decay Detection ===" << std::endl;
    std::cout << "Total processed events: " << processedEvents.size() << std::endl;
    std::cout << "Number of spills: " << spillGroups.size() << std::endl;
    
    int totalPairsChecked = 0;
    int temporalMatches = 0;
    int spatialMatches = 0;
    
    // Check each spill for muon decay candidates
    for (auto& spillPair : spillGroups) {
        int spillTag = spillPair.first;
        auto& eventsInSpill = spillPair.second;
        
        if (eventsInSpill.size() < 2) continue;  // Need at least 2 events
        
        // Sort events by minimum time within spill
        std::sort(eventsInSpill.begin(), eventsInSpill.end(), 
                 [](const ProcessedEvent* a, const ProcessedEvent* b) {
                     return a->minTime < b->minTime;
                 });
        
        std::cout << "Spill " << spillTag << ": " << eventsInSpill.size() << " events" << std::endl;
        
        // Check all pairs within this spill
        for (size_t i = 0; i < eventsInSpill.size(); i++) {
            for (size_t j = i + 1; j < eventsInSpill.size(); j++) {
                ProcessedEvent* event1 = eventsInSpill[i];  // Earlier event (potential muon)
                ProcessedEvent* event2 = eventsInSpill[j];  // Later event (potential decay electron)
                
                totalPairsChecked++;
                
                // Check temporal matching first (faster)
                if (!checkTemporalMatching(*event1, *event2)) {
                    continue;
                }
                temporalMatches++;
                
                // Check spatial matching
                auto spatialMatchList = checkSpatialMatching(*event1, *event2, 0.5);
                if (spatialMatchList.empty()) {
                    continue;
                }
                spatialMatches++;
                
                // Create muon decay candidate
                MuonDecayCandidate candidate;
                candidate.event1 = *event1;
                candidate.event2 = *event2;
                candidate.timeSeparation = event2->minTime - event1->maxTime;
                candidate.spatialMatches = spatialMatchList;
                candidate.spatialOverlapScore = (double)spatialMatchList.size() / 
                                              std::max(event1->nHits3D, event2->nHits3D);
                
                candidates.push_back(candidate);
                
                std::cout << "  Found decay candidate: Events " << event1->eventNumber 
                          << " -> " << event2->eventNumber 
                          << " (ΔT=" << candidate.timeSeparation << " ticks, " 
                          << spatialMatchList.size() << " spatial matches)" << std::endl;
            }
        }
    }
    
    std::cout << "\n=== Decay Detection Summary ===" << std::endl;
    std::cout << "Total event pairs checked: " << totalPairsChecked << std::endl;
    std::cout << "Temporal matches (400-1200 ticks): " << temporalMatches << std::endl;
    std::cout << "Spatial matches (overlapping voxels): " << spatialMatches << std::endl;
    std::cout << "Muon decay candidates found: " << candidates.size() << std::endl;
    
    return candidates;
}

// Modified main entry point function for batch processing with muon decay detection
void EventDisplays_3D_Side_Batch() {
    gROOT->SetBatch(kTRUE);
    
    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    calibManager.LoadCalibration("../../src/tools/calibration_files/gain_calibration_strict_68813_time3.txt");
    // Check arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_events.root>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /path/to/events.root" << std::endl;
        std::cerr << "This will process events and detect muon decay candidates" << std::endl;
        return;
    }
    
    TString inputFile = argv[1];
    
    std::cout << "=== Batch Processing Mode with Muon Decay Detection ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Filters:" << std::endl;
    std::cout << "  - Minimum 2D hits per event: 5" << std::endl;
    std::cout << "  - Minimum high-confidence 3D hits for PCA: 3" << std::endl;
    std::cout << "  - Valid trajectory fitting required" << std::endl;
    std::cout << "Muon decay criteria:" << std::endl;
    std::cout << "  - Temporal: 400-1200 ticks between event end and start" << std::endl;
    std::cout << "  - Spatial: At least one overlapping high-confidence voxel" << std::endl;
    
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
    
    Long64_t nEvents = tree->GetEntries();
    std::cout << "Total events in file: " << nEvents << std::endl;
    
    // Set up event reading
    Event* event = new Event();
    tree->SetBranchAddress("Event", &event);
    
    // Statistics tracking
    int processedEvents = 0;
    int skippedEvents = 0;
    std::vector<ProcessedEvent*> allProcessedEvents;
    
    auto startTime = std::chrono::steady_clock::now();
    
    std::cout << "\n=== Phase 1: Processing Individual Events ===" << std::endl;
    
    // Phase 1: Process all events and store in memory
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
        tree->GetEntry(iEvent);
        
        // Progress reporting
        if (iEvent % 500 == 0 || iEvent == nEvents - 1) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double progress = 100.0 * iEvent / nEvents;
            double rate = iEvent > 0 ? (double)iEvent / elapsed : 0;
            
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% "
                      << "(" << iEvent << "/" << nEvents << ") "
                      << "Rate: " << std::setprecision(1) << rate << " events/sec "
                      << "Processed: " << processedEvents << " Skipped: " << skippedEvents;
            std::cout.flush();
        }
        
        // Process this event
        ProcessedEvent* processedEvent = ProcessEventForDecayDetection(event, iEvent);
        
        if (processedEvent) {
            allProcessedEvents.push_back(processedEvent);
            processedEvents++;
            
            // // Save individual event file
            // SaveIndividualEvent(*processedEvent, inputFile);
        } else {
            skippedEvents++;
        }
    }
    
    auto phase1EndTime = std::chrono::steady_clock::now();
    auto phase1Time = std::chrono::duration_cast<std::chrono::seconds>(phase1EndTime - startTime).count();
    
    std::cout << "\n\n=== Phase 1 Complete ===" << std::endl;
    std::cout << "Processed events: " << processedEvents << std::endl;
    std::cout << "Skipped events: " << skippedEvents << std::endl;
    std::cout << "Processing time: " << phase1Time << " seconds" << std::endl;
    
    // Phase 2: Detect muon decay candidates
    std::cout << "\n=== Phase 2: Detecting Muon Decay Candidates ===" << std::endl;
    
    std::vector<MuonDecayCandidate> decayCandidates = detectMuonDecays(allProcessedEvents);
    
    // Phase 3: Save muon decay candidates
    std::cout << "\n=== Phase 3: Saving Muon Decay Candidates ===" << std::endl;
    
    for (size_t i = 0; i < decayCandidates.size(); i++) {
        SaveMuonDecayCandidate(decayCandidates[i], inputFile, i);
        std::cout << "Saved decay candidate " << i << ": Events " 
                  << decayCandidates[i].event1.eventNumber << " -> " 
                  << decayCandidates[i].event2.eventNumber << std::endl;
    }
    
    auto endTime = std::chrono::steady_clock::now();
    auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    
    std::cout << "\n=== Processing Complete ===" << std::endl;
    std::cout << "Total processing time: " << totalTime << " seconds" << std::endl;
    std::cout << "Total events: " << nEvents << std::endl;
    std::cout << "Processed events (PCA-viable): " << processedEvents << std::endl;
    std::cout << "Skipped events: " << skippedEvents << std::endl;
    std::cout << "Muon decay candidates found: " << decayCandidates.size() << std::endl;
    std::cout << "Processing efficiency: " << std::fixed << std::setprecision(1) 
              << (100.0 * processedEvents / nEvents) << "%" << std::endl;
    std::cout << "Decay detection rate: " << std::fixed << std::setprecision(3) 
              << (100.0 * decayCandidates.size() / processedEvents) << "%" << std::endl;
    
    // Create summary file
    TString summaryFile = inputFile;
    summaryFile.ReplaceAll(".root", "_MuonDecay_summary.txt");
    std::ofstream summary(summaryFile.Data());
    summary << "Muon Decay Detection Summary\n";
    summary << "===========================\n";
    summary << "Input file: " << inputFile << "\n";
    summary << "Total events: " << nEvents << "\n";
    summary << "Processed events (PCA-viable): " << processedEvents << "\n";
    summary << "Skipped events: " << skippedEvents << "\n";
    summary << "Processing efficiency: " << std::fixed << std::setprecision(1) << (100.0 * processedEvents / nEvents) << "%\n";
    summary << "Muon decay candidates found: " << decayCandidates.size() << "\n";
    summary << "Decay detection rate: " << std::fixed << std::setprecision(3) << (100.0 * decayCandidates.size() / processedEvents) << "%\n";
    summary << "Total processing time: " << totalTime << " seconds\n";
    summary << "\nSelection criteria:\n";
    summary << "- Minimum 2D hits: 5\n";
    summary << "- Minimum high-confidence 3D hits for PCA: 3\n";
    summary << "- Valid trajectory fitting required\n";
    summary << "\nMuon decay criteria:\n";
    summary << "- Temporal separation: 400-1200 ticks\n";
    summary << "- Spatial overlap: ≥1 matching high-confidence voxel\n";
    summary << "\nOutput directories:\n";
    summary << "- Individual events: " << createOutputDirectory(inputFile, "") << "\n";
    summary << "- Muon decay candidates: " << createOutputDirectory(inputFile, "decay") << "\n";
    
    // List all decay candidates
    summary << "\nDecay candidates found:\n";
    for (size_t i = 0; i < decayCandidates.size(); i++) {
        const auto& candidate = decayCandidates[i];
        summary << "Candidate " << i << ": Events " << candidate.event1.eventNumber 
                << " -> " << candidate.event2.eventNumber 
                << " (ΔT=" << candidate.timeSeparation << " ticks, "
                << candidate.spatialMatches.size() << " spatial matches, spill " 
                << candidate.event1.spillTag << ")\n";
    }
    
    summary.close();
    
    std::cout << "Summary written to: " << summaryFile << std::endl;
    std::cout << "Individual events directory: " << createOutputDirectory(inputFile, "") << std::endl;
    std::cout << "Muon decay candidates directory: " << createOutputDirectory(inputFile, "decay") << std::endl;
    std::cout << "To view the 3D histograms:" << std::endl;
    std::cout << "TCanvas *c1 = new TCanvas(\"c1\", \"Canvas 1\", 800, 600);"<< std::endl;
    std::cout << "TCanvas *c2 = new TCanvas(\"c2\", \"Canvas 2\", 800, 600);"<< std::endl;
    std::cout << "TCanvas *c3 = new TCanvas(\"c3\", \"Canvas 3\", 800, 600);"<< std::endl;
    std::cout << "c1->cd()"<< std::endl;
    std::cout << "h3D_Combined->Draw(\"BOX2Z\")"<< std::endl;
    std::cout << "c2->cd()"<< std::endl;
    std::cout << "h3D_Muon->Draw(\"BOX2Z\")"<< std::endl;
    std::cout << "c3->cd()"<< std::endl;
    std::cout << "h3D_Decay->Draw(\"BOX2Z\")"<< std::endl;
    std::cout << ""<< std::endl;
    // Clean up memory
    for (auto* processedEvent : allProcessedEvents) {
        delete processedEvent;
    }
    
    file->Close();
    delete event;
    exit(0);
}