// Modified EventDisplays_3D_Side for batch processing with gamma/beta filtering
// Saves reconstructed hits to TH3 histograms with perimeter filtering
#define THIS_NAME EventDisplays_3D_Side_Batch_Filtered
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

struct PerimeterResult {
    bool hasPerimeterVoxels;
    int nPerimeterVoxels;
    std::vector<std::string> perimeterTypes;  // Which boundaries have hits
    
    // New specific flags
    bool hasOnlyX23Perimeter;  // Only X=23 perimeter hits, no other boundaries
    bool hasX23;               // Has any X=23 hits
    bool hasOtherPerimeter;    // Has perimeter hits on other boundaries
    
    // New filtering flags
    bool isGammaEvent;         // All hits z > 23 AND only X=23 perimeter
    bool isBetaEvent;          // All hits z <= 23 AND only X=23 perimeter
    bool passesFilter;         // Event passes one of the two criteria
    
    PerimeterResult() : hasPerimeterVoxels(false), nPerimeterVoxels(0), 
                       hasOnlyX23Perimeter(false), hasX23(false), hasOtherPerimeter(false),
                       isGammaEvent(false), isBetaEvent(false), passesFilter(false) {}
};

// Function to analyze perimeter voxels and apply filtering criteria
PerimeterResult analyzePerimeterVoxels(const std::vector<Hit3D>& hits3D) {
    PerimeterResult result;
    
    // Track which boundaries have hits
    bool hasXMin = false, hasXMax = false;
    bool hasYMin = false, hasYMax = false; 
    bool hasZMin = false, hasZMax = false;
    
    // Track Z-coordinate distribution of X=23 perimeter hits specifically
    bool hasX23HitsLeftOfZ23 = false;   // X=23 hits with z <= 23
    bool hasX23HitsRightOfZ23 = false;  // X=23 hits with z > 23
    
    for (const auto& hit : hits3D) {
        // Convert to bin indices once
        int x = (int)(hit.x - 0.5);
        int y = (int)(hit.y - 0.5);
        int z = (int)(hit.z - 0.5);
        
        // Check if this hit is on any perimeter boundary
        bool isPerimeter = (x == 0 || x == 23 || y == 0 || y == 7 || z == 0 || z == 47);
        
        if (isPerimeter) {
            result.hasPerimeterVoxels = true;
            result.nPerimeterVoxels++;
            
            // Determine which specific boundary(ies) this hit is on
            if (x == 0) hasXMin = true;
            if (x == 23) {
                hasXMax = true;
                // Track Z-coordinate distribution ONLY for X=23 perimeter hits
                if ( z >= 10 && z <= 13) {
                    hasX23HitsLeftOfZ23 = true;
                } else {
                    hasX23HitsRightOfZ23 = true;
                }
            }
            if (y == 0) hasYMin = true;
            if (y == 7) hasYMax = true;
            if (z == 0) hasZMin = true;
            if (z == 47) hasZMax = true;
        }
    }
    
    // Set perimeter flags
    result.hasX23 = hasXMax;
    result.hasOtherPerimeter = hasXMin || hasYMin || hasYMax || hasZMin || hasZMax;
    result.hasOnlyX23Perimeter = hasXMax && !result.hasOtherPerimeter;
    
    // Apply NEW filtering criteria - based on Z-coordinates of X=23 hits only
    // Criterion 1: Only X=23 perimeter hits AND X=23 hits have z <= 23
    result.isBetaEvent = result.hasOnlyX23Perimeter && hasX23HitsLeftOfZ23 && !hasX23HitsRightOfZ23;
    
    // Criterion 2: Only X=23 perimeter hits AND X=23 hits have z > 23  
    result.isGammaEvent = result.hasOnlyX23Perimeter && hasX23HitsRightOfZ23 && !hasX23HitsLeftOfZ23;
    
    // Event passes filter if it satisfies either criterion
    result.passesFilter = result.isGammaEvent || result.isBetaEvent;
    
    // Build list of perimeter types for debugging
    if (hasXMin) result.perimeterTypes.push_back("X=0");
    if (hasXMax) result.perimeterTypes.push_back("X=23");
    if (hasYMin) result.perimeterTypes.push_back("Y=0");
    if (hasYMax) result.perimeterTypes.push_back("Y=7");
    if (hasZMin) result.perimeterTypes.push_back("Z=0");
    if (hasZMax) result.perimeterTypes.push_back("Z=47");
    
    return result;
}

TString createFilteredOutputDirectory(const TString& inputFile, const PerimeterResult& perimeterResult) {
    // Extract directory path from input file
    TString inputDir = gSystem->DirName(inputFile);
    
    // Create output directory based on event type
    TString outputDir;
    if (perimeterResult.isGammaEvent) {
        outputDir = inputDir + "/GammaEvents_OnlyX23_ZRight";
    } else if (perimeterResult.isBetaEvent) {
        outputDir = inputDir + "/BetaEvents_OnlyX23_ZLeft";
    } else {
        // This should not happen if we're filtering correctly
        return "";  // Return empty string to indicate rejection
    }
    
    // Check if directory exists, create if not
    if (gSystem->AccessPathName(outputDir)) {
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            return "";
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
        return candidateHits[spatialCandidates[0]];
    }
    
    // Multiple candidates with same x distance - use time as tie-breaker
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
    double residualPenalty = exp(-0.5 * pow((result.rmsResidual)/0.70, 2));  // approaches 1 for small residuals
    
    result.confidenceScore = linearityScore * residualPenalty * (1.0 - 0.5 * planarityPenalty);
    
    result.valid = true;
    
    return result;
}

// Function to check if event is consistent with single track
bool isSingleTrackEvent(const PCATrajectoryResult& pcaResult, 
                       double linearityThreshold = 0.95,
                       double planarityThreshold = 0.1,
                       double confidenceThreshold = 0.7) {
    
    if (!pcaResult.valid) return false;
    
    // Check multiple criteria
    bool highLinearity = pcaResult.linearity > linearityThreshold;
    bool lowPlanarity = pcaResult.planarity < planarityThreshold;
    bool highConfidence = pcaResult.confidenceScore > confidenceThreshold;
    bool goodExplainedVariance = pcaResult.explainedVariance[0] > 0.9;  // PC1 explains >90% variance
    
    return highLinearity && lowPlanarity && highConfidence && goodExplainedVariance;
}

TrajectoryParams fitTrajectory(const std::vector<Hit3D>& allHits) {
    TrajectoryParams params;
    params.valid = false;
    
    // Perform PCA analysis
    PCATrajectoryResult pcaResult = performPCA(allHits, 0.5);
    
    if (!pcaResult.valid) {
        return params;
    }
    
    // Convert PCA result to TrajectoryParams for compatibility
    params.x0 = pcaResult.x0;
    params.y0 = pcaResult.y0;
    params.z0 = pcaResult.z0;
    params.vx = pcaResult.vx;
    params.vy = pcaResult.vy;
    params.vz = pcaResult.vz;
    params.theta = pcaResult.theta;
    params.phi = pcaResult.phi;
    params.valid = pcaResult.valid;
    
    return params;
}

// Function to print summary
void printSummary(const std::vector<Hit3D>& hits3D, int eventID) {
    int nFull = 0, nPartial = 0, nGroupResolved = 0;
    std::map<std::string, int> typeCounts;
    std::map<double, int> confidenceCounts;
    
    for (const auto& hit : hits3D) {
        if (hit.nViews == 3 && hit.confidence >= 1.0) nFull++;
        else if (hit.matchType.find("Resolved") != std::string::npos) nGroupResolved++;
        else nPartial++;
        
        typeCounts[hit.matchType]++;
        
        // Round confidence to nearest 0.1 for grouping
        double roundedConf = round(hit.confidence / 0.1) * 0.1;
        confidenceCounts[roundedConf]++;
    }
}

// Function to save individual event (like first code)
void SaveFilteredEvent(Event* event, int eventNumber, const std::vector<Hit3D>& hits3D, 
                      const TrajectoryParams& trajectory, const PCATrajectoryResult& pcaResult,
                      const PerimeterResult& perimeterResult, const TString& inputFile) {
    
    // Create filtered output directory
    TString outputDir = createFilteredOutputDirectory(inputFile, perimeterResult);
    if (outputDir.IsNull()) {
        std::cerr << "Error: Could not create output directory for event " << eventNumber << std::endl;
        return;
    }
    
    TString inputBaseName = gSystem->BaseName(inputFile);
    inputBaseName.ReplaceAll(".root", "");
    
    // Create output filename with new naming scheme
    TString outputFile;
    TString eventType = perimeterResult.isGammaEvent ? "Gamma" : "Beta";
    bool isSingleTrack = isSingleTrackEvent(pcaResult);
    TString trackTag = isSingleTrack ? "SingleTrack" : "MultiTrack";

    outputFile = Form("%s/%s_3D_event%d_%s_%s.root", 
                    outputDir.Data(), inputBaseName.Data(), eventNumber, 
                    eventType.Data(), trackTag.Data());
    
    std::cout << "Accepted " << eventType << " event " << eventNumber 
              << ". Saving output to: " << outputFile << std::endl;
    
    // Create output file
    TFile* fOutput = new TFile(outputFile, "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
        return;
    }
    
    // Create TH3 histograms for different hit types
    TH3F* h3D_Full = new TH3F("h3D_Full", 
                          Form("Event %d - Full 3D Matches (3 views);Z [cm];Y [cm];X [cm]", eventNumber),
                          48, 0, 48, 8, 0, 8, 24, 0, 24);

    TH3F* h3D_GroupResolved = new TH3F("h3D_GroupResolved", 
                                     Form("Event %d - Group-Resolved 3D Hits;Z [cm];Y [cm];X [cm]", eventNumber),
                                     48, 0, 48, 8, 0, 8, 24, 0, 24);

    TH3F* h3D_Partial = new TH3F("h3D_Partial", 
                             Form("Event %d - Partial 3D Matches (2 views);Z [cm];Y [cm];X [cm]", eventNumber),
                             48, 0, 48, 8, 0, 8, 24, 0, 24);

    TH3F* h3D_All = new TH3F("h3D_All", 
                         Form("Event %d - All 3D Matches;Z [cm];Y [cm];X [cm]", eventNumber),
                         48, 0, 48, 8, 0, 8, 24, 0, 24);

    TH3F* h3D_Alpha = new TH3F("h3D_Alpha", 
                           Form("Event %d - Hit Confidence Levels;Z [cm];Y [cm];X [cm]", eventNumber),
                           48, 0, 48, 8, 0, 8, 24, 0, 24);

    // Fill histograms with charge as weight and confidence as alpha
    for (const auto& hit : hits3D) {
        int binZ = h3D_All->GetXaxis()->FindBin(hit.z);
        int binY = h3D_All->GetYaxis()->FindBin(hit.y);
        int binX = h3D_All->GetZaxis()->FindBin(hit.x);
        
        h3D_All->SetBinContent(binZ, binY, binX, hit.charge);
        h3D_Alpha->SetBinContent(binZ, binY, binX, hit.confidence);
        
        if (hit.confidence >= 1.0) {
            h3D_Full->SetBinContent(binZ, binY, binX, hit.charge);
        } else if (hit.matchType.find("Resolved") != std::string::npos) {
            h3D_GroupResolved->SetBinContent(binZ, binY, binX, hit.charge);
        } else {
            h3D_Partial->SetBinContent(binZ, binY, binX, hit.charge);
        }
    }

    // Set up color palette for charge visualization
    gStyle->SetPalette(kViridis);

    // Configure histogram display properties
    h3D_All->SetOption("BOX2Z");
    h3D_Full->SetOption("BOX2Z");
    h3D_GroupResolved->SetOption("BOX2Z");
    h3D_Partial->SetOption("BOX2Z");

    h3D_All->SetMinimum(0);
    h3D_Full->SetMinimum(0);
    h3D_GroupResolved->SetMinimum(0);
    h3D_Partial->SetMinimum(0);

    // Create trajectory visualization objects
    TPolyLine3D* trajectoryLine = nullptr;
    TF1* trajX_vs_Z = nullptr;
    TF1* trajY_vs_Z = nullptr;
    TPolyLine3D* thetaZeroLine = nullptr;
    TPolyLine3D* phiZeroLine = nullptr;
    
    // Create reference lines
    thetaZeroLine = new TPolyLine3D(2);
    thetaZeroLine->SetPoint(0, 24, 4, 0);
    thetaZeroLine->SetPoint(1, 24, 4, 24);
    thetaZeroLine->SetLineColor(kBlue);
    thetaZeroLine->SetLineWidth(2);
    thetaZeroLine->SetLineStyle(2);
    
    phiZeroLine = new TPolyLine3D(2);
    phiZeroLine->SetPoint(0, 48, 4, 12);
    phiZeroLine->SetPoint(1, 0, 4, 12);
    phiZeroLine->SetLineColor(kGreen);
    phiZeroLine->SetLineWidth(2);
    phiZeroLine->SetLineStyle(2);
    
    if (trajectory.valid) {
        double zMin = 0;
        double zMax = 48;
        
        for (const auto& hit : hits3D) {
            if (hit.confidence >= 0.5) {
                zMin = std::min(zMin, hit.z - 2.0);
                zMax = std::max(zMax, hit.z + 2.0);
            }
        }
        
        if (std::abs(trajectory.vz) > 0.01) {
            trajX_vs_Z = new TF1("trajX_vs_Z", 
                                Form("%f + (x - %f) * %f / %f", 
                                     trajectory.x0, trajectory.z0, trajectory.vx, trajectory.vz),
                                zMin, zMax);
            trajX_vs_Z->SetLineColor(kRed);
            trajX_vs_Z->SetLineWidth(2);
            
            trajY_vs_Z = new TF1("trajY_vs_Z", 
                                Form("%f + (x - %f) * %f / %f", 
                                     trajectory.y0, trajectory.z0, trajectory.vy, trajectory.vz),
                                zMin, zMax);
            trajY_vs_Z->SetLineColor(kRed);
            trajY_vs_Z->SetLineWidth(2);
        }
        
        // Create polyline for 3D visualization
        int nPoints = 100;
        trajectoryLine = new TPolyLine3D(nPoints);
        
        double tMin = 1e9, tMax = -1e9;
        for (const auto& hit : hits3D) {
            if (hit.confidence >= 0.8) {
                double t = (hit.x - trajectory.x0) * trajectory.vx + 
                          (hit.y - trajectory.y0) * trajectory.vy + 
                          (hit.z - trajectory.z0) * trajectory.vz;
                tMin = std::min(tMin, t);
                tMax = std::max(tMax, t);
            }
        }
        
        double extension = 3.0;
        tMin -= extension;
        tMax += extension;
        
        for (int i = 0; i < nPoints; i++) {
            double t = tMin + (tMax - tMin) * i / (nPoints - 1);
            double x = trajectory.x0 + t * trajectory.vx;
            double y = trajectory.y0 + t * trajectory.vy;
            double z = trajectory.z0 + t * trajectory.vz;
            
            trajectoryLine->SetPoint(i, z, y, x);
        }
        
        trajectoryLine->SetLineColor(kRed);
        trajectoryLine->SetLineWidth(3);
    }

    // Write histograms
    h3D_Full->Write();
    h3D_GroupResolved->Write();
    h3D_Partial->Write();
    h3D_All->Write();
    h3D_Alpha->Write();
    
    if (trajectoryLine) {
        trajectoryLine->Write("trajectory");
    }
    if (trajX_vs_Z) {
        trajX_vs_Z->Write();
    }
    if (trajY_vs_Z) {
        trajY_vs_Z->Write();
    }
    
    thetaZeroLine->Write("thetaZeroLine");
    phiZeroLine->Write("phiZeroLine");
    
    // Save tree with individual view charges
    TTree* tree3D = new TTree("Hits3D", "Reconstructed 3D Hits with Individual View Charges");
    
    Double_t t_x, t_y, t_z, t_charge, t_confidence;
    Double_t t_timeXY, t_timeXZ, t_timeZY;
    Double_t t_chargeXY, t_chargeXZ, t_chargeZY;
    Double_t t_originalChargeXY, t_originalChargeXZ, t_originalChargeZY;
    Int_t t_nViews, t_redistributed;
    Char_t t_matchType[100];
    
    tree3D->Branch("x", &t_x, "x/D");
    tree3D->Branch("y", &t_y, "y/D");
    tree3D->Branch("z", &t_z, "z/D");
    tree3D->Branch("charge", &t_charge, "charge/D");
    tree3D->Branch("confidence", &t_confidence, "confidence/D");
    tree3D->Branch("timeXY", &t_timeXY, "timeXY/D");
    tree3D->Branch("timeXZ", &t_timeXZ, "timeXZ/D");
    tree3D->Branch("timeZY", &t_timeZY, "timeZY/D");
    tree3D->Branch("chargeXY", &t_chargeXY, "chargeXY/D");
    tree3D->Branch("chargeXZ", &t_chargeXZ, "chargeXZ/D");
    tree3D->Branch("chargeZY", &t_chargeZY, "chargeZY/D");
    tree3D->Branch("originalChargeXY", &t_originalChargeXY, "originalChargeXY/D");
    tree3D->Branch("originalChargeXZ", &t_originalChargeXZ, "originalChargeXZ/D");
    tree3D->Branch("originalChargeZY", &t_originalChargeZY, "originalChargeZY/D");
    tree3D->Branch("nViews", &t_nViews, "nViews/I");
    tree3D->Branch("redistributed", &t_redistributed, "redistributed/I");
    tree3D->Branch("matchType", t_matchType, "matchType/C");
    
    for (const auto& hit : hits3D) {
        t_x = hit.x;
        t_y = hit.y;
        t_z = hit.z;
        t_charge = hit.charge;
        t_confidence = hit.confidence;
        t_timeXY = hit.timeXY;
        t_timeXZ = hit.timeXZ;
        t_timeZY = hit.timeZY;
        t_chargeXY = hit.chargeXY;
        t_chargeXZ = hit.chargeXZ;
        t_chargeZY = hit.chargeZY;
        t_originalChargeXY = hit.originalChargeXY;
        t_originalChargeXZ = hit.originalChargeXZ;
        t_originalChargeZY = hit.originalChargeZY;
        t_nViews = hit.nViews;
        t_redistributed = hit.redistributed ? 1 : 0;
        strcpy(t_matchType, hit.matchType.c_str());
        
        tree3D->Fill();
    }
    
    tree3D->Write();
    
    // Save trajectory parameters with full PCA analysis
    TTree* trajTree = new TTree("Trajectory", "Fitted Trajectory Parameters with PCA Analysis");
    
    // Original trajectory parameters
    Double_t tr_x0, tr_y0, tr_z0, tr_vx, tr_vy, tr_vz, tr_theta, tr_phi;
    Int_t tr_valid, tr_nHits, tr_isSingleTrack;
    
    // PCA parameters 
    Double_t pca_centroid_x, pca_centroid_y, pca_centroid_z;
    Double_t pca_pc1_x, pca_pc1_y, pca_pc1_z;
    Double_t pca_pc2_x, pca_pc2_y, pca_pc2_z;
    Double_t pca_pc3_x, pca_pc3_y, pca_pc3_z;
    Double_t pca_eigenvalues[3], pca_explainedVariance[3], pca_cumulativeVariance[3];
    Double_t pca_linearity, pca_planarity, pca_trackWidth, pca_trackThickness, pca_totalVariance;
    Double_t pca_rmsResidual, pca_maxResidual, pca_confidenceScore;
    Int_t pca_nHighConfidenceHits, pca_nOutliers;
    Double_t pca_weightedCentroidDistance, pca_outlierThreshold;
    Double_t pca_meanLongitudinalDistance, pca_meanTransverseDistance, pca_meanNormalDistance;
    Int_t pca_valid;
    
    // Branch the original parameters
    trajTree->Branch("x0", &tr_x0, "x0/D");
    trajTree->Branch("y0", &tr_y0, "y0/D");
    trajTree->Branch("z0", &tr_z0, "z0/D");
    trajTree->Branch("vx", &tr_vx, "vx/D");
    trajTree->Branch("vy", &tr_vy, "vy/D");
    trajTree->Branch("vz", &tr_vz, "vz/D");
    trajTree->Branch("theta", &tr_theta, "theta/D");
    trajTree->Branch("phi", &tr_phi, "phi/D");
    trajTree->Branch("valid", &tr_valid, "valid/I");
    trajTree->Branch("nHits", &tr_nHits, "nHits/I");
    trajTree->Branch("isSingleTrack", &tr_isSingleTrack, "isSingleTrack/I");
    
    // Branch the PCA parameters
    trajTree->Branch("pca_valid", &pca_valid, "pca_valid/I");
    trajTree->Branch("pca_centroid_x", &pca_centroid_x, "pca_centroid_x/D");
    trajTree->Branch("pca_centroid_y", &pca_centroid_y, "pca_centroid_y/D");
    trajTree->Branch("pca_centroid_z", &pca_centroid_z, "pca_centroid_z/D");
    trajTree->Branch("pca_pc1_x", &pca_pc1_x, "pca_pc1_x/D");
    trajTree->Branch("pca_pc1_y", &pca_pc1_y, "pca_pc1_y/D");
    trajTree->Branch("pca_pc1_z", &pca_pc1_z, "pca_pc1_z/D");
    trajTree->Branch("pca_pc2_x", &pca_pc2_x, "pca_pc2_x/D");
    trajTree->Branch("pca_pc2_y", &pca_pc2_y, "pca_pc2_y/D");
    trajTree->Branch("pca_pc2_z", &pca_pc2_z, "pca_pc2_z/D");
    trajTree->Branch("pca_pc3_x", &pca_pc3_x, "pca_pc3_x/D");
    trajTree->Branch("pca_pc3_y", &pca_pc3_y, "pca_pc3_y/D");
    trajTree->Branch("pca_pc3_z", &pca_pc3_z, "pca_pc3_z/D");
    trajTree->Branch("pca_eigenvalues", pca_eigenvalues, "pca_eigenvalues[3]/D");
    trajTree->Branch("pca_explainedVariance", pca_explainedVariance, "pca_explainedVariance[3]/D");
    trajTree->Branch("pca_cumulativeVariance", pca_cumulativeVariance, "pca_cumulativeVariance[3]/D");
    trajTree->Branch("pca_linearity", &pca_linearity, "pca_linearity/D");
    trajTree->Branch("pca_planarity", &pca_planarity, "pca_planarity/D");
    trajTree->Branch("pca_trackWidth", &pca_trackWidth, "pca_trackWidth/D");
    trajTree->Branch("pca_trackThickness", &pca_trackThickness, "pca_trackThickness/D");
    trajTree->Branch("pca_totalVariance", &pca_totalVariance, "pca_totalVariance/D");
    trajTree->Branch("pca_rmsResidual", &pca_rmsResidual, "pca_rmsResidual/D");
    trajTree->Branch("pca_maxResidual", &pca_maxResidual, "pca_maxResidual/D");
    trajTree->Branch("pca_confidenceScore", &pca_confidenceScore, "pca_confidenceScore/D");
    trajTree->Branch("pca_nHighConfidenceHits", &pca_nHighConfidenceHits, "pca_nHighConfidenceHits/I");
    trajTree->Branch("pca_nOutliers", &pca_nOutliers, "pca_nOutliers/I");
    trajTree->Branch("pca_weightedCentroidDistance", &pca_weightedCentroidDistance, "pca_weightedCentroidDistance/D");
    trajTree->Branch("pca_outlierThreshold", &pca_outlierThreshold, "pca_outlierThreshold/D");
    trajTree->Branch("pca_meanLongitudinalDistance", &pca_meanLongitudinalDistance, "pca_meanLongitudinalDistance/D");
    trajTree->Branch("pca_meanTransverseDistance", &pca_meanTransverseDistance, "pca_meanTransverseDistance/D");
    trajTree->Branch("pca_meanNormalDistance", &pca_meanNormalDistance, "pca_meanNormalDistance/D");
    
    // Fill the tree with both trajectory and PCA data
    if (trajectory.valid) {
        // Fill trajectory parameters
        tr_x0 = trajectory.x0;
        tr_y0 = trajectory.y0;
        tr_z0 = trajectory.z0;
        tr_vx = trajectory.vx;
        tr_vy = trajectory.vy;
        tr_vz = trajectory.vz;
        tr_theta = trajectory.theta;
        tr_phi = trajectory.phi;
        tr_valid = 1;
        tr_isSingleTrack = isSingleTrackEvent(pcaResult) ? 1 : 0;
        
        int highConfHits = 0;
        for (const auto& hit : hits3D) {
            if (hit.confidence >= 0.5) highConfHits++;
        }
        tr_nHits = highConfHits;
        
        // Fill PCA parameters
        if (pcaResult.valid) {
            pca_valid = 1;
            pca_centroid_x = pcaResult.centroid_x;
            pca_centroid_y = pcaResult.centroid_y;
            pca_centroid_z = pcaResult.centroid_z;
            pca_pc1_x = pcaResult.pc1_x;
            pca_pc1_y = pcaResult.pc1_y;
            pca_pc1_z = pcaResult.pc1_z;
            pca_pc2_x = pcaResult.pc2_x;
            pca_pc2_y = pcaResult.pc2_y;
            pca_pc2_z = pcaResult.pc2_z;
            pca_pc3_x = pcaResult.pc3_x;
            pca_pc3_y = pcaResult.pc3_y;
            pca_pc3_z = pcaResult.pc3_z;
            
            for (int i = 0; i < 3; i++) {
                pca_eigenvalues[i] = pcaResult.eigenvalues[i];
                pca_explainedVariance[i] = pcaResult.explainedVariance[i];
                pca_cumulativeVariance[i] = pcaResult.cumulativeVariance[i];
            }
            
            pca_linearity = pcaResult.linearity;
            pca_planarity = pcaResult.planarity;
            pca_trackWidth = pcaResult.trackWidth;
            pca_trackThickness = pcaResult.trackThickness;
            pca_totalVariance = pcaResult.totalVariance;
            pca_rmsResidual = pcaResult.rmsResidual;
            pca_maxResidual = pcaResult.maxResidual;
            pca_confidenceScore = pcaResult.confidenceScore;
            pca_nHighConfidenceHits = pcaResult.nHighConfidenceHits;
            pca_nOutliers = pcaResult.nOutliers;
            pca_weightedCentroidDistance = pcaResult.weightedCentroidDistance;
            pca_outlierThreshold = pcaResult.outlierThreshold;
            pca_meanLongitudinalDistance = pcaResult.meanLongitudinalDistance;
            pca_meanTransverseDistance = pcaResult.meanTransverseDistance;
            pca_meanNormalDistance = pcaResult.meanNormalDistance;
        } else {
            // Set PCA parameters to default values if PCA failed
            pca_valid = 0;
            pca_centroid_x = pca_centroid_y = pca_centroid_z = 0;
            pca_pc1_x = pca_pc1_y = pca_pc1_z = 0;
            pca_pc2_x = pca_pc2_y = pca_pc2_z = 0;
            pca_pc3_x = pca_pc3_y = pca_pc3_z = 0;
            for (int i = 0; i < 3; i++) {
                pca_eigenvalues[i] = 0;
                pca_explainedVariance[i] = 0;
                pca_cumulativeVariance[i] = 0;
            }
            pca_linearity = pca_planarity = pca_trackWidth = pca_trackThickness = 0;
            pca_totalVariance = pca_rmsResidual = pca_maxResidual = pca_confidenceScore = 0;
            pca_nHighConfidenceHits = pca_nOutliers = 0;
            pca_weightedCentroidDistance = pca_outlierThreshold = 0;
            pca_meanLongitudinalDistance = pca_meanTransverseDistance = pca_meanNormalDistance = 0;
        }
        
        trajTree->Fill();
    }
    
    trajTree->Write();
    
    // Save summary information
    TNamed eventInfo("EventInfo", Form("Event %d: %zu 3D hits reconstructed (%d high-confidence hits used for trajectory)", 
                                      eventNumber, hits3D.size(), 
                                      trajectory.valid ? tr_nHits : 0));
    eventInfo.Write();
    
    if (trajectory.valid) {
        TNamed trajInfo("TrajectoryInfo", Form("Theta=%.1f deg, Phi=%.1f deg, Entry=(%.1f,%.1f,%.1f)", 
                                              trajectory.theta, trajectory.phi, 
                                              trajectory.x0, trajectory.y0, trajectory.z0));
        trajInfo.Write();
    }
    
    fOutput->Close();
    delete fOutput;
}

// Main entry point function for batch processing with gamma/beta filtering
void EventDisplays_3D_Side_Batch_Filtered() {
    gROOT->SetBatch(kTRUE);
    
    // Get command line arguments
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    calibManager.LoadCalibration("../../src/tools/calibration_files/gain_calibration_strict_68813_time3.txt");

    // Check arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_events.root>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /path/to/events.root" << std::endl;
        std::cerr << "This will process all events and save only filtered gamma/beta events" << std::endl;
        return;
    }
    
    TString inputFile = argv[1];
    
    std::cout << "=== Batch Processing Mode with Gamma/Beta Filtering ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Filters:" << std::endl;
    std::cout << "  - Gamma events: Only X=23 perimeter hits AND X=23 hits have Z > 23" << std::endl;
    std::cout << "  - Beta events: Only X=23 perimeter hits AND X=23 hits have Z ≤ 23" << std::endl;
    std::cout << "  - Minimum 2D hits per event: 5" << std::endl;
    std::cout << "  - Valid trajectory fitting required" << std::endl;
    
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
    int skippedTooFewHits = 0;
    int skippedNoTrajectory = 0;
    int skippedPerimeterFilter = 0;
    int gammaEvents = 0;
    int betaEvents = 0;
    
    auto startTime = std::chrono::steady_clock::now();
    
    std::cout << "\n=== Processing Events ===" << std::endl;
    
    // Process all events
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
        tree->GetEntry(iEvent);
        
        // Progress reporting
        if (iEvent % 100 == 0 || iEvent == nEvents - 1) {
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double progress = 100.0 * iEvent / nEvents;
            double rate = iEvent > 0 ? (double)iEvent / elapsed : 0;
            
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% "
                      << "(" << iEvent << "/" << nEvents << ") "
                      << "Rate: " << std::setprecision(1) << rate << " events/sec "
                      << "Accepted: " << processedEvents << " (γ:" << gammaEvents 
                      << ", β:" << betaEvents << ")";
            std::cout.flush();
        }
        
        // Filter 1: Skip events with < 5 hits
        if (event->GetNHits() < 5) {
            skippedTooFewHits++;
            continue;
        }
        
        // Reconstruct 3D hits
        std::vector<Hit3D> hits3D = reconstructHits(event);
        
        // Filter 2: Check perimeter criteria
        PerimeterResult perimeterResult = analyzePerimeterVoxels(hits3D);
        if (!perimeterResult.passesFilter) {
            skippedPerimeterFilter++;
            continue;
        }
        
        // Fit trajectory using high-confidence hits
        TrajectoryParams trajectory = fitTrajectory(hits3D);
        PCATrajectoryResult pcaResult = performPCA(hits3D, 0.5);
        
        // Filter 3: Valid trajectory required
        if (!trajectory.valid) {
            skippedNoTrajectory++;
            continue;
        }
        
        // Event passed all filters - save it
        SaveFilteredEvent(event, iEvent, hits3D, trajectory, pcaResult, perimeterResult, inputFile);
        
        processedEvents++;
        if (perimeterResult.isGammaEvent) {
            gammaEvents++;
        } else if (perimeterResult.isBetaEvent) {
            betaEvents++;
        }
        
        // Print summary for accepted events
        printSummary(hits3D, event->GetEventID());
    }
    
    auto endTime = std::chrono::steady_clock::now();
    auto totalTime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    
    std::cout << "\n\n=== Processing Complete ===" << std::endl;
    std::cout << "Total processing time: " << totalTime << " seconds" << std::endl;
    std::cout << "Total events: " << nEvents << std::endl;
    std::cout << "Accepted events: " << processedEvents << std::endl;
    std::cout << "  - Gamma events: " << gammaEvents << std::endl;
    std::cout << "  - Beta events: " << betaEvents << std::endl;
    std::cout << "Rejected events: " << (nEvents - processedEvents) << std::endl;
    std::cout << "  - Too few hits (< 5): " << skippedTooFewHits << std::endl;
    std::cout << "  - Failed perimeter filter: " << skippedPerimeterFilter << std::endl;
    std::cout << "  - No valid trajectory: " << skippedNoTrajectory << std::endl;
    std::cout << "Processing efficiency: " << std::fixed << std::setprecision(1) 
              << (100.0 * processedEvents / nEvents) << "%" << std::endl;
    
    // Create summary file
    TString summaryFile = inputFile;
    summaryFile.ReplaceAll(".root", "_FilteredEvents_summary.txt");
    std::ofstream summary(summaryFile.Data());
    summary << "Filtered Events Processing Summary\n";
    summary << "=================================\n";
    summary << "Input file: " << inputFile << "\n";
    summary << "Total events: " << nEvents << "\n";
    summary << "Accepted events: " << processedEvents << "\n";
    summary << "  - Gamma events: " << gammaEvents << "\n";
    summary << "  - Beta events: " << betaEvents << "\n";
    summary << "Rejected events: " << (nEvents - processedEvents) << "\n";
    summary << "  - Too few hits (< 5): " << skippedTooFewHits << "\n";
    summary << "  - Failed perimeter filter: " << skippedPerimeterFilter << "\n";
    summary << "  - No valid trajectory: " << skippedNoTrajectory << "\n";
    summary << "Processing efficiency: " << std::fixed << std::setprecision(1) << (100.0 * processedEvents / nEvents) << "%\n";
    summary << "Total processing time: " << totalTime << " seconds\n";
    summary << "\nSelection criteria:\n";
    summary << "- Minimum 2D hits: 5\n";
    summary << "- Valid trajectory fitting required\n";
    summary << "- Gamma events: Only X=23 perimeter hits AND all hits Z > 23\n";
    summary << "- Beta events: Only X=23 perimeter hits AND all hits Z ≤ 23\n";
    summary << "\nOutput directories:\n";
    summary << "- Gamma events: <input_dir>/GammaEvents_OnlyX23_ZRight/\n";
    summary << "- Beta events: <input_dir>/BetaEvents_OnlyX23_ZLeft/\n";
    summary << "\nFile naming convention:\n";
    summary << "- <basename>_3D_event<N>_<Gamma|Beta>_<SingleTrack|MultiTrack>.root\n";
    
    summary.close();
    
    std::cout << "Summary written to: " << summaryFile << std::endl;
    std::cout << "\nOutput directories created:" << std::endl;
    std::cout << "- Gamma events: " << gSystem->DirName(inputFile) << "/GammaEvents_OnlyX23_ZRight/" << std::endl;
    std::cout << "- Beta events: " << gSystem->DirName(inputFile) << "/BetaEvents_OnlyX23_ZLeft/" << std::endl;
    std::cout << "\nTo view the 3D histograms in ROOT:" << std::endl;
    std::cout << "root <output_file>.root" << std::endl;
    std::cout << "h3D_All->Draw(\"BOX2Z\")  // All hits" << std::endl;
    std::cout << "h3D_Full->Draw(\"BOX2Z\") // Perfect 3-view hits only" << std::endl;
    std::cout << "h3D_Alpha->Draw(\"BOX2Z\") // Hit confidence levels" << std::endl;
    std::cout << "\nTo view with trajectory overlay:" << std::endl;
    std::cout << "h3D_All->Draw(\"BOX2Z\")" << std::endl;
    std::cout << "TPolyLine3D *traj = (TPolyLine3D*)gDirectory->Get(\"trajectory\")" << std::endl;
    std::cout << "if(traj) { traj->SetLineColor(2); traj->SetLineWidth(3); traj->Draw(); }" << std::endl;
    std::cout << "TPolyLine3D *theta0 = (TPolyLine3D*)gDirectory->Get(\"thetaZeroLine\")" << std::endl;
    std::cout << "if(theta0) theta0->Draw()" << std::endl;
    std::cout << "TPolyLine3D *phi0 = (TPolyLine3D*)gDirectory->Get(\"phiZeroLine\")" << std::endl;
    std::cout << "if(phi0) phi0->Draw()" << std::endl;
    
    // Clean up
    file->Close();
    delete event;
    exit(0);
}