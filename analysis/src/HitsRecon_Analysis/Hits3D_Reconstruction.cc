#ifndef HIT3D_RECONSTRUCTION_CC
#define HIT3D_RECONSTRUCTION_CC

#include "Hit3DStructures.hh"
#include "../tools/global_header.hh"
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>

// CalibrationData implementation
struct CalibrationData {
    std::map<std::string, std::map<std::pair<int,int>, double>> gainMap;
    
    bool LoadCalibration(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "Error: Cannot open calibration file " << filename << std::endl;
            return false;
        }
        
        std::string view;
        int coord1, coord2;
        double gain;
        
        while (file >> view >> coord1 >> coord2 >> gain) {
            if (gain > 0) {
                gainMap[view][std::make_pair(coord1, coord2)] = gain;
            }
        }
        file.close();
        
        std::cout << "Loaded calibration for " << gainMap.size() << " views" << std::endl;
        return true;
    }
    
    double GetGain(const std::string& view, int coord1, int coord2) {
        auto viewIt = gainMap.find(view);
        if (viewIt != gainMap.end()) {
            auto coordIt = viewIt->second.find(std::make_pair(coord1, coord2));
            if (coordIt != viewIt->second.end()) {
                return coordIt->second;
            }
        }
        return -1.0;
    }
};

// Helper function to get calibrated charge
double getCalibratedCharge(Hit* hit, CalibrationData& calibData) {
    std::string viewName;
    int coord1, coord2;

    if (hit->GetView() == 0) {
        viewName = "XY";
        coord1 = hit->GetX();
        coord2 = hit->GetY();
    }
    else if (hit->GetView() == 1) {
        viewName = "XZ";
        coord1 = hit->GetX();
        coord2 = hit->GetZ();
    }
    else if (hit->GetView() == 2) {
        viewName = "ZY";
        coord1 = hit->GetZ();
        coord2 = hit->GetY();
    }
    else {
        return hit->GetPE();
    }

    double gain = calibData.GetGain(viewName, coord1, coord2);
    if (gain > 0) {
        double originalCharge = hit->GetPE();
        double calibratedCharge = originalCharge / (3*gain);
        return calibratedCharge;
    } else {
        return hit->GetPE();
    }
}

// Utility functions
double getViewCharge(const std::string& view, double x, double y, double z, 
                    const ViewHits& viewXY, const ViewHits& viewXZ, const ViewHits& viewZY) {
    if (view == "XY") {
        int targetX = (int)x;
        int targetY = (int)y;
        for (size_t k = 0; k < viewXY.positions.size(); k++) {
            if (viewXY.positions[k].first == targetX && viewXY.positions[k].second == targetY) {
                return viewXY.charges[k];
            }
        }
    } else if (view == "XZ") {
        int targetX = (int)x;
        int targetZ = (int)z;
        for (size_t k = 0; k < viewXZ.positions.size(); k++) {
            if (viewXZ.positions[k].first == targetX && viewXZ.positions[k].second == targetZ) {
                return viewXZ.charges[k];
            }
        }
    } else if (view == "ZY") {
        int targetZ = (int)z;
        int targetY = (int)y;
        for (size_t k = 0; k < viewZY.positions.size(); k++) {
            if (viewZY.positions[k].first == targetZ && viewZY.positions[k].second == targetY) {
                return viewZY.charges[k];
            }
        }
    }
    return 0.0;
}

double getViewTiming(const std::string& view, double x, double y, double z, 
                    const ViewHits& viewXY, const ViewHits& viewXZ, const ViewHits& viewZY) {
    if (view == "XY") {
        int targetX = (int)x;
        int targetY = (int)y;
        for (size_t k = 0; k < viewXY.positions.size(); k++) {
            if (viewXY.positions[k].first == targetX && viewXY.positions[k].second == targetY) {
                return viewXY.times[k];
            }
        }
    } else if (view == "XZ") {
        int targetX = (int)x;
        int targetZ = (int)z;
        for (size_t k = 0; k < viewXZ.positions.size(); k++) {
            if (viewXZ.positions[k].first == targetX && viewXZ.positions[k].second == targetZ) {
                return viewXZ.times[k];
            }
        }
    } else if (view == "ZY") {
        int targetZ = (int)z;
        int targetY = (int)y;
        for (size_t k = 0; k < viewZY.positions.size(); k++) {
            if (viewZY.positions[k].first == targetZ && viewZY.positions[k].second == targetY) {
                return viewZY.times[k];
            }
        }
    }
    return -1.0;
}

// Forward declaration
std::vector<Hit3D> resolveAdjacentHits(const std::vector<Hit3D>& initialHits, 
                                      const ViewHits& viewXY, 
                                      const ViewHits& viewXZ, 
                                      const ViewHits& viewZY);

// Main reconstruction function
std::vector<Hit3D> reconstructHits(Event* event, CalibrationData& calibData) {
    std::vector<Hit3D> hits3D;
    
    ViewHits viewXY, viewXZ, viewZY;
    
    TClonesArray* hits = event->GetHits();
    if (!hits) return hits3D;
    
    // Sort hits into views
    for (int i = 0; i < hits->GetEntries(); i++) {
        Hit* hit = (Hit*)hits->At(i);
        if (!hit) continue;

        double calibratedCharge = getCalibratedCharge(hit, calibData);

        if (hit->GetView() == 0) {
            viewXY.indices.push_back(i);
            viewXY.positions.push_back({hit->GetX(), hit->GetY()});
            viewXY.charges.push_back(calibratedCharge);
            viewXY.times.push_back(hit->GetTfromSpill());
        }
        else if (hit->GetView() == 1) {
            viewXZ.indices.push_back(i);
            viewXZ.positions.push_back({hit->GetX(), hit->GetZ()});
            viewXZ.charges.push_back(calibratedCharge);
            viewXZ.times.push_back(hit->GetTfromSpill());
        }
        else if (hit->GetView() == 2) {
            viewZY.indices.push_back(i);
            viewZY.positions.push_back({hit->GetZ(), hit->GetY()});
            viewZY.charges.push_back(calibratedCharge);
            viewZY.times.push_back(hit->GetTfromSpill());
        }
    }
    
    std::vector<bool> usedXY(viewXY.indices.size(), false);
    std::vector<bool> usedXZ(viewXZ.indices.size(), false);
    std::vector<bool> usedZY(viewZY.indices.size(), false);
    
    // First pass: Find 3-view matches
    for (size_t ixy = 0; ixy < viewXY.positions.size(); ixy++) {
        if (usedXY[ixy]) continue;
        
        int x_xy = viewXY.positions[ixy].first;
        int y_xy = viewXY.positions[ixy].second;
        bool foundMatch = false;
        
        for (size_t ixz = 0; ixz < viewXZ.positions.size(); ixz++) {
            if (usedXZ[ixz]) continue;
            
            int x_xz = viewXZ.positions[ixz].first;
            int z_xz = viewXZ.positions[ixz].second;
            
            if (x_xy != x_xz) continue;
            
            for (size_t izy = 0; izy < viewZY.positions.size(); izy++) {
                if (usedZY[izy]) continue;
                
                int z_zy = viewZY.positions[izy].first;
                int y_zy = viewZY.positions[izy].second;
                
                if (y_xy == y_zy && z_xz == z_zy) {
                    Hit3D hit3d;
                    hit3d.x = x_xy + 0.5;
                    hit3d.y = y_xy + 0.5;
                    hit3d.z = z_xz + 0.5;
                    hit3d.charge = (viewXY.charges[ixy] + viewXZ.charges[ixz] + viewZY.charges[izy]);
                    hit3d.timeXY = viewXY.times[ixy];
                    hit3d.timeXZ = viewXZ.times[ixz];
                    hit3d.timeZY = viewZY.times[izy];
                    hit3d.chargeXY = viewXY.charges[ixy];
                    hit3d.chargeXZ = viewXZ.charges[ixz];
                    hit3d.chargeZY = viewZY.charges[izy];
                    hit3d.originalChargeXY = viewXY.charges[ixy];
                    hit3d.originalChargeXZ = viewXZ.charges[ixz];
                    hit3d.originalChargeZY = viewZY.charges[izy];
                    hit3d.nViews = 3;
                    hit3d.matchType = "Full";
                    hit3d.confidence = 1.0;
                    hit3d.redistributed = false;
                    
                    hits3D.push_back(hit3d);
                    usedXY[ixy] = true;
                    usedXZ[ixz] = true;
                    usedZY[izy] = true;
                    foundMatch = true;
                    break;
                }
            }
            if (foundMatch) break;
        }
    }
    
    // Second pass: Find 2-view matches (similar pattern for all combinations)
    // ... [Include all the 2-view matching code here]
    for (size_t ixy = 0; ixy < viewXY.positions.size(); ixy++) {
        if (usedXY[ixy]) continue;
        
        int x_xy = viewXY.positions[ixy].first;
        int y_xy = viewXY.positions[ixy].second;
        
        for (size_t ixz = 0; ixz < viewXZ.positions.size(); ixz++) {
            if (usedXZ[ixz]) continue;
            
            int x_xz = viewXZ.positions[ixz].first;
            int z_xz = viewXZ.positions[ixz].second;
            
            if (x_xy == x_xz) {
                Hit3D hit3d;
                hit3d.x = x_xy + 0.5;
                hit3d.y = y_xy + 0.5;
                hit3d.z = z_xz + 0.5;
                hit3d.charge = (viewXY.charges[ixy] + viewXZ.charges[ixz]);
                hit3d.timeXY = viewXY.times[ixy];
                hit3d.timeXZ = viewXZ.times[ixz];
                hit3d.timeZY = -1;
                // Initialize individual view charges
                hit3d.chargeXY = viewXY.charges[ixy];
                hit3d.chargeXZ = viewXZ.charges[ixz];
                hit3d.chargeZY = -1;  // Missing view
                // Initialize original charges
                hit3d.originalChargeXY = viewXY.charges[ixy];
                hit3d.originalChargeXZ = viewXZ.charges[ixz];
                hit3d.originalChargeZY = -1;  // Never had this view
                hit3d.nViews = 2;
                hit3d.matchType = "Partial-XY-XZ";
                hit3d.confidence = 0.1;
                hit3d.redistributed = false;
                
                hits3D.push_back(hit3d);
                usedXY[ixy] = true;
                usedXZ[ixz] = true;
                break;
            }
        }
    }
    
    // XY-ZY matches
    for (size_t ixy = 0; ixy < viewXY.positions.size(); ixy++) {
        if (usedXY[ixy]) continue;
        
        int x_xy = viewXY.positions[ixy].first;
        int y_xy = viewXY.positions[ixy].second;
        
        for (size_t izy = 0; izy < viewZY.positions.size(); izy++) {
            if (usedZY[izy]) continue;
            
            int z_zy = viewZY.positions[izy].first;
            int y_zy = viewZY.positions[izy].second;
            
            if (y_xy == y_zy) {
                Hit3D hit3d;
                hit3d.x = x_xy + 0.5;
                hit3d.y = y_xy + 0.5;
                hit3d.z = z_zy + 0.5;
                hit3d.charge = (viewXY.charges[ixy] + viewZY.charges[izy]);
                hit3d.timeXY = viewXY.times[ixy];
                hit3d.timeXZ = -1;
                hit3d.timeZY = viewZY.times[izy];
                // Initialize individual view charges
                hit3d.chargeXY = viewXY.charges[ixy];
                hit3d.chargeXZ = -1;  // Missing view
                hit3d.chargeZY = viewZY.charges[izy];
                // Initialize original charges
                hit3d.originalChargeXY = viewXY.charges[ixy];
                hit3d.originalChargeXZ = -1;  // Never had this view
                hit3d.originalChargeZY = viewZY.charges[izy];
                hit3d.nViews = 2;
                hit3d.matchType = "Partial-XY-ZY";
                hit3d.confidence = 0.1;
                hit3d.redistributed = false;
                
                hits3D.push_back(hit3d);
                usedXY[ixy] = true;
                usedZY[izy] = true;
                break;
            }
        }
    }
    
    // XZ-ZY matches
    for (size_t ixz = 0; ixz < viewXZ.positions.size(); ixz++) {
        if (usedXZ[ixz]) continue;
        
        int x_xz = viewXZ.positions[ixz].first;
        int z_xz = viewXZ.positions[ixz].second;

        for (size_t izy = 0; izy < viewZY.positions.size(); izy++) {
            if (usedZY[izy]) continue;
            
            int z_zy = viewZY.positions[izy].first;
            int y_zy = viewZY.positions[izy].second;
            
            if (z_xz == z_zy) {
                Hit3D hit3d;
                hit3d.x = x_xz + 0.5;
                hit3d.y = y_zy + 0.5;
                hit3d.z = z_xz + 0.5;
                hit3d.charge = (viewXZ.charges[ixz] + viewZY.charges[izy]);
                hit3d.timeXY = -1;
                hit3d.timeXZ = viewXZ.times[ixz];
                hit3d.timeZY = viewZY.times[izy];
                // Initialize individual view charges
                hit3d.chargeXY = -1;  // Missing view
                hit3d.chargeXZ = viewXZ.charges[ixz];
                hit3d.chargeZY = viewZY.charges[izy];
                // Initialize original charges
                hit3d.originalChargeXY = -1;  // Never had this view
                hit3d.originalChargeXZ = viewXZ.charges[ixz];
                hit3d.originalChargeZY = viewZY.charges[izy];
                hit3d.nViews = 2;
                hit3d.matchType = "Partial-XZ-ZY";
                hit3d.confidence = 0.1;
                hit3d.redistributed = false;
                
                hits3D.push_back(hit3d);
                usedXZ[ixz] = true;
                usedZY[izy] = true;
                break;
            }
        }
    }
    // Apply resolution
    hits3D = resolveAdjacentHits(hits3D, viewXY, viewXZ, viewZY);
    
    return hits3D;
}

// Resolution function
std::vector<Hit3D> resolveAdjacentHits(const std::vector<Hit3D>& initialHits, 
                                      const ViewHits& viewXY, 
                                      const ViewHits& viewXZ, 
                                      const ViewHits& viewZY) {
    // ... [Include the full resolution algorithm here]
    std::vector<Hit3D> resolvedHits = initialHits;
    // Step 1: Find all overlapping pairs and form groups
    struct OverlapGroup {
        std::vector<int> hitIndices;
        std::string sharedView;  // The view where all hits in this group overlap
        int x, y, z;  // The shared coordinates
    };
    
    std::vector<OverlapGroup> finalGroups;
    
    // For each unique shared position, collect all hits
    std::map<std::string, OverlapGroup> positionGroups;
    
    for (size_t i = 0; i < initialHits.size(); i++) {
        const Hit3D& hit = initialHits[i];
        
        // Create keys for each view position
        std::string keyXY = "XY_" + std::to_string((int)hit.x) + "_" + std::to_string((int)hit.y);
        std::string keyXZ = "XZ_" + std::to_string((int)hit.x) + "_" + std::to_string((int)hit.z);
        std::string keyZY = "ZY_" + std::to_string((int)hit.z) + "_" + std::to_string((int)hit.y);
        
        // Add to XY group
        if (positionGroups.find(keyXY) == positionGroups.end()) {
            positionGroups[keyXY].sharedView = "XY";
            positionGroups[keyXY].x = (int)hit.x;
            positionGroups[keyXY].y = (int)hit.y;
            positionGroups[keyXY].z = -1;
        }
        positionGroups[keyXY].hitIndices.push_back(i);
        
        // Add to XZ group
        if (positionGroups.find(keyXZ) == positionGroups.end()) {
            positionGroups[keyXZ].sharedView = "XZ";
            positionGroups[keyXZ].x = (int)hit.x;
            positionGroups[keyXZ].y = -1;
            positionGroups[keyXZ].z = (int)hit.z;
        }
        positionGroups[keyXZ].hitIndices.push_back(i);
        
        // Add to ZY group
        if (positionGroups.find(keyZY) == positionGroups.end()) {
            positionGroups[keyZY].sharedView = "ZY";
            positionGroups[keyZY].x = -1;
            positionGroups[keyZY].y = (int)hit.y;
            positionGroups[keyZY].z = (int)hit.z;
        }
        positionGroups[keyZY].hitIndices.push_back(i);
    }
    
    // Convert map to vector and filter out single-hit groups
    for (const auto& pair : positionGroups) {
        if (pair.second.hitIndices.size() > 1) {
            finalGroups.push_back(pair.second);
        }
    }
    
    // Step 2: Process each group that contains at least one full (3-view) hit
    int groupsProcessed = 0;
    int hitsResolved = 0;
    
    for (size_t groupIdx = 0; groupIdx < finalGroups.size(); groupIdx++) {
        const auto& group = finalGroups[groupIdx];
        
        // Check if group contains at least one full hit
        bool hasFullHit = false;
        std::vector<int> fullHitIndices;
        std::vector<int> twoViewHitIndices;
        
        for (int hitIdx : group.hitIndices) {
            const Hit3D& hit = initialHits[hitIdx];
            if (hit.matchType == "Full" || hit.confidence >= 1.0) {
                hasFullHit = true;
                fullHitIndices.push_back(hitIdx);
            } else if (hit.nViews == 2) {
                twoViewHitIndices.push_back(hitIdx);
            }
        }
        
        if (!hasFullHit || twoViewHitIndices.empty()) {
            continue;  // Skip groups without full hits or without 2-view hits to resolve
        }
        
        groupsProcessed++;
        
        // Get timing from the first full hit
        const Hit3D& referenceHit = initialHits[fullHitIndices[0]];
        double sharedTiming = -1;
        if (group.sharedView == "XY") sharedTiming = referenceHit.timeXY;
        else if (group.sharedView == "XZ") sharedTiming = referenceHit.timeXZ;
        else if (group.sharedView == "ZY") sharedTiming = referenceHit.timeZY;
        
        // Process 2-view hits in this group
        for (int hitIdx : twoViewHitIndices) {
            Hit3D& hit = resolvedHits[hitIdx];
            
            // Check if the missing view matches the shared view
            bool missingXY = (hit.timeXY < 0);
            bool missingXZ = (hit.timeXZ < 0);
            bool missingZY = (hit.timeZY < 0);
            
            bool canResolve = false;
            if (group.sharedView == "XY" && missingXY) canResolve = true;
            else if (group.sharedView == "XZ" && missingXZ) canResolve = true;
            else if (group.sharedView == "ZY" && missingZY) canResolve = true;
            
            if (canResolve) {
                // Copy timing information
                if (group.sharedView == "XY") hit.timeXY = sharedTiming;
                else if (group.sharedView == "XZ") hit.timeXZ = sharedTiming;
                else if (group.sharedView == "ZY") hit.timeZY = sharedTiming;
                
                hit.nViews = 3;
                hit.matchType += "-Resolved";
                hitsResolved++;
            }
        }
        
        // Step 3: Redistribute charge for all hits in the group
        double totalSharedCharge = getViewCharge(group.sharedView, 
                                                resolvedHits[group.hitIndices[0]].x,
                                                resolvedHits[group.hitIndices[0]].y,
                                                resolvedHits[group.hitIndices[0]].z,
                                                viewXY, viewXZ, viewZY);
        
        if (totalSharedCharge <= 0) {
            continue;
        }
        
        // Store original charges before redistribution for all hits in group
        for (int hitIdx : group.hitIndices) {
            Hit3D& hit = resolvedHits[hitIdx];
            hit.redistributed = true;  // Mark this hit as having gone through redistribution
            
            // Store original charges for the shared view (before redistribution)
            if (group.sharedView == "XY" && hit.chargeXY > 0) {
                hit.originalChargeXY = hit.chargeXY;
            } else if (group.sharedView == "XZ" && hit.chargeXZ > 0) {
                hit.originalChargeXZ = hit.chargeXZ;
            } else if (group.sharedView == "ZY" && hit.chargeZY > 0) {
                hit.originalChargeZY = hit.chargeZY;
            }
        }
        
        // Calculate average charge from unshared views for each hit
        std::vector<double> unsharedAverages;
        double totalWeight = 0;
        
        for (int hitIdx : group.hitIndices) {
            const Hit3D& hit = resolvedHits[hitIdx];
            double unsharedSum = 0;
            int unsharedCount = 0;
            
            // Get charges from non-shared views
            if (group.sharedView != "XY" && hit.timeXY >= 0) {
                double charge = getViewCharge("XY", hit.x, hit.y, hit.z, viewXY, viewXZ, viewZY);
                if (charge > 0) {
                    unsharedSum += charge;
                    unsharedCount++;
                }
            }
            if (group.sharedView != "XZ" && hit.timeXZ >= 0) {
                double charge = getViewCharge("XZ", hit.x, hit.y, hit.z, viewXY, viewXZ, viewZY);
                if (charge > 0) {
                    unsharedSum += charge;
                    unsharedCount++;
                }
            }
            if (group.sharedView != "ZY" && hit.timeZY >= 0) {
                double charge = getViewCharge("ZY", hit.x, hit.y, hit.z, viewXY, viewXZ, viewZY);
                if (charge > 0) {
                    unsharedSum += charge;
                    unsharedCount++;
                }
            }
            
            double avgUnshared = (unsharedCount > 0) ? unsharedSum / unsharedCount : 1.0;
            unsharedAverages.push_back(avgUnshared);
            totalWeight += avgUnshared;
        }
        
        // Redistribute shared charge proportionally and update individual view charges
        if (totalWeight > 0) {
            for (size_t i = 0; i < group.hitIndices.size(); i++) {
                int hitIdx = group.hitIndices[i];
                double fraction = unsharedAverages[i] / totalWeight;
                double assignedCharge = totalSharedCharge * fraction;
                
                Hit3D& hit = resolvedHits[hitIdx];
                
                // Update the redistributed charge for the shared view
                if (group.sharedView == "XY") {
                    hit.chargeXY = assignedCharge;
                } else if (group.sharedView == "XZ") {
                    hit.chargeXZ = assignedCharge;
                } else if (group.sharedView == "ZY") {
                    hit.chargeZY = assignedCharge;
                }
                
                // Recalculate total charge for the hit using final redistributed charges
                double totalCharge = 0;
                int viewCount = 0;
                
                if (hit.timeXY >= 0) {
                    totalCharge += hit.chargeXY;
                    viewCount++;
                }
                if (hit.timeXZ >= 0) {
                    totalCharge += hit.chargeXZ;
                    viewCount++;
                }
                if (hit.timeZY >= 0) {
                    totalCharge += hit.chargeZY;
                    viewCount++;
                }
                
                hit.charge = (viewCount > 0) ? totalCharge : hit.charge;
            }
        }
    }
    
    // Step 4: Assign confidence levels based on timing correlation
    for (auto& hit : resolvedHits) {
        // Skip hits that already have perfect confidence
        if (hit.matchType == "Full" && hit.confidence >= 1.0) continue;
        
        // For resolved hits (now 3-view), calculate max time difference
        if (hit.matchType.find("Resolved") != std::string::npos && hit.nViews == 3) {
            std::vector<double> times;
            if (hit.timeXY >= 0) times.push_back(hit.timeXY);
            if (hit.timeXZ >= 0) times.push_back(hit.timeXZ);
            if (hit.timeZY >= 0) times.push_back(hit.timeZY);
            
            if (times.size() >= 2) {
                double maxTimeDiff = 0;
                for (size_t i = 0; i < times.size(); i++) {
                    for (size_t j = i + 1; j < times.size(); j++) {
                        maxTimeDiff = std::max(maxTimeDiff, std::abs(times[i] - times[j]));
                    }
                }
                
                // Assign confidence based on max time difference
                if (maxTimeDiff < 5.0) {
                    hit.confidence = 0.9;
                } else if (maxTimeDiff < 8.0) {
                    hit.confidence = 0.7;
                } else if (maxTimeDiff < 11.0) {
                    hit.confidence = 0.5;
                } else {
                    hit.confidence = 0.3;
                }
            }
        }
        // For unresolved 2-view hits, set confidence to 0.1
        else if (hit.nViews == 2 && hit.matchType.find("Resolved") == std::string::npos) {
            hit.confidence = 0.1;
        }
    }
    return resolvedHits;
}

// Entry point finding
Hit3D findEntryPoint(const std::vector<Hit3D>& allHits) {
    // ... [Include the full entry point finding algorithm]
    if (allHits.empty()) {
        Hit3D dummy = {};
        dummy.x = 0.0;
        dummy.y = 0.0;
        dummy.z = 0.0;
        dummy.timeXY = -1.0;
        dummy.timeXZ = -1.0;
        dummy.timeZY = -1.0;
        dummy.chargeXY = -1.0;
        dummy.chargeXZ = -1.0;
        dummy.chargeZY = -1.0;
        dummy.matchType = "Dummy";
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

// PCA analysis
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
    double residualPenalty = exp(-0.5 * pow((result.rmsResidual)/0.70, 2));  // approaches 1 for small residuals
    // double outlierRatio = (double)result.nOutliers / (double)result.nHighConfidenceHits;
    // double outlierPenalty = exp(-outlierRatio / 0.3);;  // exponential penalty for outliers
    
    result.confidenceScore = linearityScore * residualPenalty * (1.0 - 0.5 * planarityPenalty);
    
    result.valid = true;
    return result;
}

// Single track check
bool isSingleTrackEvent(const PCATrajectoryResult& pcaResult, 
                       double linearityThreshold = 0.95,
                       double planarityThreshold = 0.1,
                       double confidenceThreshold = 0.7) {
    if (!pcaResult.valid) return false;
    
    bool highLinearity = pcaResult.linearity > linearityThreshold;
    bool lowPlanarity = pcaResult.planarity < planarityThreshold;
    bool highConfidence = pcaResult.confidenceScore > confidenceThreshold;
    bool goodExplainedVariance = pcaResult.explainedVariance[0] > 0.9;
    
    return highLinearity && lowPlanarity && highConfidence && goodExplainedVariance;
}

// Trajectory fitting
TrajectoryParams fitTrajectory(const std::vector<Hit3D>& allHits) {
    TrajectoryParams params;
    params.valid = false;
    
    PCATrajectoryResult pcaResult = performPCA(allHits, 0.5);
    
    if (!pcaResult.valid) {
        return params;
    }
    
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

// Other utility functions
TString createOutputDirectory(const TString& inputFile, bool isSingleTrack) {
    TString inputDir = gSystem->DirName(inputFile);
    
    TString outputDir;
    if (isSingleTrack) {
        outputDir = inputDir + "/3D_Display_root_SingleTrack";
    } else {
        outputDir = inputDir + "/3D_Display_root_Other";
    }
    
    if (gSystem->AccessPathName(outputDir)) {
        if (gSystem->mkdir(outputDir, kTRUE) != 0) {
            std::cerr << "Warning: Could not create directory " << outputDir << std::endl;
            std::cerr << "Using current directory instead." << std::endl;
            if (isSingleTrack) {
                return "./3D_Display_root_SingleTrack";
            } else {
                return "./3D_Display_root_MultiTrack";
            }
        }
    }
    
    return outputDir;
}

void printSummary(const std::vector<Hit3D>& hits3D, int eventID) {
    int nFull = 0, nPartial = 0, nGroupResolved = 0;
    std::map<std::string, int> typeCounts;
    std::map<double, int> confidenceCounts;
    
    for (const auto& hit : hits3D) {
        if (hit.nViews == 3 && hit.confidence >= 1.0) nFull++;
        else if (hit.matchType.find("Resolved") != std::string::npos) nGroupResolved++;
        else nPartial++;
        
        typeCounts[hit.matchType]++;
        
        double roundedConf = round(hit.confidence / 0.1) * 0.1;
        confidenceCounts[roundedConf]++;
    }
    
    std::cout << "\nConfidence level distribution:\n";
    for (const auto& p : confidenceCounts) {
        std::cout << "    " << std::fixed << std::setprecision(1) << p.first << ": " << p.second << " hits" << std::endl;
    }
}

#endif // HIT3D_RECONSTRUCTION_CC