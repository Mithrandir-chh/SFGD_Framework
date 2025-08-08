#include "Hit3DReconstructor.hh"
#include <iostream>
#include <map>
#include <algorithm>
#include <limits>

// Static helper functions implementation
template<typename Derived>
double Hit3DReconstructorCRTP<Derived>::GetViewCharge(const std::string& view, double x, double y, double z, 
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

template<typename Derived>
double Hit3DReconstructorCRTP<Derived>::GetViewTiming(const std::string& view, double x, double y, double z, 
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

template<typename Derived>
Hit3D Hit3DReconstructorCRTP<Derived>::FindEntryPoint(const std::vector<Hit3D>& allHits) {
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

// Placeholder implementations for the complex algorithm functions
// TODO: Move your actual algorithm implementations here

template<typename Derived>
std::vector<Hit3D> Hit3DReconstructorCRTP<Derived>::FindThreeViewMatches(const ViewHits& viewXY, 
                                                                         const ViewHits& viewXZ, 
                                                                         const ViewHits& viewZY,
                                                                         std::vector<bool>& usedXY, 
                                                                         std::vector<bool>& usedXZ, 
                                                                         std::vector<bool>& usedZY) {
    std::vector<Hit3D> hits3D;
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
                    // Use bin centers, not edges (add 0.5 to get center of unit bin)
                    hit3d.x = x_xy + 0.5;
                    hit3d.y = y_xy + 0.5;
                    hit3d.z = z_xz + 0.5;
                    hit3d.charge = (viewXY.charges[ixy] + viewXZ.charges[ixz] + viewZY.charges[izy]);
                    hit3d.timeXY = viewXY.times[ixy];
                    hit3d.timeXZ = viewXZ.times[ixz];
                    hit3d.timeZY = viewZY.times[izy];
                    // Initialize individual view charges
                    hit3d.chargeXY = viewXY.charges[ixy];
                    hit3d.chargeXZ = viewXZ.charges[ixz];
                    hit3d.chargeZY = viewZY.charges[izy];
                    // Initialize original charges (same as current for full hits initially)
                    hit3d.originalChargeXY = viewXY.charges[ixy];
                    hit3d.originalChargeXZ = viewXZ.charges[ixz];
                    hit3d.originalChargeZY = viewZY.charges[izy];
                    hit3d.nViews = 3;
                    hit3d.matchType = "Full";
                    hit3d.confidence = 1.0;  // Highest confidence - unambiguous 3-view matches
                    hit3d.redistributed = false;  // Initially not redistributed
                    
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
    return hits3D;
}

template<typename Derived>
std::vector<Hit3D> Hit3DReconstructorCRTP<Derived>::FindTwoViewMatches(const ViewHits& viewXY, 
                                                                       const ViewHits& viewXZ, 
                                                                       const ViewHits& viewZY,
                                                                       std::vector<bool>& usedXY, 
                                                                       std::vector<bool>& usedXZ, 
                                                                       std::vector<bool>& usedZY) {
    std::vector<Hit3D> hits3D;
    
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
    return hits3D;
}

template<typename Derived>
std::vector<Hit3D> Hit3DReconstructorCRTP<Derived>::ResolveAdjacentHits(const std::vector<Hit3D>& initialHits, 
                                                                        const ViewHits& viewXY, 
                                                                        const ViewHits& viewXZ, 
                                                                        const ViewHits& viewZY) {
    std::vector<Hit3D> resolvedHits = initialHits;  // Start with initial hits
    
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
        double totalSharedCharge = GetViewCharge(group.sharedView, 
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
                double charge = GetViewCharge("XY", hit.x, hit.y, hit.z, viewXY, viewXZ, viewZY);
                if (charge > 0) {
                    unsharedSum += charge;
                    unsharedCount++;
                }
            }
            if (group.sharedView != "XZ" && hit.timeXZ >= 0) {
                double charge = GetViewCharge("XZ", hit.x, hit.y, hit.z, viewXY, viewXZ, viewZY);
                if (charge > 0) {
                    unsharedSum += charge;
                    unsharedCount++;
                }
            }
            if (group.sharedView != "ZY" && hit.timeZY >= 0) {
                double charge = GetViewCharge("ZY", hit.x, hit.y, hit.z, viewXY, viewXZ, viewZY);
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
    
    return resolvedHits;  // Placeholder return
}