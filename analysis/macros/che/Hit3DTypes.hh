#ifndef HIT3D_TYPES_H
#define HIT3D_TYPES_H

#include <vector>
#include <string>
#include <utility>

// Structure to hold 3D hit information
struct Hit3D {
    double x, y, z;
    double charge;
    double timeXY, timeXZ, timeZY;
    // Individual view charges (new addition)
    double chargeXY, chargeXZ, chargeZY;
    // Store original charges before redistribution for overlapped hits
    double originalChargeXY, originalChargeXZ, originalChargeZY;
    int nViews;
    std::string matchType;
    double confidence;  // confidence level 0.0 to 1.0
    bool redistributed;  // Flag to track if this hit went through redistribution
};

// Structure to organize hits by view
struct ViewHits {
    std::vector<int> indices;
    std::vector<std::pair<int, int>> positions;
    std::vector<double> charges;
    std::vector<double> times;
};

// Structure for trajectory parameters
struct TrajectoryParams {
    double x0, y0, z0;  // Entry point
    double vx, vy, vz;  // Direction vector
    double theta, phi;  // Angles
    bool valid;
};

struct PCATrajectoryResult {
    double x0, y0, z0, vx, vy, vz, theta, phi;
    bool valid;
    
    // New PCA data
    double centroid_x, centroid_y, centroid_z;
    double pc1_x, pc1_y, pc1_z, pc2_x, pc2_y, pc2_z, pc3_x, pc3_y, pc3_z;
    double eigenvalues[3], explainedVariance[3], cumulativeVariance[3];
    double linearity, planarity, trackWidth, trackThickness, totalVariance;
    double rmsResidual, maxResidual, confidenceScore;
    int nHighConfidenceHits;
    double weightedCentroidDistance;
    double meanLongitudinalDistance, meanTransverseDistance, meanNormalDistance;
    int nOutliers;
    double outlierThreshold;
};

#endif // HIT3D_TYPES_H