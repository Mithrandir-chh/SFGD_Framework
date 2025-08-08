#ifndef CALIBRATION_MANAGER_H
#define CALIBRATION_MANAGER_H

#include <map>
#include <string>
#include <utility>

// Forward declaration from the framework
class Hit;

class CalibrationManager {
public:
    CalibrationManager();
    ~CalibrationManager();
    
    // Load calibration data from file
    bool LoadCalibration(const std::string& filename);
    
    // Get gain for a specific view and coordinates
    double GetGain(const std::string& view, int coord1, int coord2);
    
    // Get calibrated charge from a hit object
    double GetCalibratedCharge(Hit* hit);
    
    // Check if calibration data is loaded
    bool IsCalibrationLoaded() const;
    
private:
    std::map<std::string, std::map<std::pair<int,int>, double>> gainMap;
    bool calibrationLoaded;
};

#endif // CALIBRATION_MANAGER_H