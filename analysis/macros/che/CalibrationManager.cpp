#include "CalibrationManager.hh"
#include <iostream>
#include <fstream>

CalibrationManager::CalibrationManager() : calibrationLoaded(false) {
}

CalibrationManager::~CalibrationManager() {
}

bool CalibrationManager::LoadCalibration(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error: Cannot open calibration file " << filename << std::endl;
        return false;
    }
    
    gainMap.clear();
    
    std::string view;
    int coord1, coord2;
    double gain;
    
    while (file >> view >> coord1 >> coord2 >> gain) {
        // Skip invalid gains (marked as -1.000000)
        if (gain > 0) {
            gainMap[view][std::make_pair(coord1, coord2)] = gain;
        }
    }
    file.close();
    
    calibrationLoaded = true;
    // std::cout << "Loaded calibration for " << gainMap.size() << " views" << std::endl;
    return true;
}

double CalibrationManager::GetGain(const std::string& view, int coord1, int coord2) {
    auto viewIt = gainMap.find(view);
    if (viewIt != gainMap.end()) {
        auto coordIt = viewIt->second.find(std::make_pair(coord1, coord2));
        if (coordIt != viewIt->second.end()) {
            return coordIt->second;
        }
    }
    return -1.0; // Return -1 if calibration not found
}

bool CalibrationManager::IsCalibrationLoaded() const {
    return calibrationLoaded;
}