message(STATUS "building che workspace... ")
if    (ROOT_FOUND)
  ADD_EXECUTABLE(EventDisplayLite                   EventDisplayLite.cpp)
  TARGET_LINK_LIBRARIES(EventDisplayLite            SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(3D_EventTrajectory                   3D_EventTrajectory.cpp)
  TARGET_LINK_LIBRARIES(3D_EventTrajectory           SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(3D_EventTrajectorySide                   
    3D_EventTrajectorySide.cpp
    CalibrationManager.cpp
    )
  TARGET_LINK_LIBRARIES(3D_EventTrajectorySide           SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(Theta_phi                   Theta_phi.cpp)
  TARGET_LINK_LIBRARIES(Theta_phi           SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(SingleTrackFilter                   SingleTrackFilter.cpp)
  TARGET_LINK_LIBRARIES(SingleTrackFilter          SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(VoxelTrackLength                     VoxelTrackLength.cpp)
  TARGET_LINK_LIBRARIES(VoxelTrackLength             SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(GainCalculations                     GainCalculations.cpp)
  TARGET_LINK_LIBRARIES(GainCalculations             SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(CalibPlots                     CalibPlots.cpp)
  TARGET_LINK_LIBRARIES(CalibPlots             SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(MuonDecayEvents                     
    MuonDecayEvents.cpp
    CalibrationManager.cpp
    )
  TARGET_LINK_LIBRARIES(MuonDecayEvents             SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(3D_EventGamma                     
    3D_EventGamma.cpp
    CalibrationManager.cpp
    )
  TARGET_LINK_LIBRARIES(3D_EventGamma             SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(DecayMuonDistribution                     DecayMuonDistribution.cpp)
  TARGET_LINK_LIBRARIES(DecayMuonDistribution              SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(BetaDistribution                     BetaDistribution.cpp)
  TARGET_LINK_LIBRARIES(BetaDistribution              SuperFGD_Analysis ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(GammaDistribution                     GammaDistribution.cpp)
  TARGET_LINK_LIBRARIES(GammaDistribution              SuperFGD_Analysis ${ROOT_LIBRARIES})
endif (ROOT_FOUND)
