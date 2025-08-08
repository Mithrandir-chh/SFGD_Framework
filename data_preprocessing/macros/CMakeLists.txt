message(STATUS "examples ... ")

ADD_EXECUTABLE(dump         DataDump.cpp)
TARGET_LINK_LIBRARIES(dump  SuperFGDunpack)

if    (ROOT_FOUND)
  ADD_EXECUTABLE(unpack                             Unpack.cpp)
  TARGET_LINK_LIBRARIES(unpack                      SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(ufe-monitor                        ufe_monitor.cpp)
  TARGET_LINK_LIBRARIES(ufe-monitor                 SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(TDMunpack                          TDMUnpack.cpp)
  TARGET_LINK_LIBRARIES(TDMunpack                   SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(GTrigTagSearch                     GTrigTagSearch.cpp)
  TARGET_LINK_LIBRARIES(GTrigTagSearch              SuperFGDunpack ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(HG_vs_LG_calib                     HG_LG_calib.cpp)
  TARGET_LINK_LIBRARIES(HG_vs_LG_calib              SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(HG_vs_ToT_calib                    HG_vs_ToT_calib.cpp)
  TARGET_LINK_LIBRARIES(HG_vs_ToT_calib             SuperFGDunpack ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(LG_vs_ToT_calib                    LG_vs_ToT_calib.cpp)
  TARGET_LINK_LIBRARIES(LG_vs_ToT_calib             SuperFGDunpack ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(Calibration	                    Calibration.cpp)
  TARGET_LINK_LIBRARIES(Calibration                 SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(EventStructure	        	    EventStructure.cpp)
  TARGET_LINK_LIBRARIES(EventStructure              SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(EventStructure_cosmics                     EventStructure_cosmics.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_cosmics              SuperFGDunpack ${ROOT_LIBRARIES})  

  ADD_EXECUTABLE(EventStructure_MC                     EventStructure_MC.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_MC              SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(EventStructure_MC_update                     EventStructure_MC_update.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_MC_update              SuperFGDunpack ${ROOT_LIBRARIES})  

  ADD_EXECUTABLE(EventStructure_MC2                     EventStructure_MC2.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_MC2              SuperFGDunpack ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(EventStructure_MC_pion                     EventStructure_MC_pion.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_MC_pion              SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(EventStructure_MC_proton                     EventStructure_MC_proton.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_MC_proton              SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(EventStructure_MC_highAngle                     EventStructure_MC_highAngle.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_MC_highAngle              SuperFGDunpack ${ROOT_LIBRARIES})
  
  ADD_EXECUTABLE(EventStructure_12free                     EventStructure_12free.cpp)
  TARGET_LINK_LIBRARIES(EventStructure_12free              SuperFGDunpack ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(auto_dq                        auto_dq.cpp)
  TARGET_LINK_LIBRARIES(auto_dq                 SuperFGDunpack ${ROOT_LIBRARIES})

endif (ROOT_FOUND)
