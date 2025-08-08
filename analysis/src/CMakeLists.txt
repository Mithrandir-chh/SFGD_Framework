message(STATUS "src ...")

AUX_SOURCE_DIRECTORY(./classes                 SOURCE_FILES)

# REQUIRED TO GIVE RELATIVE PATHS TO ROOT_GENERATE
set (EVENT_DICT_HEADERFILES
    classes/Event.hh
    classes/Hit.hh
    classes/ND280SFGDEvent.hh
    classes/ND280SFGDHit.hh
    classes/ND280SFGDTrack.hh
    classes/ND280SFGDVoxel.hh
    classes/ND280SFGDVoxelSet.hh
)

set( EVENT_HEADERFILES
    ${EVENT_DICT_HEADERFILES}
)

##################################################
# Dictionary

set( EVENT_LINKDEF_HEADERFILE LinkDef.h )
set( EVENT_SOURCEFILES ${EVENT_SOURCEFILES} ${CMAKE_CURRENT_BINARY_DIR}/EventDict.cxx )
set( EVENT_DICT_PCMFILE ${CMAKE_CURRENT_BINARY_DIR}/EventDict_rdict.pcm )

ROOT_GENERATE_DICTIONARY(  EventDict ${EVENT_DICT_HEADERFILES} LINKDEF ${EVENT_LINKDEF_HEADERFILE} )



ADD_LIBRARY(SuperFGD_Analysis      SHARED      ${SOURCE_FILES}  EventDict.cxx )
target_link_libraries(SuperFGD_Analysis ${ROOT_LIBRARIES})


add_custom_command(
    TARGET SuperFGD_Analysis
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install
    )
