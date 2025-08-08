message(STATUS "src ... ")

AUX_SOURCE_DIRECTORY(.            SOURCE_FILES)
#AUX_SOURCE_DIRECTORY(../../analysis/src/classes/       SOURCE_FILES)

set (EVENT_DICT_HEADERFILES
    ../../analysis/src/classes/Event.hh
    ../../analysis/src/classes/Hit.hh
)

set( EVENT_HEADERFILES
    ${EVENT_DICT_HEADERFILES}
)

ROOT_GENERATE_DICTIONARY(EventDict ${EVENT_DICT_HEADERFILES} LINKDEF LinkDef.h )

ADD_LIBRARY(SuperFGDunpack      SHARED      ${SOURCE_FILES}  ../../analysis/src/classes/Event.cc ../../analysis/src/classes/Hit.cc EventDict.cxx )
target_link_libraries(SuperFGDunpack ${ROOT_LIBRARIES})


add_custom_command(
    TARGET SuperFGDunpack
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install
    )  



