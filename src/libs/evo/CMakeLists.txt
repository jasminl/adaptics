CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(evo_source_files CFRR2.cpp CRR2.cpp Population.cpp
RoyalRoad.cpp Stats.cpp ESDna.h GrayCode.cpp
Seq.cpp
)

link_directories(${MAINFOLDER}/lib)
SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES})

include_directories(${MAINFOLDER}/thirdparty/include/boxmuller)

ADD_LIBRARY(evo SHARED ${evo_source_files})
TARGET_LINK_LIBRARIES(evo ${project_LIBS})
#SET_TARGET_PROPERTIES(evo PROPERTIES VERSION "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}" OUTPUT_NAME ${project_BIN} CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS evo DESTINATION lib)

