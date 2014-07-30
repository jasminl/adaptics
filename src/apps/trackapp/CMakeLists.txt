CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(trackapp_source_files main.cpp track.cpp trackFeature.cpp trackMatch.cpp trackSAC.cpp
fastmap.h TrackMSTargetBase.cpp trackMSTargetColor.cpp trackMSTargetRGB.cpp
trackMSTargetSIFT.h trackFlow.cpp trackParticleFilter.h trackMeanShift.cpp)

#Add SYSTEM to disable warnings related to third party tnt library
include_directories(SYSTEM ${MAINFOLDER}/thirdparty/src/vlfeat) 

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(trackapp ${trackapp_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES} ${MAINFOLDER}/thirdparty/src/vlfeat/bin/g64/libvl.so)
TARGET_LINK_LIBRARIES(trackapp ${project_LIBS})