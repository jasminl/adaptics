CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(trackapp_source_files main.cpp Track.cpp TrackFeature.cpp TrackMatch.cpp TrackSAC.cpp
fastmap.h TrackMSTargetBase.cpp TrackMSTargetColor.cpp TrackMSTargetRGB.cpp
TrackMSTargetSIFT.h TrackFlow.cpp TrackParticleFilter.h TrackMeanShift.cpp)

#Add SYSTEM to disable warnings related to third party tnt library
include_directories(SYSTEM ${MAINFOLDER}/thirdparty/src/vlfeat) 

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(trackapp ${trackapp_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES} ${MAINFOLDER}/thirdparty/src/vlfeat/bin/g64/libvl.so)
TARGET_LINK_LIBRARIES(trackapp ${project_LIBS})