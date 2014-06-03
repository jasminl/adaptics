CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(visio_source_files
Layer.cpp KernelLayer.cpp Buf.cpp
Level2.cpp Level3.cpp Level4.cpp Level5A.cpp 
Level5B.cpp Level6.cpp Parameters.cpp)

link_directories(${MAINFOLDER}/lib)
SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES} fftw3.so)
include_directories(${MAINFOLDER}/thirdparty/include ${MAINFOLDER}/thirdparty/include/tnt)

#Disable unknown warnings in TNT/JAMA libs, add support for openmp
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -fopenmp")

ADD_LIBRARY(visio SHARED ${visio_source_files})
TARGET_LINK_LIBRARIES(visio ${project_LIBS})
INSTALL(TARGETS visio DESTINATION lib)

