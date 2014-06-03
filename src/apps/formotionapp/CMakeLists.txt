CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(formotionapp_source_files main1.cpp xmlUtils.cpp DirManager.cpp)

#Add SYSTEM to disable warnings related to third party tnt library
include_directories(SYSTEM ${MAINFOLDER}/thirdparty/include 
${MAINFOLDER}/thirdparty/include/tnt)

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(formotionapp ${formotionapp_source_files})

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -g")

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES} )
TARGET_LINK_LIBRARIES(formotionapp ${MAINFOLDER}/thirdparty/lib/libtinyxml.so visio.so fftw3.so ${project_LIBS})