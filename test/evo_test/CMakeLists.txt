CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(evo_test_source_files EvoTest.cpp)

link_directories(${MAINFOLDER}/lib)
include_directories(${MAINFOLDER}/thirdparty/include/boxmuller)

ADD_EXECUTABLE(EvoTest ${evo_test_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES})
TARGET_LINK_LIBRARIES(EvoTest ${project_LIBS} evo.so)