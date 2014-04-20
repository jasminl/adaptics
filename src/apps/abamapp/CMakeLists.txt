CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(abamapp_source_files main.cpp)

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(abamapp ${abamapp_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES})
TARGET_LINK_LIBRARIES(abamapp ${project_LIBS} learn.so)