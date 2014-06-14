CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(graycodeapp_source_files main.cpp)

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(graycodeapp ${graycodeapp_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES} libevo.so)
TARGET_LINK_LIBRARIES(graycodeapp ${project_LIBS})