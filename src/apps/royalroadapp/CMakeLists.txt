CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(royalroadapp_source_files main.cpp)

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(royalroadapp ${royalroadapp_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES})
TARGET_LINK_LIBRARIES(royalroadapp ${project_LIBS} evo.so)