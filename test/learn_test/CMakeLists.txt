CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

SET(learn_test_source_files LearnTest.cpp)

link_directories(${MAINFOLDER}/lib)

ADD_EXECUTABLE(learntest ${learn_test_source_files})

SET (project_LIBS ${Boost_LIBRARIES} ${Log4Cxx_LIBRARIES})
TARGET_LINK_LIBRARIES(learntest ${project_LIBS} learn.so)