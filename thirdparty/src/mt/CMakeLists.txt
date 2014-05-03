CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

include_directories(${MAINFOLDER}/thirdparty/include/mt)

SET(mt_source_files mt19937-3.c
${MAINFOLDER}/thirdparty/include/mt/mt19937-3.h)

link_directories(${MAINFOLDER}/thirdparty/lib)

ADD_LIBRARY(mt SHARED ${mt_source_files})
SET(LIBRARY_OUTPUT_PATH ${MAINFOLDER}/thirdparty/lib)
