CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

include_directories(${MAINFOLDER}/thirdparty/include/boxmuller)

SET(bm_source_files boxmuller.c 
	${MAINFOLDER}/thirdparty/include/boxmuller/boxmuller.h)

link_directories(${MAINFOLDER}/thirdparty/lib)

ADD_LIBRARY(boxmuller SHARED ${bm_source_files})
SET(LIBRARY_OUTPUT_PATH ${MAINFOLDER}/thirdparty/lib)
