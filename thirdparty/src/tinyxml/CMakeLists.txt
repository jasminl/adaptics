CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (Adaptics)

include_directories(${MAINFOLDER}/thirdparty/include/tinyxml )

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

SET(txml_source_files tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp 
	${MAINFOLDER}/thirdparty/include/tinyxml/tinyxml.h)

link_directories(${MAINFOLDER}/thirdparty/lib)

ADD_LIBRARY(tinyxml SHARED ${txml_source_files})
SET(LIBRARY_OUTPUT_PATH ${MAINFOLDER}/thirdparty/lib)
