cmake_minimum_required(VERSION 3.6)
project(CoreModule2_PeptideFiltering_Windows-x86_64)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)

SET (BOOST_ROOT "/home/miin/boost_1_60_0/boost")
SET (BOOST_INCLUDEDIR "/home/miin/boost_1_60_0/boost/include")
SET (BOOST_LIBRARYDIR "/home/miin/boost_1_60_0/boost/lib")


set(SOURCE_FILES main.cpp)

file(GLOB SOURCES
        *.h
        *.cpp)

#set(Boost_USE_STATIC_LIBS   ON)
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBRARIES OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(CoreModule2_PeptideFiltering_Windows-x86_64 ${SOURCE_FILES} ${SOURCES})
