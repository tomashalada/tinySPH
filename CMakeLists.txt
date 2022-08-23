cmake_minimum_required(VERSION 3.2)

project(CompactNSearch)

# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (UNIX)
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)
endif (UNIX)

OPTION(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" OFF)
if (BUILD_AS_SHARED_LIBS)
	if(WIN32)
		add_definitions(-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")
	endif(WIN32)
endif (BUILD_AS_SHARED_LIBS)

SET(CMAKE_DEBUG_POSTFIX "_d")

set (HEADER_FILES
		libs/CompactNSearch/include/Config.h
		libs/CompactNSearch/include/CompactNSearch.h
		libs/CompactNSearch/include/PointSet.h
		libs/CompactNSearch/include/DataStructures.h)

if (BUILD_AS_SHARED_LIBS)
	add_library(CompactNSearch SHARED
	    ${HEADER_FILES}
		libs/CompactNSearch/src/CompactNSearch.cpp
)
else()
	add_library(CompactNSearch
	    ${HEADER_FILES}
		libs/CompactNSearch/src/CompactNSearch.cpp
)
endif ()

target_include_directories(CompactNSearch PUBLIC libs/CompactNSearch/include)

install(FILES "include/CompactNSearch" ${HEADER_FILES}
	DESTINATION libs/CompactNSearch/include/)

install(TARGETS CompactNSearch
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

OPTION(USE_DOUBLE_PRECISION "Use double precision"	ON)
if (USE_DOUBLE_PRECISION)
	add_definitions( -DUSE_DOUBLE)
endif (USE_DOUBLE_PRECISION)

### Testcases
# Dambreak
option(BUILD_DAMBREAK2D_DBC "Build example of how to use this library."
		ON)
if(BUILD_DAMBREAK2D_DBC)
	add_subdirectory(cases/dambreak2D_withDBC dambreak2D_withDBC)
endif(BUILD_DAMBREAK2D_DBC)

option(BUILD_DAMBREAK2D_DBCsimple "Build example of how to use this library."
		ON)
if(BUILD_DAMBREAK2D_DBCsimple)
	add_subdirectory(cases/dambreak2D_withDBCsimple dambreak2D_withDBCsimple)
endif(BUILD_DAMBREAK2D_DBCsimple)

option(BUILD_DAMBREAK2D_BI "Build example of how to use this library."
		ON)
if(BUILD_DAMBREAK2D_BI)
	add_subdirectory(cases/dambreak2D_withBI dambreak2D_withBI)
endif(BUILD_DAMBREAK2D_BI)

option(BUILD_DAMBREAK2D_MDBC "Build example of how to use this library."
		ON)
if(BUILD_DAMBREAK2D_MDBC)
    add_subdirectory(cases/dambreak2D_withMDBC dambreak2D_withMDBC)
endif(BUILD_DAMBREAK2D_MDBC)

# Stillwater
option(BUILD_STILLWATER2D_DBC "Build example of how to use this library."
		ON)
if(BUILD_STILLWATER2D_DBC)
    add_subdirectory(cases/stillwater2D_withDBC stillwater2D_withDBC)
endif(BUILD_STILLWATER2D_DBC)

option(BUILD_STILLWATER2D_MDBC "Build example of how to use this library."
		ON)
if(BUILD_STILLWATER2D_MDBC)
    add_subdirectory(cases/stillwater2D_withMDBC stillwater2D_withMDBC)
endif(BUILD_STILLWATER2D_MDBC)

option(BUILD_STILLWATER2D_MDBC_PERIODIC "Build example of how to use this library."
		ON)
if(BUILD_STILLWATER2D_MDBC_PERIODIC)
    add_subdirectory(cases/stillwater2D_withMDBC_periodic stillwater2D_withMDBC_periodic)
endif(BUILD_STILLWATER2D_MDBC_PERIODIC)
