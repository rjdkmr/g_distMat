cmake_minimum_required(VERSION 2.8)
project(g_distMat C)

set(g_distMat_MAJOR_VERSION 1)
set(g_distMat_g_distMat_MINOR_VERSION 0)
set(g_distMat_PATCH_VERSION 0)
set(g_distMat_VERSION
  ${g_distMat_MAJOR_VERSION}.${g_distMat_g_distMat_MINOR_VERSION}.${g_distMat_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
	set(DEF_INSTALL_CMAKE_DIR CMake)
else()
set(DEF_INSTALL_CMAKE_DIR lib/CMake/g_distMat)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_${p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

##############################################################
################### CHECKING GROMACS PATH ####################
##############################################################
find_package(GROMACS COMPONENTS "libgromacs")
if (NOT GROMACS_FOUND)
	message(STATUS "Continue checking for other Gromacs version library file")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	find_package(GROMACS COMPONENTS "libgmx")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs\n  or  \n export CMAKE_PREFIX_PATH=\"/path/to/gromacs\"\n" )
endif(NOT GROMACS_FOUND)

message(STATUS "Found Gromacs include directory: " ${GROMACS_INCLUDE_DIRS})

if(${GROMACS_VERSION} VERSION_GREATER "4.5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.0")
	set(HAVE_GROMACS45 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS45")
endif()

if(${GROMACS_VERSION} VERSION_GREATER "5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.1")
	set(HAVE_GROMACS50 TRUE)
	add_definitions("-DHAVE_GROMACS50")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_GREATER "5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.1")

if(${GROMACS_VERSION} VERSION_GREATER "5.1" AND ${GROMACS_VERSION} VERSION_LESS "5.2" )
	set(HAVE_GROMACS51 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS51")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_GREATER "5.1" AND ${GROMACS_VERSION} VERSION_LESS "5.2")

if(${GROMACS_VERSION} VERSION_GREATER "2016" AND ${GROMACS_VERSION} VERSION_LESS "2017" )
	set(HAVE_GROMACS2016 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS2016")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_GREATER "2016" AND ${GROMACS_VERSION} VERSION_LESS "2017")

include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )
##############################################################

##############################################################
# After gromacs2016, it uses c++11, and therefore it must be compiled with c++11 compiler
if( HAVE_GROMACS2016 )
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
	if(COMPILER_SUPPORTS_CXX11)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	elseif(COMPILER_SUPPORTS_CXX0X)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	else()
		message(FATAL_ERROR "This version of GROMACS requires C++11-compatible standard library. Please use a newer compiler, or a newer standard library, or use the GROMACS 5.1.x release. See the installation guide for details.")
	endif()

endif()
##############################################################

if( HAVE_GROMACS45 )
  add_executable(g_distMat g_distMat.c)
  target_link_libraries(g_distMat ${GROMACS_LIBRARIES} )
endif()

if(HAVE_GROMACS50)
  add_executable(g_distMat g_distMat.c)
  target_link_libraries(g_distMat ${GROMACS_LIBRARIES} )
  SET_TARGET_PROPERTIES(g_distMat PROPERTIES LINKER_LANGUAGE CXX)
endif(HAVE_GROMACS50)

if(HAVE_GROMACS51)
  add_executable(g_distMat g_distMat.c)
  target_link_libraries(g_distMat ${GROMACS_LIBRARIES} )
  SET_TARGET_PROPERTIES(g_distMat PROPERTIES LINKER_LANGUAGE CXX)
endif(HAVE_GROMACS51)

if(HAVE_GROMACS2016)
  add_executable(g_distMat g_distMat_2016.cpp)
  target_link_libraries(g_distMat ${GROMACS_LIBRARIES} )
  SET_TARGET_PROPERTIES(g_distMat PROPERTIES LINKER_LANGUAGE CXX)
endif(HAVE_GROMACS2016)

install(TARGETS g_distMat DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
