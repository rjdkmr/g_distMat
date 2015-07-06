cmake_minimum_required(VERSION 2.8)
project(g_distMat C)

set(g_distMat_MAJOR_VERSION 1)
set(g_distMat_g_distMat_MINOR_VERSION 1)
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

if(${GROMACS_VERSION} VERSION_LESS "4.5.0")
	message(FATAL_ERROR, "\n Not implemented for Gromacs version: ${GROMACS_VERSION} \n")
endif()

if(${GROMACS_VERSION} VERSION_EQUAL "5.0" OR ${GROMACS_VERSION} VERSION_GREATER "5.0")
	set(HAVE_GROMACS50 TRUE)
	add_definitions("-DHAVE_GROMACS50")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_EQUAL "5.0" OR ${GROMACS_VERSION} VERSION_GREATER "5.0")

include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )
##############################################################

add_executable(g_distMat g_distMat.c)
target_link_libraries(g_distMat ${GROMACS_LIBRARIES})
if(HAVE_GROMACS50)
  SET_TARGET_PROPERTIES(g_distMat PROPERTIES LINKER_LANGUAGE CXX)
endif(HAVE_GROMACS50)

install(TARGETS g_distMat DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
