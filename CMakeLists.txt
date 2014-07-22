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

# set up include-directories
include_directories(
	"${g_distMat_SOURCE_DIR}"   # to find src/*.h
	"${g_distMat_BINARY_DIR}")  # to find src/config.h	 

# Checking GROMACS PAth
if (NOT DEFINED GMX_PATH AND (NOT DEFINED GMX_LIB OR NOT DEFINED GMX_INCLUDE) )
message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs \nOR\n -DGMX_LIB=/path/to/gromacs/lib -DGMX_INCLUDE=/path/to/gromacs/include")
endif()


if(DEFINED GMX_PATH)
	set(GMX_INCLUDE "${GMX_PATH}/include")
	set(GMX_LIB "${GMX_PATH}/lib")
endif()

include_directories("${GMX_INCLUDE}/gromacs")

# Searching FFTW library
find_library(LIBFFTW fftw3f HINTS ${FFTW_LIB})
if(DEFINED LIBFFTW-NOTFOUND)
message(FATAL_ERROR "\nFFTW library file libfftw3f.so or libfftw3f.a not found at standard location\nUSE:\n -DFFTW_LIB=/path/to/fftw3/lib")
endif()

# Searching GROMACS Library
find_library(LIBGMX gmx HINTS ${GMX_LIB})
find_library(LIBGMXANA gmxana HINTS ${GMX_LIB})
find_library(LIBMD md HINTS ${GMX_LIB})

add_executable(g_distMat g_distMat.c)


target_link_libraries(g_distMat m ${LIBFFTW} ${LIBGMX} ${LIBGMXANA} ${LIBMD} pthread dl)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

install(TARGETS g_distMat DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)