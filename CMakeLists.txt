cmake_minimum_required(VERSION 3.16)

# set the project name
project(pdb-redo VERSION 1.1.0 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

enable_testing()

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

# Build shared libraries by default (not my cup of tea, but hey)
option(BUILD_SHARED_LIBS "Build a shared library instead of a static one" ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# set(CMAKE_DEBUG_POSTFIX d)

if(MSVC)
    # make msvc standards compliant...
    add_compile_options(/permissive-)

	macro(get_WIN32_WINNT version)
		if (WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})

	# On Windows, do not install in the system location
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(WARNING "The library and auxiliary files will be installed in $ENV{LOCALAPPDATA}/${PROJECT_NAME}")
		set(CMAKE_INSTALL_PREFIX "$ENV{LOCALAPPDATA}/${PROJECT_NAME}" CACHE PATH "..." FORCE)
	endif()

	# Find out the processor type for the target
	if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
		set(COFF_TYPE "x64")
	elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i386")
		set(COFF_TYPE "x86")
	elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ARM64")
		set(COFF_TYPE "arm64")
	else()
		message(FATAL_ERROR "Unsupported or unknown processor type ${CMAKE_SYSTEM_PROCESSOR}")
	endif()	

	set(COFF_SPEC "--coff=${COFF_TYPE}")
endif()

if(UNIX AND NOT APPLE)
	# On Linux, install in the $HOME/.local folder by default
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(WARNING "The library and auxiliary files will be installed in $ENV{HOME}/.local")
		set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "..." FORCE)
	endif()
endif()

# Optionally use mrc to create resources
find_program(MRC mrc HINTS "$ENV{LOCALAPPDATA}/mrc" "$ENV{HOME}/.local/bin" "${CMAKE_INSTALL_PREFIX}/../mrc" "/usr/local/bin")

if(MRC)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(WARNING "Not using resources since mrc was not found")
endif()

if(USE_RSRC STREQUAL "ON")
	set(USE_RSRC 1)

	message("Using resources compiled with ${MRC}")
	add_compile_definitions(USE_RSRC)
endif()

# Libraries

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

set (Boost_DETAILED_FAILURE_MSG ON)
find_package(Boost 1.71.0 REQUIRED COMPONENTS system iostreams regex date_time program_options)

find_package(ZLIB)
find_package(BZip2)
find_package(GSL)
find_package(zeep REQUIRED)
find_package(cifpp 1.1 REQUIRED)

find_library(NEWUOA_LIBRARY newuoa REQUIRED)
find_path(NEWUOA_INCLUDE_DIR newuoa.h)
check_include_file( "newuoa.h" HAVE_NEWUOA )

find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-ccp4 clipper-mmdb
                           clipper-contrib clipper-minimol)

include_directories(${Boost_INCLUDE_DIR} ${NEWUOA_INCLUDE_DIR} ${ZEEP_INCLUDE_DIR} ${CIFPP_INCLUDE_DIR})
link_libraries(${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CCP4_LIBRARIES} cifpp::cifpp zeep::zeep ${NEWUOA_LIBRARY})

# Create a revision file, containing the current git version info

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
	include(GetGitRevisionDescription)
    get_git_head_revision(REFSPEC COMMITHASH)

	# Generate our own version string
	git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)
else()
    message(WARNING "Git not found, cannot set version info")

    SET(BUILD_VERSION_STRING "unknown")
endif()

# generate version.h
include_directories(${CMAKE_BINARY_DIR} PRIVATE)
string(TIMESTAMP BUILD_DATE_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
configure_file("${CMAKE_SOURCE_DIR}/src/revision.hpp.in" "${CMAKE_BINARY_DIR}/revision.hpp" @ONLY)

# Sources

set(project_sources 
${PROJECT_SOURCE_DIR}/src/AtomShape.cpp
${PROJECT_SOURCE_DIR}/src/ClipperWrapper.cpp
${PROJECT_SOURCE_DIR}/src/DistanceMap.cpp
${PROJECT_SOURCE_DIR}/src/MapMaker.cpp
${PROJECT_SOURCE_DIR}/src/Ramachandran.cpp
${PROJECT_SOURCE_DIR}/src/ResolutionCalculator.cpp
${PROJECT_SOURCE_DIR}/src/SkipList.cpp
${PROJECT_SOURCE_DIR}/src/Statistics.cpp
${PROJECT_SOURCE_DIR}/src/Symmetry-2.cpp
)

set(project_headers 
${PROJECT_SOURCE_DIR}/include/pdb-redo/AtomShape.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/ClipperWrapper.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/DistanceMap.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/MapMaker.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/Ramachandran.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/ResolutionCalculator.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/SkipList.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/Statistics.hpp
${PROJECT_SOURCE_DIR}/include/pdb-redo/Symmetry-2.hpp
)

add_library(pdb-redo ${project_sources} ${project_headers})
set_target_properties(pdb-redo PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(pdb-redo
	PUBLIC
	"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_link_options(pdb-redo PRIVATE -undefined dynamic_lookup)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR} )
set(LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR} )
set(SHARE_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/libpdb-redo)

# Fix cache dir
add_compile_definitions(
	# CACHE_DIR="${CMAKE_INSTALL_PREFIX}/${SHARE_INSTALL_DIR}"
	DATA_DIR="${CMAKE_INSTALL_PREFIX}/${SHARE_INSTALL_DIR}" )

generate_export_header(pdb-redo
	EXPORT_FILE_NAME pdb-redo/pdb-redo_export.hpp)

# Install rules

install(TARGETS pdb-redo
	EXPORT pdb-redoTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT pdb-redoTargets
	FILE "pdb-redoTargets.cmake"
	NAMESPACE pdb-redo::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pdb-redo
)

install(
	DIRECTORY include/pdb-redo
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
)

install(
	FILES "${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redo_export.hpp"
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pdb-redo
	COMPONENT Devel
)

configure_package_config_file(Config.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pdb-redo
	PATH_VARS INCLUDE_INSTALL_DIR LIBRARY_INSTALL_DIR SHARE_INSTALL_DIR
)

install(FILES
		"${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfig.cmake"
		"${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfigVersion.cmake"
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pdb-redo
	COMPONENT Devel
)

set(pdb-redo_MAJOR_VERSION 1)
set_target_properties(pdb-redo PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION 1
	INTERFACE_pdb-redo_MAJOR_VERSION 1)

set_property(TARGET pdb-redo APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING pdb-redo_MAJOR_VERSION
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/pdb-redo/pdb-redoConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Unit tests

option(PDB_REDO_BUILD_TESTS "Build test exectuables" OFF)

if(PDB_REDO_BUILD_TESTS)

	if(USE_RSRC)
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pdb_redo_test_rsrc.obj
			COMMAND ${MRC} -o ${CMAKE_CURRENT_BINARY_DIR}/pdb_redo_test_rsrc.obj ${CIFPP_SHARE_DIR}/mmcif_pdbx_v50.dic ${COFF_SPEC}
		)
		set(PDB_REDO_TEST_RESOURCE ${CMAKE_CURRENT_BINARY_DIR}/pdb_redo_test_rsrc.obj)
	endif()

	list(APPEND PDB_REDO_tests unit)

	foreach(PDB_REDO_TEST IN LISTS PDB_REDO_tests)
		set(PDB_REDO_TEST "${PDB_REDO_TEST}-test")
		set(PDB_REDO_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/test/${PDB_REDO_TEST}.cpp")

		add_executable(${PDB_REDO_TEST} ${PDB_REDO_TEST_SOURCE} ${PDB_REDO_TEST_RESOURCE})

		target_include_directories(${PDB_REDO_TEST} PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${CMAKE_CURRENT_BINARY_DIR}  # for config.h
		)

		target_link_libraries(${PDB_REDO_TEST} Threads::Threads ${Boost_LIBRARIES} pdb-redo cifpp::cifpp)

		if(${ZLIB_FOUND})
			target_link_libraries(${PDB_REDO_TEST} ZLIB::ZLIB)
		endif()

		if(${BZip2_FOUND})
			target_link_libraries(${PDB_REDO_TEST} BZip2::BZip2)
		endif()
	
		if(MSVC)
			# Specify unwind semantics so that MSVC knowns how to handle exceptions
			target_compile_options(${PDB_REDO_TEST} PRIVATE /EHsc)
		endif()

		add_custom_target("run-${PDB_REDO_TEST}" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${PDB_REDO_TEST}.touch ${PDB_REDO_TEST})

		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${PDB_REDO_TEST}.touch
			COMMAND $<TARGET_FILE:${PDB_REDO_TEST}>
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)

		add_test(NAME ${PDB_REDO_TEST}
			COMMAND $<TARGET_FILE:${PDB_REDO_TEST}>
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)

	endforeach()
endif()

message("Still to do: libpdb-redo.pc, package, tools, create symop_Data")


# include_directories(
#     ${PROJECT_SOURCE_DIR}/include
#     ${PROJECT_SOURCE_DIR}/src
# )

# add_library(${PROJECT_NAME} ${project_sources} ${project_headers})

# configure_file(${PROJECT_SOURCE_DIR}/src/config.hpp.cmakein ${CMAKE_BINARY_DIR}/config.hpp)
# target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")
# install(FILES ${CMAKE_BINARY_DIR}/config.hpp DESTINATION ${INCLUDE_INSTALL_DIR}/pdb-redo/)

# if(BUILD_SHARED_LIBS)
# set(prefix      ${CMAKE_INSTALL_PREFIX})
# set(exec_prefix ${CMAKE_INSTALL_PREFIX})
# set(libdir      ${CMAKE_INSTALL_PREFIX}/lib)
# set(includedir  ${CMAKE_INSTALL_PREFIX}/include)
# configure_file(${PROJECT_SOURCE_DIR}/libpdb-redo.pc.cmakein
#    ${CMAKE_BINARY_DIR}/libpdb-redo.pc.in @ONLY)
# file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libpdb-redo.pc
#    INPUT ${CMAKE_BINARY_DIR}/libpdb-redo.pc.in)
# install(FILES ${CMAKE_BINARY_DIR}/libpdb-redo.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)
# endif(BUILD_SHARED_LIBS)

# target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")

# configure_file(${PROJECT_SOURCE_DIR}/include/pdb-redo/revision.hpp.in revision.hpp @ONLY)
# install(FILES ${CMAKE_BINARY_DIR}/revision.hpp DESTINATION ${INCLUDE_INSTALL_DIR}/pdb-redo/)

# #if(Boost_FOUND)
#     include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
#     target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
# #endif()
# #if(CCP4_FOUND)
#     include_directories(${PROJECT_NAME} PUBLIC ${CCP4_INCLUDE_DIRS})
#     target_link_libraries(${PROJECT_NAME} ${CCP4_LIBRARIES})
# #endif()
# #if(ZEEP_FOUND)
#     target_link_libraries(${PROJECT_NAME} ${ZEEP_LIBRARY})
#     include_directories(${PROJECT_NAME} PUBLIC ${ZEEP_INCLUDE_DIR})
# #endif()
# #if(NEWUOA_FOUND)
#     target_link_libraries(${PROJECT_NAME} ${NEWUOA_LIBRARY})
#     include_directories(${PROJECT_NAME} PUBLIC ${NEWUOA_INCLUDE_DIR})
# #endif()
# #if(NEWUOA_FOUND)
#     target_link_libraries(${PROJECT_NAME} ${CIFPP_LIBRARY})
#     include_directories(${PROJECT_NAME} PUBLIC ${CIFPP_INCLUDE_DIR})
# #endif()
# #if(GSL_FOUND)
#     target_link_libraries(${PROJECT_NAME} GSL::gsl GSL::gslcblas)
#     include_directories(${PROJECT_NAME} PUBLIC GSL_INCLUDE_DIRS)
# #endif()

# install(DIRECTORY include/pdb-redo DESTINATION ${INCLUDE_INSTALL_DIR})

# install(TARGETS ${PROJECT_NAME}
#     LIBRARY DESTINATION ${LIB_INSTALL_DIR}
#     ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
#     RUNTIME DESTINATION ${BIN_INSTALL_DIR}
# )

# target_compile_definitions(${PROJECT_NAME}
# 	PUBLIC "" CACHE_DIR="${SHARE_INSTALL_DIR}/libpdb-redo" DATA_DIR="${SHARE_INSTALL_DIR}/libpdb-redo" )

# if(MSVC)
#     # make msvc standards compliant...
#     target_compile_options(${PROJECT_NAME} PRIVATE /permissive-)
# endif()

