find_package(Catch2 3 QUIET)

if(NOT Catch2_FOUND)
	include(FetchContent)

	FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG GIT_TAG v3.5.2)

	FetchContent_MakeAvailable(Catch2)
endif()

list(APPEND PDB_REDO_tests pdb-redo-unit rsr restr2ccd)

foreach(PDB_REDO_TEST IN LISTS PDB_REDO_tests)
	set(PDB_REDO_TEST "${PDB_REDO_TEST}-test")
	set(PDB_REDO_TEST_SOURCE "${PROJECT_SOURCE_DIR}/test/${PDB_REDO_TEST}.cpp")

	add_executable(${PDB_REDO_TEST} ${PDB_REDO_TEST_SOURCE})
	target_link_libraries(${PDB_REDO_TEST} PRIVATE pdb-redo Catch2::Catch2)

	if(MSVC)
		# Specify unwind semantics so that MSVC knowns how to handle exceptions
		target_compile_options(${PDB_REDO_TEST} PRIVATE /EHsc)
	endif()

	add_custom_target("run-${PDB_REDO_TEST}" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${PDB_REDO_TEST}.touch ${PDB_REDO_TEST})

	add_custom_command(
		OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${PDB_REDO_TEST}.touch
		COMMAND $<TARGET_FILE:${PDB_REDO_TEST}> -- ${PROJECT_SOURCE_DIR}/test)

	add_test(NAME ${PDB_REDO_TEST}
		COMMAND $<TARGET_FILE:${PDB_REDO_TEST}> --data-dir ${PROJECT_SOURCE_DIR}/test)
endforeach()
