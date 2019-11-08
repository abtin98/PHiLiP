set(TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

configure_file(1d_collocation_test.prm 1d_collocation_test.prm COPYONLY)
add_test(
	NAME 1D_COLLOCATION_TEST
  COMMAND mpirun -n 1 ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_1D -i ${CMAKE_CURRENT_BINARY_DIR}/1d_collocation_test.prm
  WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(2d_collocation_test.prm 2d_collocation_test.prm COPYONLY)
add_test(
	NAME 2D_COLLOCATION_TEST
  COMMAND mpirun -n 1 ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_2D -i ${CMAKE_CURRENT_BINARY_DIR}/2d_collocation_test.prm
  WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(3d_collocation_test.prm 3d_collocation_test.prm COPYONLY)
add_test(
	NAME 3D_COLLOCATION_TEST
  COMMAND mpirun -n 1 ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/3d_collocation_test.prm
  WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

