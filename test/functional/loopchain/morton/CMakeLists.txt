list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-morton.cpp.in test-morton-${BACKEND}.cpp)

	raja_add_test( NAME test-morton-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-morton-${BACKEND}.cpp)

	target_include_directories(test-morton-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

