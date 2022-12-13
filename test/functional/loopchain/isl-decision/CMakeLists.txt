list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-isl-decision.cpp.in test-isl-decision-${BACKEND}.cpp)

	raja_add_test( NAME test-isl-decision-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-isl-decision-${BACKEND}.cpp)

	target_include_directories(test-isl-decision-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

