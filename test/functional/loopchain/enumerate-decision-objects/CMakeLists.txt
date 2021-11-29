list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-enumerate-decision-objects.cpp.in test-enumerate-decision-objects-${BACKEND}.cpp)

	raja_add_test( NAME test-enumerate-decision-objects-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-enumerate-decision-objects-${BACKEND}.cpp)

	target_include_directories(test-enumerate-decision-objects-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

