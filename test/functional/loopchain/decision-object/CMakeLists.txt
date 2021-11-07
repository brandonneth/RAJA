list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-decision-object.cpp.in test-decision-object-${BACKEND}.cpp)

	raja_add_test( NAME test-decision-object-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-decision-object-${BACKEND}.cpp)

	target_include_directories(test-decision-object-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

