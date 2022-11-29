list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-format-decisions.cpp.in test-format-decisions-${BACKEND}.cpp)

	raja_add_test( NAME test-format-decisions-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-format-decisions-${BACKEND}.cpp)

	target_include_directories(test-format-decisions-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

