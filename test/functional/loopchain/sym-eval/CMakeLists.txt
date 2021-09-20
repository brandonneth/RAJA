list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-sym-eval.cpp.in test-sym-eval-${BACKEND}.cpp)

	raja_add_test( NAME test-sym-eval-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-sym-eval-${BACKEND}.cpp)

	target_include_directories(test-sym-eval-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

