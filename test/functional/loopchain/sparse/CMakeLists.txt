list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-sparse.cpp.in test-sparse-${BACKEND}.cpp)

	raja_add_test( NAME test-sparse-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-sparse-${BACKEND}.cpp)

	target_include_directories(test-sparse-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-sparse-demo.cpp.in test-sparse-demo-${BACKEND}.cpp)

	raja_add_test( NAME test-sparse-demo-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-sparse-demo-${BACKEND}.cpp)

	target_include_directories(test-sparse-demo-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

