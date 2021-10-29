list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-auto-layout.cpp.in test-auto-layout-${BACKEND}.cpp)

	raja_add_test( NAME test-auto-layout-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-auto-layout-${BACKEND}.cpp)

	target_include_directories(test-auto-layout-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

