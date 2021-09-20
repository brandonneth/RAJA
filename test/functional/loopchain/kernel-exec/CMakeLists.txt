list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-kernel-exec.cpp.in test-kernel-exec-${BACKEND}.cpp)

	raja_add_test( NAME test-kernel-exec-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-kernel-exec-${BACKEND}.cpp)

	target_include_directories(test-kernel-exec-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

