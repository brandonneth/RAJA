list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-apply-shift-and-fuse.cpp.in test-apply-shift-and-fuse-${BACKEND}.cpp)

	raja_add_test( NAME test-apply-shift-and-fuse-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-apply-shift-and-fuse-${BACKEND}.cpp)

	target_include_directories(test-apply-shift-and-fuse-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

