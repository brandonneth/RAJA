list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-apply-overlapped-tile.cpp.in test-apply-overlapped-tile-${BACKEND}.cpp)

	raja_add_test( NAME test-apply-overlapped-tile-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-apply-overlapped-tile-${BACKEND}.cpp)

	target_include_directories(test-apply-overlapped-tile-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

