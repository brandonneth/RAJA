list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-symbolic-segment.cpp.in test-symbolic-segment-${BACKEND}.cpp)

	raja_add_test( NAME test-symbolic-segment-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-symbolic-segment-${BACKEND}.cpp)

	target_include_directories(test-symbolic-segment-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)
