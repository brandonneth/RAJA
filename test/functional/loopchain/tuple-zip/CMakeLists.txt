list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-tuple-zip.cpp.in test-tuple-zip-${BACKEND}.cpp)

	raja_add_test( NAME test-tuple-zip-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-tuple-zip-${BACKEND}.cpp)

	target_include_directories(test-tuple-zip-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

