list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-format-converter.cpp.in test-format-converter-${BACKEND}.cpp)

	raja_add_test( NAME test-format-converter-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-format-converter-${BACKEND}.cpp)

	target_include_directories(test-format-converter-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

