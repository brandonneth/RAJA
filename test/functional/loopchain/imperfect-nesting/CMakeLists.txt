list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-imperfect-nesting.cpp.in test-imperfect-nesting-${BACKEND}.cpp)

	raja_add_test( NAME test-imperfect-nesting-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-imperfect-nesting-${BACKEND}.cpp)

	target_include_directories(test-imperfect-nesting-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

