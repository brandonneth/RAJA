list(APPEND FORALL_BACKENDS Sequential)

if(RAJA_ENABLE_OPENMP)
  list(APPEND FORALL_BACKENDS OpenMP)
endif()

foreach( BACKEND ${FORALL_BACKENDS})

	configure_file(test-policy-to-nesting-order.cpp.in test-policy-to-nesting-order-${BACKEND}.cpp)

	raja_add_test( NAME test-policy-to-nesting-order-${BACKEND}
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-policy-to-nesting-order-${BACKEND}.cpp)

	target_include_directories(test-policy-to-nesting-order-${BACKEND}.exe
                           	   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
endforeach()


unset(FORALL_BACKENDS)

