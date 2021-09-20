raja_add_test( NAME test-forall-exec
               SOURCES ${CMAKE_CURRENT_BINARY_DIR}/test-forall-exec.cpp)

target_include_directories(test-forall-exec.exe
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)

