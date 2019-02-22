add_executable(test-ebasi ebasi.cpp)
target_link_libraries(test-ebasi gtest gtest_main)
gtest_add_tests(TARGET test-ebasi)