add_executable(test-foo foo.cpp foo.h)
target_link_libraries(test-foo gtest gtest_main)

add_executable(test-bar bar.cpp)
target_link_libraries(test-bar gtest gtest_main)

gtest_add_tests(TARGET test-foo)
gtest_add_tests(TARGET test-bar)

