add_executable(rz-test-foo foo.cpp foo.h)
target_link_libraries(rz-test-foo gtest gtest_main)

add_executable(rz-test-bar bar.cpp)
target_link_libraries(rz-test-bar gtest gtest_main)

#add_executable(rz-test rz-testmain.cpp)
#target_link_libraries(rz-test gtest gtest_main
#        rz-test-foo
#        rz-test-bar)
#include(GoogleTest)
gtest_add_tests(TARGET rz-test-foo)
gtest_add_tests(TARGET rz-test-bar)
#include(CTest)

#add_test(NAME hui COMMAND rz-test-foo && rz-test-bar)
