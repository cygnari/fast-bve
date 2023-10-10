add_executable(general_util_tests general_util_tests.cpp)
target_link_libraries(general_util_tests PUBLIC SrcLib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(general_util_tests)
