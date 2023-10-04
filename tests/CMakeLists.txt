add_executable(CountNanTest countnantest.cpp)
target_link_libraries(CountNanTest PUBLIC SrcLib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(CountNanTest)
