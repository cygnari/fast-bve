add_executable(general_util_tests general_util_tests.cpp)
target_link_libraries(general_util_tests PUBLIC fastbve GTest::gtest_main)
target_include_directories(general_util_tests PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

include(GoogleTest)
gtest_discover_tests(general_util_tests)
