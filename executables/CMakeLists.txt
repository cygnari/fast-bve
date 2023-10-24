add_executable(driver driver.cpp)
target_link_libraries(driver PUBLIC fastbve)
target_include_directories(driver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(single_rhs single_rhs.cpp)
target_link_libraries(single_rhs PUBLIC fastbve)
target_include_directories(single_rhs PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

add_executable(lpm_interface_example lpm_interface_example.cpp)
target_link_libraries(lpm_interface_example PUBLIC fastbve)
target_include_directories(lpm_interface_example PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

if (BUILD_BARYTREE)
  add_executable(bltc_rhs bltc_rhs.cpp)
  target_link_libraries(bltc_rhs PUBLIC fastbve)
  target_include_directories(bltc_rhs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>)
  target_include_directories(bltc_rhs PUBLIC "../../BVE-BaryTree/src/interface/")
  target_link_directories(bltc_rhs PUBLIC "../../BVE-BaryTree/build/lib/")
  target_link_libraries(bltc_rhs PUBLIC BaryTree_cpu)
endif()

install(TARGETS driver single_rhs lpm_interface_example
        EXPORT FastBVETargets
        DESTINATION ${CMAKE_INSTALL_BINDIR})
