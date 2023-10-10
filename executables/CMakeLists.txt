add_executable(driver driver.cpp)
add_executable(single_rhs single_rhs.cpp)
target_link_libraries(driver PUBLIC fastbve)
target_include_directories(driver PUBLIC INTERFACE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
           $<INSTALL_INTERFACE:include>
           $<INSTALL_INTERFACE:include>
)
target_link_libraries(single_rhs PUBLIC fastbve)
target_include_directories(single_rhs PUBLIC INTERFACE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
           $<INSTALL_INTERFACE:include>
           $<INSTALL_INTERFACE:include>
)
set(FastBVE_install_bin driver single_rhs PARENT_SCOPE)
install(TARGETS driver single_rhs
        EXPORT FastBVETargets
        DESTINATION bin
)

