add_executable(driver driver.cpp)
add_executable(single_rhs single_rhs.cpp)
target_link_libraries(driver PUBLIC SrcLib)
target_link_libraries(single_rhs PUBLIC SrcLib)
