add_library(fastbve
    # add all FastBVE library source/impl files here
    amr.cpp
    conservation_fixer.cpp
    direct_sum_utils.cpp
    fast_sum_utils.cpp
    general_utils.cpp
    green_funcs.cpp
    init_utils.cpp
    input_utils.cpp
    interp_utils.cpp
    io_utils.cpp
    mpi_utils.cpp
    rhs_utils.cpp
    vorticity_functions.cpp)

target_link_libraries(fastbve PUBLIC MPI::MPI_C LAPACK::LAPACK)
target_include_directories(fastbve PUBLIC
           $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
           $<INSTALL_INTERFACE:include>)

install(FILES
    # add all FastBVE header/include files here
    ${CMAKE_BINARY_DIR}/fastbve-config.h
    amr.hpp
    conservation_fixer.hpp
    direct_sum_utils.hpp
    fast_sum_utils.hpp
    general_utils.hpp
    green_funcs.hpp
    init_utils.hpp
    input_utils.hpp
    interp_utils.hpp
    io_utils.hpp
    mpi_utils.hpp
    rhs_utils.hpp
    structs.hpp
    vorticity_functions.hpp
    DESTINATION include)

install(TARGETS fastbve
        EXPORT FastBVETargets
        DESTINATION lib)

