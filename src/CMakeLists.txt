add_library(SrcLib amr.cpp
    conservation_fixer.cpp
    direct_sum_utils.cpp
    fast_sum_utils.cpp
    general_utils.cpp
    green_funcs.cpp
    init_utils.cpp
    input_utils.cpp
    interp_utils.cpp
    io_utils.cpp
    lpm_fast_sum.cpp
    mpi_utils.cpp
    rhs_utils.cpp
    vorticity_functions.cpp)

target_link_libraries(SrcLib PUBLIC MPI::MPI_C LAPACK::LAPACK)

install(FILES amr.hpp
    conservation_fixer.hpp
    direct_sum_utils.hpp
    fast_sum_utils.hpp
    general_utils.hpp
    green_funcs.hpp
    init_utils.hpp
    input_utils.hpp
    interp_utils.hpp
    io_utils.hpp
    lpm_fast_sum.hpp
    mpi_utils.hpp
    rhs_utils.hpp
    structs.hpp
    vorticity_functions.hpp
    DESTINATION ${PROJECT_SOURCE_DIR}/include)
