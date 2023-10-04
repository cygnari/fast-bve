#ifndef H_IO_UTILS_H
#define H_IO_UTILS_H

#include "structs.hpp"
#include <fstream>
#include <sstream>
#include <vector>

void write_state(const RunConfig &run_information,
                 const std::vector<double> &dynamics_state,
                 const std::vector<double> &dynamics_area,
                 std::ofstream &file_writer1, std::ofstream &file_writer2);

void write_triangles(
    const RunConfig &run_information,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    std::ofstream &file_writer3, std::ofstream &file_writer4);

std::string create_config(const RunConfig &run_information);

#endif
