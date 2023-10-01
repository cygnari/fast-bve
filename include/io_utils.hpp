#ifndef H_IO_UTILS_H
#define H_IO_UTILS_H

#include <fstream>
#include <sstream>
#include "structs.hpp"
#include <vector>

using namespace std;

void write_state(const run_config& run_information, const vector<double>& dynamics_state, const vector<double>& dynamics_area, ofstream& file_writer1, ofstream& file_writer2);

void write_triangles(const run_config& run_information, const vector<vector<vector<int>>>& dynamics_triangles,
            const vector<vector<bool>>& dynamics_triangles_is_leaf, ofstream& file_writer3, ofstream& file_writer4);

string create_config(const run_config& run_information);

#endif
