#ifndef H_STRUCTS_H
#define H_STRUCTS_H

#include <string>
#include <vector>

struct RunConfig {
  // primitive configuration options
  bool use_amr = false;
  bool use_remesh = false;
  bool use_fast = false;
  bool use_fixer = false;
  bool vor_fix = false;
  bool vor_limiter = false;
  bool bltc = false;
  bool same_source_targets = true;
  std::string out_path;    // ../run-output/ locally,
                           // on Derecho, /glade/derecho/scratch/achen/bve/
  int write_precision = 6; // number of decimal places, 6 for data
                           // visualization, 16 for error testing
  bool write_output = false;
  bool write_tris = false;
  bool write_stream = false; // compute and write stream function at each time step
  double radius = 1.0;
  double end_time;         // in days
  double delta_t;          // time step size
  int dynamics_levels_min; // initial icosahedron refinement levels, 0 = base icosahedron
  int dynamics_levels_max; // max icosahedron refinement if using amr
  std::string initial_vor_condition; // initial vorticity distribution
  std::string vor_forcing;           // vorticity forcing, if there is one
  int init_cond_param1;              // parameter for initial condition
  double init_cond_param2;           // parameter for initial condition
  // for RH, ICP1 is wavenumber, ICP2 is wave speed, default ICP1 = 4, ICP2 = 0
  // for GV, ICP1 is radius parameter, ICP2 is starting latitude * pi, default ICP1 = 4, ICP2 = 0.05
  // for RV, default ICP1 = 10, ICP2 = 0.05
  // for PV, no params
  int forcing_param1;    // parameter for forcing
  double forcing_param2; // parameter for forcing
  // for ssw forcing, FRP1 is forcing wavenumber, FRP2 is forcing duration [days], default FRP1 = 1, FRP2 = 11
  // for ssw_blend , FRP2 is forcing duration [days], default FRP2 = 11

  int interp_degree;      // interpolation degree
  int interp_point_count; // number of interpolation points
  int info_per_point; // how many doubles each point is, for example,
                      // storing x y z vor tracer = 5
  double amr_circ_thresh = 0.005; // threshold for circulation in amr
  double amr_vor_thresh = 0.4;    // threshold for vorticity difference in amr
  int amr_levels;                 // how many levels of AMR are permitted

  // fast sum info
  int fast_sum_cluster_thresh; // threshold for a triangle being a cluster
  int fast_sum_tree_levels;    // icosahedron levels for fast summation
  double fast_sum_theta;       // well separated threshold
  bool fast_sum_rotate = true; // whether or not to rotate the fast sum grid
  double fast_sum_rotate_alph = 0.01; // 3 rotation coefficients
  double fast_sum_rotate_beta = 0.01;
  double fast_sum_rotate_gamm = 0.01;

  // derived run config info
  int time_steps;                 // number of time steps
  int dynamics_initial_points;    // initial number of dynamics points
  int dynamics_max_points;        // maximum number of dynamics points
  int dynamics_initial_triangles; // initial number of dynamics triangles
  int dynamics_max_triangles;     // maximum number of dynamics triangles
  int dynamics_curr_point_count;  // current number of points
  int dynamics_curr_tri_count;    // current number of triangles
  int tracer_count;               // number of tracers
  int target_points;              // number of convolution target points

  // mpi info
  int mpi_P;       // total MPI ranks
  int mpi_ID;      // own MPI rank
  int particle_lb; // range of assigned particles
  int particle_ub;
  int particle_own;
  int target_lb;
  int target_ub;
  int target_own;
};

struct IcosTree {
  int tree_depth;
  std::vector<std::vector<double>> icosahedron_vertex_coords;
  std::vector<std::vector<std::vector<double>>> icosahedron_tri_centers;
  std::vector<std::vector<double>> icosahedron_tri_radii;
  std::vector<std::vector<std::vector<int>>>
      icosahedron_triangle_vertex_indices;
};

struct InteractionPair {
  int lev_target; // icosahedron level of target/source
  int lev_source;
  int curr_target; // index of target/source
  int curr_source;
  int count_target; // number of particles in target/source
  int count_source;
  int type; // pp, pc, cp, or cc
            // 0 for pp, 1, for pc, 2 for cp, 3 for cc
};

#endif
