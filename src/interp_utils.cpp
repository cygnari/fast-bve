#include "general_utils.hpp"
#include "structs.hpp"
#include <cmath>
#include <iostream>
#include <tuple>

double sbb_coeff(const int deg, const int i, const int j) {
  std::vector<double> log_vals(deg + 1, 0);
  double accum = 0;
  for (int k = 1; k < deg + 1; k++) {
    log_vals[k] = log(k);
    accum += log_vals[k];
  }
  for (int k = 1; k < i + 1; k++) {
    accum -= log_vals[k];
  }
  for (int k = 1; k < j + 1; k++) {
    accum -= log_vals[k];
  }
  for (int k = 1; k < (deg-i-j) + 1; k++) {
    accum -= log_vals[k];
  }
  return exp(accum);
}

void fekete_init(std::vector<std::vector<double>> &points, const int degree) {
  double delta_x = 1.0 / degree;
  int index;
  double a, b, c, part;
  for (int i = 0; i < degree + 1; i++) {
    a = 1 - i * delta_x;
    a = 0.5 * (1 + sin(M_PI / 2 * (2 * a - 1)));
    part = a;
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      c = j * delta_x;
      b = i*delta_x - b;
      b = 0.5 * (1 + sin(M_PI / 2 * (2 * b - 1)));
      c = 0.5 * (1 + sin(M_PI / 2 * (2 * c - 1)));
      part += b + c;
      points[index][0] = a / part;
      points[index][1] = b / part;
      points[index][2] = c / part;
    }
  }
}

void interp_mat_init(
    std::vector<double> &mat, const std::vector<std::vector<double>> &points,
    const int degree,
    const int point_count) { // sets up matrix to interpolate with fekete points
  // simple polynomial in s and t, barycentric coordinates
  // for example, for deg 2, evaluates 1, s, t, s^2, st, t^2 at interpolation points
  int index, place;
  double a, b;
  for (int i = 0; i < degree + 1; i++) {
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      for (int k = 0; k < point_count; k++) {
        a = points[k][0];
        b = points[k][1];
        place = index * point_count + k;
        mat[place] = pow(a, i - j) * pow(b, j);
      }
    }
  }
}

void interp_mat_init_sbb(
    std::vector<double> &mat, const std::vector<std::vector<double>> &points,
    const int degree,
    const int point_count) { // sets up matrix to interpolate with fekete points
  // uses spherical bezier bernstein polynomials
  // for example, for deg 2, evaluates s^2, t^2, u^2, st, su, tu at interpolation points
  int index = 0, place;
  double s, t, u, si, tj, uk, comp, val, spart = 1, factor, tpart=1;
  for (int k = 0; k < point_count; k++) {
    s = points[k][0];
    t = points[k][1];
    u = points[k][2];
    index = 0;
    spart = 1;
    // factor = t/u;
    for (int i = 0; i < degree + 1; i++) {
      // val = spart * pow(u, degree-i);
      for (int j = 0; j < degree+1-i; j++) {
        val = spart * tpart * pow(u, degree-i-j);
        place = point_count * index + k;
        mat[place] = val;
        index++;
        val *= factor;
        tpart *= t;
      }
      spart *= s;
    }
  }
}

double interp_eval(
    const std::vector<double> &alphas, const double s, const double t,
    const int degree) { // evaluate interpolation polynomial with coefficients
                        // alpha and barycentric point (s, t)
  double accum = 0;
  int index;
  for (int i = 0; i < degree + 1; i++) {
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      accum += pow(s, i - j) * pow(t, j) * alphas[index];
    }
  }
  return accum;
}

std::vector<double> interp_vals_sbb(const double s, const double t, const double u, const int degree) {
  // returns vector of SBB basis values of s, t, u
  int count = (degree + 1) * (degree + 2) / 2;
  std::vector<double> out_vals (count, 0);
  double val, factor, spart = 1;
  factor = t / u;
  int index = 0;
  for (int i = 0; i < degree + 1; i++) { // degree of s
    val = spart * pow(u, degree - i);
    for (int j = 0; j < degree + 1 - i; j++) {
      out_vals[index] = val;
      index += 1;
      val *= factor;
    }
    spart *= s;
  }
  return out_vals;
}

double interp_eval_sbb(
    const std::vector<double> &alphas, const double s, const double t, const double u,
    const int degree) { // evaluate SBB interpolation polynomial with coefficients
                        // alpha and barycentric point (s, t, u)
  double accum = 0;
  int index = 0;
  double val, factor, spart = 1;
  factor = t / u;
  for (int i = 0; i < degree + 1; i++) { // degree of s
    val = spart * pow(u, degree - i);
    for (int j = 0; j < degree + 1 - i; j++) {
      accum += val * alphas[index];
      index += 1;
      val *= factor;
    }
    spart *= s;
  }
  return accum;
}

std::vector<double> bilinear_interp(const RunConfig &run_information,
                                    const std::vector<double> &target_point,
                                    const int iv1, const int iv2, const int iv3,
                                    const std::vector<double> &dynamics_state) {
  std::vector<double> v1, v2, v3, bary_cords, out;

  v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
  v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
  v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
  bary_cords = barycoords(v1, v2, v3, target_point);
  v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1,
             run_information.info_per_point);
  v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1,
             run_information.info_per_point);
  v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1,
             run_information.info_per_point);
  scalar_mult(v1, bary_cords[0]);
  scalar_mult(v2, bary_cords[1]);
  scalar_mult(v3, bary_cords[2]);
  out = v1;
  vec_add(out, v2);
  vec_add(out, v3);
  return out;
}

std::vector<double>
biquadratic_interp(const RunConfig &run_information,
                   const std::vector<double> &target_point, const int iv1,
                   const int iv2, const int iv3, const int iv4, const int iv5,
                   const int iv6, const std::vector<double> &dynamics_state) {

  std::vector<double> v1, v2, v3, v4, v5, v6, curr_alphas, bary_cords;
  std::vector<std::vector<double>> points(6, std::vector<double>(3, 0));
  std::vector<double> output_values(run_information.info_per_point);
  std::vector<double> interp_matrix(36, 0);
  std::vector<int> ipiv(6, 0);
  std::vector<double> b_vec(6*run_information.tracer_count+6, 0);

  output_values[0] = target_point[0];
  output_values[1] = target_point[1];
  output_values[2] = target_point[2];

  v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
  v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
  v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
  v4 = slice(dynamics_state, run_information.info_per_point * iv4, 1, 3);
  v5 = slice(dynamics_state, run_information.info_per_point * iv5, 1, 3);
  v6 = slice(dynamics_state, run_information.info_per_point * iv6, 1, 3);

  points[0][0] = 1;
  points[1][1] = 1;
  points[2][2] = 1;
  points[3] = normalized_barycoords(v1, v2, v3, v4);
  points[4] = normalized_barycoords(v1, v2, v3, v5);
  points[5] = normalized_barycoords(v1, v2, v3, v6);

  b_vec[0] = dynamics_state[run_information.info_per_point * iv1 + 3];
  b_vec[1] = dynamics_state[run_information.info_per_point * iv2 + 3];
  b_vec[2] = dynamics_state[run_information.info_per_point * iv3 + 3];
  b_vec[3] = dynamics_state[run_information.info_per_point * iv4 + 3];
  b_vec[4] = dynamics_state[run_information.info_per_point * iv5 + 3];
  b_vec[5] = dynamics_state[run_information.info_per_point * iv6 + 3];

  for (int j = 0; j < run_information.tracer_count; j++) {
    b_vec[6 * (j+1)] = dynamics_state[run_information.info_per_point * iv1 + 4 + j];
    b_vec[6 * (j+1) + 1] = dynamics_state[run_information.info_per_point * iv2 + 4 + j];
    b_vec[6 * (j+1) + 2] = dynamics_state[run_information.info_per_point * iv3 + 4 + j];
    b_vec[6 * (j+1) + 3] = dynamics_state[run_information.info_per_point * iv4 + 4 + j];
    b_vec[6 * (j+1) + 4] = dynamics_state[run_information.info_per_point * iv5 + 4 + j];
    b_vec[6 * (j+1) + 5] = dynamics_state[run_information.info_per_point * iv6 + 4 + j];
  }

  interp_mat_init(interp_matrix, points, 2, 6);
  int info = linear_solve(interp_matrix, b_vec, 6, 1+run_information.tracer_count, 1);
  if (info > 0)
    throw std::runtime_error("Biquadratic interpolation linear solve failed at line 136");
  bary_cords = barycoords(v1, v2, v3, target_point);
  curr_alphas = slice(b_vec, 0, 1, 6);
  output_values[3] = interp_eval(curr_alphas, bary_cords[0], bary_cords[1], 2);
  for (int j = 0; j < run_information.tracer_count; j++) {
    curr_alphas = slice(b_vec, 6 * j + 6, 1, 6);
    output_values[4 + j] = interp_eval(curr_alphas, bary_cords[0], bary_cords[1], 2);
  }
  return output_values;
}

void remesh_points(
    const RunConfig &run_information, std::vector<double> &target_points,
    const std::vector<double> &dynamics_state,
    const std::vector<std::vector<std::vector<int>>> &dynamics_triangles,
    const std::vector<std::vector<bool>> &dynamics_triangles_is_leaf,
    const int point_count, const double omega) {
  // remesh points back to regular point distribution
  std::vector<double> curr_target;
  int iv1, iv2, iv3, iv4, iv5, iv6, curr_level, tri_loc, super_tri_loc;
  double vor1, vor2, vor3, vor4, vor5, vor6, vormax, vormin, vor;
  for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
    if ((i < run_information.particle_lb) or (i >= run_information.particle_ub)) {
      curr_target.assign(run_information.info_per_point, 0);
    } else {
      curr_target =
          slice(target_points, run_information.info_per_point * i, 1, 3);
      std::tie(curr_level, tri_loc) = find_leaf_tri(
          curr_target, dynamics_state, dynamics_triangles,
          dynamics_triangles_is_leaf, run_information.info_per_point,
          run_information.dynamics_levels_max);
      super_tri_loc = floor(tri_loc / 4.0);

      iv1 = dynamics_triangles[curr_level - 1][super_tri_loc][0];
      iv2 = dynamics_triangles[curr_level - 1][super_tri_loc][1];
      iv3 = dynamics_triangles[curr_level - 1][super_tri_loc][2];
      iv4 = dynamics_triangles[curr_level][4 * super_tri_loc + 3][0];
      iv5 = dynamics_triangles[curr_level][4 * super_tri_loc + 3][1];
      iv6 = dynamics_triangles[curr_level][4 * super_tri_loc + 3][2];

      curr_target = biquadratic_interp(run_information, curr_target, iv1, iv2,
                                       iv3, iv4, iv5, iv6, dynamics_state);

      if (run_information.use_fast and (not run_information.fast_sum_rotate)) {
        // fast sum icos not rotated
        vor1 = dynamics_state[run_information.info_per_point * iv1 + 3] +
               2 * omega * dynamics_state[run_information.info_per_point * iv1 + 2];
        vor2 = dynamics_state[run_information.info_per_point * iv2 + 3] +
               2 * omega * dynamics_state[run_information.info_per_point * iv2 + 2];
        vor3 = dynamics_state[run_information.info_per_point * iv3 + 3] +
               2 * omega * dynamics_state[run_information.info_per_point * iv3 + 2];
        vor4 = dynamics_state[run_information.info_per_point * iv4 + 3] +
               2 * omega * dynamics_state[run_information.info_per_point * iv4 + 2];
        vor5 = dynamics_state[run_information.info_per_point * iv5 + 3] +
               2 * omega * dynamics_state[run_information.info_per_point * iv5 + 2];
        vor6 = dynamics_state[run_information.info_per_point * iv6 + 3] +
               2 * omega * dynamics_state[run_information.info_per_point * iv6 + 2];

        vormax = std::max(vor1, std::max(vor2,
                     std::max(vor3, std::max(vor4, std::max(vor5, vor6)))));
        vormin = std::min(vor1, std::min(vor2,
                     std::min(vor3, std::min(vor4, std::min(vor5, vor6)))));

        if (vormax > 0) { // some leeway
          vormax *= 1.1;
        } else {
          vormax *= 0.9;
        }
        if (vormin > 0) {
          vormin *= 0.9;
        } else {
          vormin *= 1.1;
        }

        vor = curr_target[3] + 2 * omega * curr_target[2];
        if ((vor > vormax) or (vor < vormin)) {
          // violate monotonicity, do bilinear interp
          curr_target =
              slice(target_points, run_information.info_per_point * i, 1, 3);
          iv1 = dynamics_triangles[curr_level][tri_loc][0];
          iv2 = dynamics_triangles[curr_level][tri_loc][1];
          iv3 = dynamics_triangles[curr_level][tri_loc][2];
          curr_target = bilinear_interp(run_information, curr_target, iv1, iv2,
                                        iv3, dynamics_state);
        }
      }
    }

    vector_copy(target_points, curr_target, run_information.info_per_point * i,
                run_information.info_per_point);
  }
}
