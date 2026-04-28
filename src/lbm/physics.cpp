#include <lbm/physics.hpp>

#include <cassert>
#include <cstdlib>

#include <omp.h>

#include <lbm/communications.hpp>
#include <lbm/config.hpp>
#include <lbm/structures.hpp>

#if DIRECTIONS == 9 && DIMENSIONS == 2
/// Definition of the 9 base vectors used to discretize the directions on each mesh.
const Vector direction_matrix[DIRECTIONS] = {
  // clang-format off
  {+0.0, +0.0},
  {+1.0, +0.0}, {+0.0, +1.0}, {-1.0, +0.0}, {+0.0, -1.0},
  {+1.0, +1.0}, {-1.0, +1.0}, {-1.0, -1.0}, {+1.0, -1.0},
  // clang-format on
};
#else
#error Need to define adapted direction matrix.
#endif

#if DIRECTIONS == 9
/// Weigths used to compensate the differences in lenght of the 9 directional vectors.
const double equil_weight[DIRECTIONS] = {
  // clang-format off
  4.0 / 9.0,
  1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
  1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
  // clang-format on
};

/// Opposite directions for bounce back implementation
const int opposite_of[DIRECTIONS] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
#else
#error Need to define adapted equilibrium distribution function
#endif

double get_vect_norm_2(const Vector a, const Vector b) {
  double res = 0.0;
  for (size_t k = 0; k < DIMENSIONS; k++) {
    res += a[k] * b[k];
  }
  return res;
}

double get_cell_density(const Mesh* mesh, int x, int y) {
  double res = 0.0;
  for (size_t k = 0; k < DIRECTIONS; k++) {
    res += Mesh_get_cell_const(mesh, x, y, k);
  }
  return res;
}

void get_cell_velocity(Vector v, const Mesh* mesh, int x, int y, double cell_density) {
  assert(v != NULL);

  // Loop on all dimensions
  for (size_t d = 0; d < DIMENSIONS; d++) {
    v[d] = 0.0;

    // Sum all directions
    for (size_t k = 0; k < DIRECTIONS; k++) {
      v[d] += Mesh_get_cell_const(mesh, x, y, k) * direction_matrix[k][d];
    }

    // Normalize
    v[d] /= cell_density;
  }
}

double compute_equilibrium_profile(Vector velocity, double density, int direction) {
  const double v2 = get_vect_norm_2(velocity, velocity);

  // Compute `e_i * v_i / c`
  const double p  = get_vect_norm_2(direction_matrix[direction], velocity);
  const double p2 = p * p;

  // Terms without density and direction weight
  double f_eq = 1.0 + (3.0 * p) + ((9.0 / 2.0) * p2) - ((3.0 / 2.0) * v2);

  // Multiply everything by the density and direction weight
  f_eq *= equil_weight[direction] * density;

  return f_eq;
}

static inline double compute_equilibrium_profile_opt(Vector velocity, double density, int direction, double v2) {
  const double p  = get_vect_norm_2(direction_matrix[direction], velocity);
  const double p2 = p * p;
  double f_eq = 1.0 + (3.0 * p) + ((9.0 / 2.0) * p2) - ((3.0 / 2.0) * v2);
  f_eq *= equil_weight[direction] * density;
  return f_eq;
}

void compute_cell_collision(Mesh* mesh_out, const Mesh* mesh_in, int x, int y) {
  // Compute macroscopic values
  const double density = get_cell_density(mesh_in, x, y);
  Vector v;
  get_cell_velocity(v, mesh_in, x, y, density);

  const double v2 = get_vect_norm_2(v, v);

  // Loop on microscopic directions
  for (size_t k = 0; k < DIRECTIONS; k++) {
    // call avec v2
    double f_eq = compute_equilibrium_profile_opt(v, density, k, v2);
    // Compute f_out
    Mesh_get_cell(mesh_out, x, y, k) = Mesh_get_cell_const(mesh_in, x, y, k) - RELAX_PARAMETER * (Mesh_get_cell_const(mesh_in, x, y, k) - f_eq);
  }
}

void compute_bounce_back(Mesh* mesh, int x, int y) {
  double tmp[DIRECTIONS];
  for (size_t k = 0; k < DIRECTIONS; k++) {
    tmp[k] = Mesh_get_cell(mesh, x, y, opposite_of[k]);
  }
  for (size_t k = 0; k < DIRECTIONS; k++) {
    Mesh_get_cell(mesh, x, y, k) = tmp[k];
  }
}

double helper_compute_poiseuille(const size_t i, const size_t size) {
  const double y = (double)(i - 1);
  const double L = (double)(size - 1);
  return 4.0 * INFLOW_MAX_VELOCITY / (L * L) * (L * y - y * y);
}

void compute_inflow_zou_he_poiseuille_distr(const Mesh* mesh_in, Mesh* mesh_out, int x, int y, size_t id_y) {
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif

  // Set macroscopic fluid info
  // Poiseuille distribution on X and null on Y
  // We just want the norm, so `v = v_x`
  const double v = helper_compute_poiseuille(id_y, mesh_in->height);

  // Compute rho from U and inner flow on surface
  const double rho = (Mesh_get_cell_const(mesh_in, x, y, 0) + Mesh_get_cell_const(mesh_in, x, y, 2) + Mesh_get_cell_const(mesh_in, x, y, 4) + 2 * (Mesh_get_cell_const(mesh_in, x, y, 3) + Mesh_get_cell_const(mesh_in, x, y, 6) + Mesh_get_cell_const(mesh_in, x, y, 7))) / (1.0 - v);

  // Now compute unknown microscopic values
  Mesh_get_cell(mesh_out, x, y, 1) = Mesh_get_cell_const(mesh_in, x, y, 3); // + (2.0/3.0) * density * v_y <--- no velocity on Y so v_y = 0
  Mesh_get_cell(mesh_out, x, y, 5) = Mesh_get_cell_const(mesh_in, x, y, 7) - (1.0 / 2.0) * (Mesh_get_cell_const(mesh_in, x, y, 2) - Mesh_get_cell_const(mesh_in, x, y, 4))
            + (1.0 / 6.0) * (rho * v); // + (1.0/2.0) * rho * v_y    <--- no velocity on Y so v_y = 0
  Mesh_get_cell(mesh_out, x, y, 8) = Mesh_get_cell_const(mesh_in, x, y, 6) + (1.0 / 2.0) * (Mesh_get_cell_const(mesh_in, x, y, 2) - Mesh_get_cell_const(mesh_in, x, y, 4))
            + (1.0 / 6.0) * (rho * v); //- (1.0/2.0) * rho * v_y    <--- no velocity on Y so v_y = 0

  // No need to copy already known one as the value will be "loss" in the wall at propagatation time
}

void compute_outflow_zou_he_const_density(Mesh* mesh, int x, int y) {
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif

  const double rho = 1.0;
  // Compute macroscopic velocity depending on inner flow going onto the wall
  const double v = -1.0 + (1.0 / rho) * (Mesh_get_cell(mesh, x, y, 0) + Mesh_get_cell(mesh, x, y, 2) + Mesh_get_cell(mesh, x, y, 4) + 2 * (Mesh_get_cell(mesh, x, y, 1) + Mesh_get_cell(mesh, x, y, 5) + Mesh_get_cell(mesh, x, y, 8)));

  // Now can compute unknown microscopic values
  Mesh_get_cell(mesh, x, y, 3) = Mesh_get_cell(mesh, x, y, 1) - (2.0 / 3.0) * rho * v;
  Mesh_get_cell(mesh, x, y, 7) = Mesh_get_cell(mesh, x, y, 5)
            + (1.0 / 2.0) * (Mesh_get_cell(mesh, x, y, 2) - Mesh_get_cell(mesh, x, y, 4))
            // - (1.0/2.0) * (rho * v_y)    <--- no velocity on Y so v_y = 0
            - (1.0 / 6.0) * (rho * v);
  Mesh_get_cell(mesh, x, y, 6) = Mesh_get_cell(mesh, x, y, 8)
            + (1.0 / 2.0) * (Mesh_get_cell(mesh, x, y, 4) - Mesh_get_cell(mesh, x, y, 2))
            // + (1.0/2.0) * (rho * v_y)    <--- no velocity on Y so v_y = 0
            - (1.0 / 6.0) * (rho * v);
}

void special_cells(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  // Loop on all inner cells
  for (size_t i = 1; i < mesh->width - 1; i++) {
    for (size_t j = 1; j < mesh->height - 1; j++) {
      switch (*(lbm_cell_type_t_get_cell(mesh_type, i, j))) {
      case CELL_FUILD:
        break;
      case CELL_BOUNCE_BACK:
        compute_bounce_back(mesh, i, j);
        break;
      case CELL_LEFT_IN:
        compute_inflow_zou_he_poiseuille_distr(mesh, mesh, i, j, j + mesh_comm->y);
        break;
      case CELL_RIGHT_OUT:
        compute_outflow_zou_he_const_density(mesh, i, j);
        break;
      }
    }
  }
}

void collision(Mesh* mesh_out, const Mesh* mesh_in) {
  assert(mesh_in->width == mesh_out->width);
  assert(mesh_in->height == mesh_out->height);

// Loop on all inner cells
#pragma omp parallel for schedule(static)
  for (size_t j = 1; j < mesh_in->height - 1; j++) {
    #pragma omp simd
    for (size_t i = 1; i < mesh_in->width - 1; i++) {
      compute_cell_collision(mesh_out, mesh_in, i, j);
    }
  }
}

void propagation(Mesh* mesh_out, const Mesh* mesh_in) {
// Loop on all cells
#pragma omp parallel for schedule(static)
  for (size_t j = 1; j < mesh_out->height - 1; j++) {
    #pragma omp simd
    for (size_t i = 1; i < mesh_out->width - 1; i++) {
      // For all direction
      for (size_t k = 0; k < DIRECTIONS; k++) {
        // Compute destination point
        const size_t ii = i + direction_matrix[k][0];
        const size_t jj = j + direction_matrix[k][1];
        // Propagate to neighboor nodes
        if ((ii >= 0 && ii < mesh_out->width) && (jj >= 0 && jj < mesh_out->height)) {
          Mesh_get_cell(mesh_out, ii, jj, k) = Mesh_get_cell_const(mesh_in, i, j, k);
        }
      }
    }
  }
}