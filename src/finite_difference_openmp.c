/* atlc - arbitrary transmission line calculator.

   Multi-core finite-difference solver, using red-black (checkerboard)
   Gauss-Seidel SOR parallelised with OpenMP.

   This replaces the old (disabled, incorrect) pthread solver in
   finite_difference_multi_threaded.c. The per-node update is the SAME function
   (new_node_value() in relax_node.h) that the single-threaded 4-direction
   sweep uses, so the physics is identical; only the sweep ORDER differs.

   Red-black is correct here because every node's update reads only its four
   orthogonal neighbours, which on a (i+j)%2 checkerboard are all the opposite
   colour. Updating one colour at a time therefore has no data dependencies
   within the colour, so the result is independent of the number of threads.

   If atlc is built without OpenMP support, finite_difference_openmp() simply
   falls back to the single-threaded solver, and finite_difference() always
   calls the single-threaded solver. */

#include "config.h"

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#include "definitions.h"
#include "exit_codes.h"
#include "relax_node.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern int width, height;
extern double **Vij, **Er;
extern unsigned char **oddity;
extern int dielectrics_to_consider_just_now;
extern double r;
extern int coupler;
extern int number_of_workers;

double finite_difference_openmp()
{
#ifdef _OPENMP
  int i, j, sweep, colour;
  double g, energy_per_metre, capacitance_per_metre;

  /* Same SOR weight selection as update_voltage_array(). */
  if (dielectrics_to_consider_just_now == 1)
    g = r;
  else
    g = 1;

  /* ITERATIONS full red-black sweeps. One red pass + one black pass updates
     every node once, matching the 25*4 = 100 node updates the single-threaded
     solver performs per call. */
  for (sweep = 0; sweep < ITERATIONS; sweep++) {
    for (colour = 0; colour < 2; colour++) {
#pragma omp parallel for schedule(static) private(j) default(shared)
      for (i = 0; i < width; i++)
        for (j = ((i + colour) & 1); j < height; j += 2)
          Vij[i][j] = new_node_value(i, j, g, Vij);
    }
  }

  /* Energy integral - identical loop and range to the single-threaded solver
     (find_energy_per_metre() is a pure read-only function, so the reduction is
     safe). The summation order varies with thread count, giving differences at
     the ~1e-15 level only. */
  energy_per_metre = 0.0;
#pragma omp parallel for schedule(static) private(j) reduction(+:energy_per_metre)
  for (i = 0; i < width; i++)
    for (j = 0; j < height; j++)
      energy_per_metre += find_energy_per_metre(i, j);

  if (coupler == FALSE)
    capacitance_per_metre = 2 * energy_per_metre;
  else
    capacitance_per_metre = energy_per_metre;
  return capacitance_per_metre;
#else
  /* Built without OpenMP: behave like the single-threaded solver. */
  return finite_difference_single_threaded();
#endif
}

/* Solver dispatcher used throughout do_fd_calculation.c.
   number_of_workers == 0 forces the original single-threaded algorithm (kept
   as the benchmark/baseline path); any other value uses the OpenMP solver,
   which honours the thread count set by omp_set_num_threads() / OMP_NUM_THREADS
   (defaulting to all available cores when -t is not given). */
double finite_difference()
{
#ifdef _OPENMP
  if (number_of_workers == 0)
    return finite_difference_single_threaded();
  return finite_difference_openmp();
#else
  return finite_difference_single_threaded();
#endif
}
