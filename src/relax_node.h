/* relax_node.h - the per-node SOR update, shared by both finite-difference
   solvers:

     * update_voltage_array.c          (single-threaded, 4-direction sweep)
     * finite_difference_openmp.c      (multi-core, red-black sweep)

   The body below is a verbatim transcription of the per-pixel logic that used
   to live inline in update_voltage_array.c, so the numerical result of a node
   update is IDENTICAL between the two solvers - only the order in which nodes
   are swept differs (4-direction Gauss-Seidel vs red-black Gauss-Seidel). Both
   orderings converge to the same discrete solution.

   Every branch reads only the four orthogonal neighbours (i+-1,j)/(i,j+-1)
   (corner/edge branches read fixed boundary nodes). On a (i+j)%2 checkerboard
   those neighbours are always the opposite colour, which is what makes the
   red-black decomposition correct and thread-count independent.

   Includers must already have included "config.h", <stdio.h>, <stdlib.h>,
   "definitions.h" and "exit_codes.h" (for the oddity constants, fprintf/exit
   and INTERNAL_ERROR). */

#ifndef RELAX_NODE_H
#define RELAX_NODE_H

extern int width, height;
extern double **Er;
extern unsigned char **oddity;
extern int dielectrics_to_consider_just_now;

/* Return the new voltage for node (i,j), reading neighbour voltages from V.
   g is the SOR weight (r when relaxing a single dielectric, 1 otherwise).
   Conductors return their fixed potential; an unrecognised non-conductor node
   is left unchanged (matching the original fall-through behaviour). */
static inline double new_node_value(int i, int j, double g, double **V)
{
  unsigned char oddity_value = oddity[i][j];
  double Va, Vb, Vl, Vr, ERa, ERb, ERl, ERr, Vnew;

  if (oddity_value == CONDUCTOR_ZERO_V)
    return 0.0;

  else if (oddity_value == CONDUCTOR_PLUS_ONE_V)
    return 1.0;

  else if (oddity_value == CONDUCTOR_MINUS_ONE_V)
    return -1.0;

  else if (oddity_value == TOP_LEFT_CORNER) {
    Vnew = 0.5 * (V[1][0] + V[0][1]);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == TOP_RIGHT_CORNER) {
    Vnew = 0.5 * (V[width-2][0] + V[width-1][1]);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == BOTTOM_LEFT_CORNER) {
    Vnew = 0.5 * (V[0][height-2] + V[1][height-1]);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == BOTTOM_RIGHT_CORNER) {
    Vnew = 0.5 * (V[width-2][height-1] + V[width-1][height-2]);
    return g * Vnew + (1 - g) * V[i][j];
  }

  /* The four sides */

  else if (oddity_value == ORDINARY_POINT_LEFT_EDGE) {
    Vnew = 0.25 * (V[0][j-1] + V[0][j+1] + 2 * V[1][j]);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == ORDINARY_POINT_RIGHT_EDGE) {
    Vnew = 0.25 * (V[width-1][j+1] + V[width-1][j-1] + 2 * V[width-2][j]);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == ORDINARY_POINT_TOP_EDGE) {
    Vnew = 0.25 * (V[i-1][0] + V[i+1][0] + 2 * V[i][1]);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == ORDINARY_POINT_BOTTOM_EDGE) {
    Vnew = 0.25 * (V[i-1][height-1] + V[i+1][height-1] + 2 * V[i][height-2]);
    return g * Vnew + (1 - g) * V[i][j];
  }

  else if (oddity_value == ORDINARY_INTERIOR_POINT ||
           (oddity_value >= DIFFERENT_DIELECTRIC_ABOVE_AND_RIGHT &&
            oddity_value < UNDEFINED_ODDITY &&
            dielectrics_to_consider_just_now == 1)) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = (Va + Vb + Vl + Vr) / 4.0;
    return g * Vnew + (1 - g) * V[i][j];
  }

  /* Nodes with metal on one or two sides (see the discussion in the original
     update_voltage_array.c about the 4/3, 2/3 weighting). */

  else if (oddity_value == METAL_ABOVE) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Va/3 + 2*Vb/3 + Vl + Vr);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_BELOW) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vb/3 + 2*Va/3 + Vl + Vr);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_LEFT) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vl/3 + 2*Vr/3 + Va + Vb);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_RIGHT) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vr/3 + 2*Vl/3 + Va + Vb);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_ABOVE_AND_RIGHT) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vr/3 + 4*Va/3 + 2*Vl/3 + 2*Vb/3);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_ABOVE_AND_LEFT) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vl/3 + 4*Va/3 + 2*Vr/3 + 2*Vb/3);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_BELOW_AND_LEFT) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vl/3 + 4*Vb/3 + 2*Vr/3 + 2*Va/3);
    return g * Vnew + (1 - g) * V[i][j];
  }
  else if (oddity_value == METAL_BELOW_AND_RIGHT) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    Vnew = 0.25 * (4*Vb/3 + 4*Vr/3 + 2*Va/3 + 2*Vl/3);
    return g * Vnew + (1 - g) * V[i][j];
  }

  /* A change of permittivity around the node (multiple-dielectric stage). */

  else if (dielectrics_to_consider_just_now > 1) {
    Va = V[i][j-1]; Vb = V[i][j+1]; Vl = V[i-1][j]; Vr = V[i+1][j];
    ERa = Er[i][j-1]; ERb = Er[i][j+1]; ERl = Er[i-1][j]; ERr = Er[i+1][j];
    Vnew = (Va*ERa + Vb*ERb + Vl*ERl + Vr*ERr) / (ERa + ERb + ERl + ERr);
    return g * Vnew + (1 - g) * V[i][j];
  }

  else if ((dielectrics_to_consider_just_now == 1 && oddity_value == UNDEFINED_ODDITY) ||
           (dielectrics_to_consider_just_now > 1)) {
    fprintf(stderr, "Internal error in relax_node.h\n");
    fprintf(stderr, "i=%d j=%d oddity[%d][%d]=%d\n", i, j, i, j, oddity[i][j]);
    exit(INTERNAL_ERROR);
  }

  return V[i][j];   /* unrecognised non-conductor node: leave unchanged */
}

#endif /* RELAX_NODE_H */
