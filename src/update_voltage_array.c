/*
The per-node voltage update (given the voltages and permittivities at the four
adjacent points) used to live inline in this file. It has been moved verbatim
to the shared header relax_node.h so that both this single-threaded
4-direction sweep and the OpenMP red-black sweep (finite_difference_openmp.c)
compute each node identically. See relax_node.h for the per-node maths and the
comment about multiple dielectrics.
*/


#include "config.h"

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#include "definitions.h"


extern int width, height;
extern double **Er;
extern unsigned char  **oddity;
extern int dielectrics_to_consider_just_now;
extern double r;
extern int coupler;

#include "exit_codes.h"
#include "relax_node.h"

/* The following function updates the voltage on the matrix V_to given data about the
oddity of the location i,j and the voltages in the matrix V_from. It does this for n interations
between rows jmin and jmax inclusive and between columns imain and imax inclusive */

void update_voltage_array(int nmax, int imin, int imax, int jmin, int jmax, double **V_from, double **V_to)
{
  int k, i, j, n;
  double g;

  if (dielectrics_to_consider_just_now==1)
    g=r;
  else
    g=1;
  for(n=0; n  < nmax; ++n)
    for(k=0; k < 4; ++k)
      for (i = k&1 ? imax : imin;   k&1 ? i >=imin : i <= imax ;  k&1 ? i-- : i++)
        for (j = (k==0 || k ==3) ? jmin : jmax; (k ==0 || k == 3)  ? j <= jmax : j >= jmin ; (k == 0 || k ==3) ?  j++ : j--){
          V_to[i][j] = new_node_value(i, j, g, V_from);
        } /* end of j loop */
}
