#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from:

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
void init_flag(char **flag, int imax, int jmax, double delx, double dely,
    int *ibound)
{
    int i, j;

    // Constants for nozzle design
    double nozzle_start                 = 0.1 * imax;
    double nozzle_middle                = 0.2 * imax;
    double nozzle_end                   = 0.4 * imax;
    double nozzle_start_size            = 0.4 * jmax;
    double nozzle_middle_size           = 0.1 * jmax;
    double nozzle_end_size              = 0.5 * jmax;
    double nozzle_converging_length     = nozzle_middle - nozzle_start;
    double nozzle_diverging_length      = nozzle_end - nozzle_middle;

    // Diameter calculation of nozzle
    double diameter, diameter_start_ratio;

    for (i=1;i<=imax;i++) {
      for (j=1;j<=jmax;j++) {
        if (i >= nozzle_start && i <= nozzle_end) {
          if (i < nozzle_middle) {
            diameter_start_ratio = (i - nozzle_start) / nozzle_converging_length;
            diameter = (nozzle_middle_size * diameter_start_ratio) +
                       (nozzle_start_size * (1 - diameter_start_ratio));
          } else {
            diameter_start_ratio = (i - nozzle_middle) / nozzle_diverging_length;
            diameter = (nozzle_end_size * diameter_start_ratio) +
                       (nozzle_middle_size * (1 - diameter_start_ratio));
          }

          // Check if value is outside the diameter
          if (j < (jmax / 2) - diameter / 2 || j > (jmax / 2) + diameter / 2) {
            flag[i][j] = C_B;
          } else {
            flag[i][j] = C_F;
          }

          continue;
        }

        flag[i][j] = C_F;
      }
    }

    /* Mark the north & south boundary cells */
    for (i=0; i<=imax+1; i++) {
        flag[i][0]      = C_B;
        flag[i][jmax+1] = C_B;
    }
    /* Mark the east and west boundary cells */
    for (j=1; j<=jmax; j++) {
        flag[0][j]      = C_B;
        flag[imax+1][j] = C_B;
    }

    /* flags for boundary cells */
    *ibound = 0;
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (!(flag[i][j] & C_F)) {
                (*ibound)++;
                if (flag[i-1][j] & C_F) flag[i][j] |= B_W;
                if (flag[i+1][j] & C_F) flag[i][j] |= B_E;
                if (flag[i][j-1] & C_F) flag[i][j] |= B_S;
                if (flag[i][j+1] & C_F) flag[i][j] |= B_N;
            }
        }
    }
}
