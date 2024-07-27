/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "src-local/contact-fixed.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

double thetaA = 60.0*pi/180.0, thetaR = 40.0*pi/180.0;

// scalar f[];
vector h[];
char filename[80];

double contact_angle_tanh(double u, double thetaA, double thetaB, double width) {
    return thetaA + 0.5 * (thetaB - thetaA) * (1 + tanh(u / width));
}

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  #if TREE
    f.prolongation = fraction_refine;
  #endif
  f.height = h;
  h.t[left] = contact_angle(contact_angle_tanh(u.y[], thetaA, thetaR, 0.1));
  boundary((scalar *){f});
  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}
