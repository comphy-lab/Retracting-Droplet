/*
# Author: Aman Bhargava & Vatsal Sanjay
# Physics of Fluids
# Last Updated: Jul 27, 2024
version 1.0
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "src-local/contact-fixed.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

double contact_angle_tanh(double u, double thetaA, double thetaB, double width) {
    return thetaA + 0.5 * (thetaB - thetaA) * (1 + tanh(u / width));
}

double thetaA = 60.0*pi/180.0, thetaR = 40.0*pi/180.0;

vector h[];
double theta0;
char filename[80];

int main(int argc, char const *argv[]) {
  sprintf (filename, "%s", argv[1]);
  restore (file = filename);
  #if TREE
    f.prolongation = fraction_refine;
  #endif
  f.height = h;
  h.t[left] = contact_angle(contact_angle_tanh(u.y[], thetaA, thetaR, 0.1));
  boundary((scalar *){f,h});

  FILE * fp = ferr;
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2){
        double x0 = x + (segment[0].x+segment[1].x)*Delta/2.;
        double y0 = y + (segment[0].y+segment[1].y)*Delta/2.;
        // knowing n.x and n.y normal vector, find the tangent vector
        t.x = -n.y/sqrt(sq(n.x) + sq(n.y));
        t.y = n.x/sqrt(sq(n.x) + sq(n.y));
        fprintf(fp, "%g %g %g %g\n", x0, y0, t.x, t.y);
      }
    }
  }
  fflush (fp);
  fclose (fp);

}