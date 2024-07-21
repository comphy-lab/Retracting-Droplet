#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"


vector h[];
double theta0;
char filename[80];

int main(int argc, char const *argv[]) {
  sprintf (filename, "%s", argv[1]);
  restore (file = filename);
  theta0 = atof(argv[2]);

  h.t[left] = contact_angle (theta0*pi/180.);
  heights (f, h);
  f.height = h;

  FILE * fp = ferr;
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2){
        double x0 = x + (segment[0].x+segment[1].x)*Delta/2.;
        double y0 = y + (segment[0].y+segment[1].y)*Delta/2.;
        fprintf(fp, "%g %g %g %g\n", x0, y0, n.x, n.y);
      }
    }
  }
  fflush (fp);
  fclose (fp);

}