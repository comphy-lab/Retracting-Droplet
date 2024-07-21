/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

scalar KAPPA[];
char filename[80];
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  scalar pos[];
  position (f, pos, {0,1});
  double ymax = statsf(pos).max;
  //fprintf(ferr, "%g %g\n", t, xmax);

  face vector s = {{-1}};
  // to get the entire interface in a parametric way: (x, y, KAPPA[])
  // foreach(){
  //   // if (interfacial(point, f))
  //   // if (KAPPA[] != nodata)
  //   if (f[] > 1e-6 && f[] < 1. - 1e-6) 
  //   {
  //     coord n = facet_normal (point, f, s);
  //     double alpha = plane_alpha (f[], n);
  //     coord segment[2];
  //     if (facets (n, alpha, segment) == 2){
  //       double xInt = x + (segment[0].x*Delta+segment[1].x*Delta)/2.;
  //       double yInt = y + (segment[0].y*Delta+segment[1].y*Delta)/2.;
  //       fprintf (ferr, "%g %g %g\n", xInt, yInt, KAPPA[]);
  //     }
  //   }
  // }
  // dump("newFile");

  double xAtMAX = 0., yMAX = -HUGE;
  foreach(){
    // if (interfacial(point, f))
    // if (KAPPA[] != nodata)
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && x < 1e-2) 
    {
      coord n = facet_normal (point, f, s);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2){
        double xInt, yInt;
        if (segment[0].y > segment[1].y){
          xInt = x + segment[0].x*Delta;
          yInt = y + segment[0].y*Delta;
        } else {
          xInt = x + segment[1].x*Delta;
          yInt = y + segment[1].y*Delta;
        }
        // fprintf (ferr, "%g %g %g\n", xInt, yInt, KAPPA[]);
        if (yInt > yMAX){
          yMAX = yInt;
          xAtMAX = xInt;
        }
      }
    }
  }
  fprintf(ferr, "%g %g\n", t, yMAX);

}
