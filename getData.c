/* Title: getting Data from simulation snapshot
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

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay, mu1, mu2, rho1, rho2;
double thetaA = 60.0*pi/180.0, thetaR = 40.0*pi/180.0;
scalar D2c[], vel[], mom[];
vector h[];
scalar * list = NULL;

double contact_angle_tanh(double u, double thetaA, double thetaB, double width) {
    return thetaA + 0.5 * (thetaB - thetaA) * (1 + tanh(u / width));
}

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);
  
  // fprintf(ferr, "xmin %g, ymin %g, xmax %g, ymax %g, ny %d\n", xmin, ymin, xmax, ymax, ny);

  list = list_add (list, D2c);
  list = list_add (list, vel);
  list = list_add (list, mom);

  mu1 = 7.6e-3;
  mu2 = 1.52e-4;

  rho1 = 1e0;
  rho2 = 1e-3;

  f.height = h;

  // boundary conditions
  u.t[top] = neumann(0.);
  uf.t[top] = neumann(0.);
  u.n[top] = neumann(0.);
  uf.n[top] = neumann(0.);
  p[top] = dirichlet(0.);

  u.t[right] = neumann(0.);
  uf.t[right] = neumann(0.);
  u.n[right] = neumann(0.);
  uf.n[right] = neumann(0.);
  p[right] = dirichlet(0.);

  u.t[left] = dirichlet(0.);
  uf.t[left] = dirichlet(0.);
  u.n[left] = dirichlet(0.);
  uf.n[left] = dirichlet(0.);
  p[left] = neumann(0.);
  h.t[left] = contact_angle(contact_angle_tanh(u.y[], thetaA, thetaR, 0.1));

  /*
  Actual run and codes!
  */
  restore (file = filename);

  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/y);
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13);
    D2c[] = 2*(clamp(f[], 0., 1.) * mu1 + clamp((1.-f[]), 0., 1.) * mu2)*D2;

    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10);
    } else {
      D2c[] = -10;
    }
    vel[] = sqrt(sq(u.x[])+sq(u.y[]));
    mom[] =(clamp(f[], 0., 1.) * rho1 + clamp((1.-f[]), 0., 1.) * rho2)*vel[]; 
  }

  FILE * fp = ferr;
  Deltay = (double)((ymax-ymin)/(ny));
  // fprintf(ferr, "Deltay %g\n", Deltay);
  nx = (int)(xmax - xmin)/Deltay;
  // fprintf(ferr, "nx %d\n", nx);
  Deltax = (double)((xmax-xmin)/(nx));
  // fprintf(ferr, "%g\n", Deltax);
  len = list_len(list);
  // fprintf(ferr, "%d\n", len);
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);
}
