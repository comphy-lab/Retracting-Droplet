/** Title: Drop retraction
 * Initial condition: a pancake shaped drop of heigh h0.
 * using a tanh transition function
# Author: Aman Bhargava & Vatsal Sanjay
# Physics of Fluids
# Last Updated: Jul 27, 2024
version 1.0
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "../src-local/contact-fixed.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

//

#define tmax 100.0                                    // maximum time

#define tsnap1 (1e-2)
#define tsnap2 (1e-3)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-3)                            // error tolerances in velocity

#define h0 1e0
#define RhoR (1e-3)
#define MuR (2e-2)

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

double Oh1, Oh2, Bo, Gamma, Ldomain, DeltaMin, thetaA = 60.0*pi/180.0, thetaR = 40.0*pi/180.0;
int MAXlevel;

double contact_angle_tanh(double u, double thetaA, double thetaB, double width) {
    return thetaA + 0.5 * (thetaB - thetaA) * (1 + tanh(u / width));
}

vector h[];
// h.t[left] = contact_angle((u.y[0,0] > 0.0001) ? thetaA : (u.y[0,0] < -0.0001) ? thetaR : (thetaA + thetaR)/2);
// h.t[left] = contact_angle((u.y[0,0] > 0.0001) ? thetaA : thetaR);

int main(int argc, char const *argv[]) {

  Oh1 = 7.6e-3;
  Oh2 = Oh1*MuR;
  Bo = 5.7e-3;
  Gamma = 6;

  // H = 1.79*pow(3*Gamma*Gamma+1, -1.0/3.0)*1e-3;
  Ldomain = Gamma+3;
  MAXlevel = 9;

  //Gamma = atof(argv[1]);
  //Ldomain = atof(argv[2]);
  //MAXlevel = atoi(argv[3]);

  L0=Ldomain;
  X0=0.0; Y0=0.;
  init_grid (1 << (MAXlevel));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1e0; mu1 = Oh1;
  rho2 = RhoR; mu2 = Oh2;
  G.x = -Bo;

  f.sigma = 1.0;
  f.height = h;

  
  DeltaMin = L0/pow(2,MAXlevel);

  fprintf(ferr, "Level %d tmax %g. Ohf %3.2e, Gamma %g, Bo %3.2e, L0 %g, MaxLevel %d\n", MAXlevel, tmax, Oh1, Gamma, Bo, L0, MAXlevel);
  run();

}



event init(t = 0){
  if(!restore (file = "dump")){
    /**
    We can now initialize the volume fractions in the domain. */
    refine(y<Gamma && x < h0 + 10*DeltaMin && level<MAXlevel);

    fraction(f, y > Gamma-h0 ? h0-sqrt((sq(x)+sq(y-(Gamma-h0)))) : h0-x); //pancake
    // fraction(f, y < Gamma*sqrt((1-sq(x/h0)))); //ellipse

    h.t[left] = contact_angle(contact_angle_tanh(0.0, thetaA, thetaR, 0.1));
  }
}

event properties(i++){
  h.t[left] = contact_angle(contact_angle_tanh(u.y[], thetaA, thetaR, 0.1));
}

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MAXlevel-4);
}

// Outputs
event writingFiles (t = 0; t += tsnap1; t <= tmax + tsnap1) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0; t += tsnap2; t <= tmax + tsnap1) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f[]);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

}
