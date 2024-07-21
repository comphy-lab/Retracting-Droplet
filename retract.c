#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "adapt_wavelet_limited_v2.h"
#include "reduced.h"
// #include "contact.h"
#include "contact-fixed.h"
#include "axi.h"
//

#define tmax 100.0                                                 // maximum time

#define tsnap1 (1e-2)
#define tsnap2 (1e-2)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-3)                            // error tolerances in velocity

#define h0 1e0
#define RhoR (1e-3)
#define MuR (2e-2)

// boundary conditions


// u.t[bottom] = neumann(0.);
// uf.t[bottom] = neumann(0.);
// u.n[bottom] = dirichlet(0.);
// uf.n[bottom] = dirichlet(0.);
// p[bottom] = neumann(0.);

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

double Oh1, Oh2, Bo, Gamma, Ldomain, DeltaMin, theta0 = 25.0*pi/180.0, thetaA = 60.0*pi/180.0, thetaR = 40.0*pi/180.0, ymax;
int MAXlevel;

// double gettingFacets(char* filename) {
//   char command[50];
//   sprintf(command, "./getFacet %s", filename);

//   FILE *fp;
//   char temp2[MAX_NUM_LINES][MAX_LINE_LENGTH];
//   int line_count = 0;

//   fp = popen(command, "r");
//   if (fp == NULL) {
//     printf("Failed to run command\n" );
//     exit(1);
//   }

//   while (fgets(temp2[line_count], MAX_LINE_LENGTH, fp) != NULL) {
//     line_count++;
//   }

//   pclose(fp);

//   double r0 = 0.0;
//   if (line_count > 100) {
//     int skip = 0;
//     for (int n1 = 0; n1 < line_count; n1++) {
//       char* temp3 = strtok(temp2[n1], " ");
//       if (strcmp(temp3, "") == 0) {
//         skip = 0;
//       } else {
//         if (!skip) {
//           char* temp4 = strtok(temp2[n1+1], " ");
//           double r1 = atof(strtok(temp3, " "));
//           double z1 = atof(strtok(NULL, " "));
//           double r2 = atof(strtok(temp4, " "));
//           double z2 = atof(strtok(NULL, " "));
//           double midpoint[2] = {-(r1 + r2) / 2, (z1 + z2) / 2};

//           if (midpoint[1] < 0.1) {
//             r0 = midpoint[0];
//             break;
//           }
//         }
//       }
//     }
//   }

//   return r0;
// }


vector h[];
// h.t[left] = contact_angle((u.y[0,0] > 0.0001) ? thetaA : (u.y[0,0] < -0.0001) ? thetaR : (thetaA + thetaR)/2);
h.t[left] = contact_angle((u.y[0,0] > 0.0001) ? thetaA : thetaR);

int main(int argc, char const *argv[]) {

  // ymax = gettingFacets("dump_test_25");
  //if (argc != 4){
    //fprintf(ferr, "Usage: %s <Gamma> <Ldomain> <MAXlevel>\n", argv[0]);
    //return 1;
  //}
  f.height = h;

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
  sprintf (comm, "mkdir -p test_angle");
  system(comm);

  rho1 = 1e0; mu1 = Oh1;
  rho2 = RhoR; mu2 = Oh2;
  G.x = -Bo;

  f.sigma = 1.0;
  
  DeltaMin = L0/pow(2,MAXlevel);

  fprintf(ferr, "Level %d tmax %g. Ohf %3.2e, Gamma %g, Bo %3.2e, L0 %g, MaxLevel %d\n", MAXlevel, tmax, Oh1, Gamma, Bo, L0, MAXlevel);
  run();

}



event init(t = 0){
  if(!restore (file = "dump_test_angle")){
    /**
    We can now initialize the volume fractions in the domain. */
    refine(y<Gamma && x < h0 + 10*DeltaMin && level<MAXlevel);

    fraction(f, y > Gamma-h0 ? h0-sqrt((sq(x)+sq(y-(Gamma-h0)))) : h0-x); //pancake
    // fraction(f, y < Gamma*sqrt((1-sq(x/h0)))); //ellipse

  }
}

int refRegion(double x, double y, double z){
  return (x < 1e-2? MAXlevel+2:
  x < 1e-1? MAXlevel+1:
  x < 1.5e0? MAXlevel:
  x < 3e0? MAXlevel-1:
  x < 6e0? MAXlevel-2:
  // x < 12e0? MAXlevel-3:
  MAXlevel-3
  );
}

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  // adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, omega},
  //   (double[]){fErr, VelErr, VelErr, KErr, OmegaErr},
  //   refRegion, MINlevel);
  adapt_wavelet_limited((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    refRegion, MAXlevel-4);
}

// Outputs
event writingFiles (t = 0; t += tsnap1; t <= tmax + tsnap1) {
  dump (file = "dump_test_angle");
  char nameOut[80];
  sprintf (nameOut, "test_angle/snapshot-%5.4f", t);
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
