#include "utils.h"
#include "output.h"
#include "fractions.h"
#include "tag.h"

char filename[80];
scalar f[];

// Structure to hold the result of the tagging function
typedef struct {
  int main_phase;
  scalar d;
} TagResult;

// Function to tag connected regions and find the largest one
TagResult tag_largest_region(scalar field, double threshold) {
  scalar d[];
  foreach() {
    d[] = (field[] > threshold);
  }

  int n = tag(d), size[n];
  for (int i = 0; i < n; i++) {
    size[i] = 0;
  }

  foreach_leaf() {
    if (d[] > 0) {
      size[((int) d[]) - 1]++;
    }
  }

  int max_size = 0;
  int main_phase = 0;
  for (int i = 0; i < n; i++) {
    if (size[i] > max_size) {
      max_size = size[i];
      main_phase = i + 1;
    }
  }

  TagResult result = {main_phase, d};
  return result;
}

int main(int argc, char const *argv[]) {
  sprintf(filename, "%s", argv[1]);
  restore(file = filename);

  double threshold = 1e-4;

  // Tag all liquid parts and find the largest connected region
  TagResult liquid_result = tag_largest_region(f, threshold);
  int main_phase_liquid = liquid_result.main_phase;
  scalar d_liquid = liquid_result.d;

  // Tag all gas parts and find the largest connected region
  scalar dg[];
  foreach() {
    dg[] = (f[] < 1. - threshold);
  }
  TagResult gas_result = tag_largest_region(dg, threshold);
  int main_phase_gas = gas_result.main_phase;
  scalar d_gas = gas_result.d;

  FILE *fp = ferr;
  foreach() {
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d_liquid[] == main_phase_liquid && d_gas[] == main_phase_gas) {
      coord n = interface_normal(point, f);
      coord t;
      double alpha = plane_alpha(f[], n);
      coord segment[2];
      if (facets(n, alpha, segment) == 2) {
        double x0 = x + (segment[0].x + segment[1].x) * Delta / 2.;
        double y0 = y + (segment[0].y + segment[1].y) * Delta / 2.;

        t.x = -n.y / sqrt(sq(n.x) + sq(n.y));
        t.y = n.x / sqrt(sq(n.x) + sq(n.y));
        fprintf(fp, "%g %g %g %g\n", x0, y0, t.x, t.y);
      }
    }
  }
  fflush(fp);
  fclose(fp);

  // dump(file = "interfaceAndAngle");

  // return 0;
}
