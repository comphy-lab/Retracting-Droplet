/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

double thetaA = 60.0*pi/180.0, thetaR = 40.0*pi/180.0;

// scalar f[];
vector h[];
char filename[80];

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}
