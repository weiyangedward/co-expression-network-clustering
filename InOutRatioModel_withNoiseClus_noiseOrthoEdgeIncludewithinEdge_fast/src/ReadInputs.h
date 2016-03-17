#ifndef __READINPUTS__
#define __READINPUTS__

#include "cs-grn.h"

using namespace std;
//using namespace __gnu_cxx;

SpeciesNetwork *ReadNetworkFile(char *fn);

Orthology *ReadOrthologyFile(char *orthfn, SpeciesNetwork **sns, int numSpc);

#endif
