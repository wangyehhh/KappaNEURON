#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _caPump1_reg();
extern void _caPump2_reg();
extern void _capulse_reg();
extern void _glupulse_reg();
extern void _napulse_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," caPump1.mod");
fprintf(stderr," caPump2.mod");
fprintf(stderr," capulse.mod");
fprintf(stderr," glupulse.mod");
fprintf(stderr," napulse.mod");
fprintf(stderr, "\n");
    }
_caPump1_reg();
_caPump2_reg();
_capulse_reg();
_glupulse_reg();
_napulse_reg();
}
