#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include "pgapack.h"

PGAInteger pop [][6] =
{{0, 1, 2, 3, 4, 5}
,{1, 3, 2, 0, 4, 5}
,{0, 2, 1, 3, 4, 5}
,{0, 1, 2, 3, 4, 5}
,{0, 1, 2, 4, 3, 5}
,{0, 1, 3, 2, 4, 5}
,{0, 1, 3, 4, 2, 5}
};
int popsize = sizeof (pop) / (sizeof (PGAInteger) * 6);

int main (int argc, char **argv)
{
    int i;
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, 6, PGA_MINIMIZE);
    PGASetCrossoverType (ctx, PGA_CROSSOVER_EDGE);
    PGASetRandomSeed    (ctx, 2);
    PGASetUp (ctx);
    for (i=0; i<popsize-1; i++) {
        memcpy (ctx->ga.oldpop [0].chrom, pop + 0,     sizeof (pop [0]));
        memcpy (ctx->ga.oldpop [1].chrom, pop + i + 1, sizeof (pop [0]));
        PGAPrintString (ctx, stdout, 0, PGA_OLDPOP);
        PGAPrintString (ctx, stdout, 1, PGA_OLDPOP);
        PGAIntegerEdgeCrossover (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
}
