/*  Sequence problems */

#include <pgapack.h>
#include "optimize.h"

double distance [36][36];

double evaluate (PGAContext *ctx, int p, int pop, double *aux)
{
    int i;
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    PGAInteger *node = (PGAInteger *)ind->chrom;
    int len = PGAGetStringLength (ctx);
    double s = 0;
    for (i=0; i<len; i++) {
        int j = (i+1) % len;
        s += distance [node [i]][node [j]];
    }
    return s;
}

int stop_cond (PGAContext *ctx)
{
    double best = PGAGetBestReport (ctx, PGA_OLDPOP, 0);
    if (best == 36) {
        return PGA_TRUE;
    }
    return PGACheckStoppingConditions (ctx);
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int popsize = 100;
    int maxiter = 2500;
    int x1, y1, x2, y2;
    int seed = 1;

    for (x1=0; x1<6; x1++) {
        for (y1=0; y1<6; y1++) {
            for (x2=0; x2<6; x2++) {
                for (y2=0; y2<6; y2++) {
                    distance [x1 + y1 * 6][x2 + y2 * 6] = sqrt
                        (pow (x1 - x2, 2) + pow (y1 - y2, 2));
                }
            }
        }
    }
#if 0
    for (x1=0; x1<36; x1++) {
        for (y1=0; y1<36; y1++) {
            printf ("(%2d, %2d): %e\n", x1, y1, distance [x1][y1]);
        }
    }
#endif
    if (argc > 1) {
        seed = atoi (argv [1]);
    }

    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, 36, PGA_MINIMIZE);
    
    PGASetRandomSeed         (ctx, seed);
    PGASetPopSize            (ctx, popsize);
    PGASetNumReplaceValue    (ctx, 10);
    PGASetSelectType         (ctx, PGA_SELECT_TOURNAMENT);
    PGASetPopReplaceType     (ctx, PGA_POPREPL_BEST);
    PGASetMaxGAIterValue     (ctx, maxiter);
    PGASetNoDuplicatesFlag   (ctx, PGA_TRUE);
    PGASetCrossoverType      (ctx, PGA_CROSSOVER_EDGE);
    PGASetUserFunction       (ctx, PGA_USERFUNCTION_STOPCOND, stop_cond);
    PGASetTournamentSize     (ctx, 1.7);
    PGASetMixingType         (ctx, PGA_MIX_MUTATE_OR_CROSS);
    //PGASetTournamentWithReplacement (ctx, PGA_FALSE);
    
    PGASetUp   (ctx);
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
