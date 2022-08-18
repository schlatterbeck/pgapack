/*  Sequence problems */

#include <pgapack.h>
#include <unistd.h>
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
    int opt;
    PGAInteger edge [][2] = {{0, 1}, {1, 2}, {2, 3}, {4, 5}};
    size_t el = sizeof (edge) / (2 * sizeof (PGAInteger));
    size_t nfix = el;

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
    while ((opt = getopt (argc, argv, "r:e:")) > 0) {
        switch (opt) {
        case 'r':
            seed = atoi (optarg);
            break;
        case 'e':
            nfix = atoi (optarg);
            if (nfix > el) {
                nfix = el;
            }
            break;
        default:
            fprintf (stderr, "Usage: %s [-r seed] [-e nfixed]\n", argv [0]);
            exit (EXIT_FAILURE);
        }
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
    PGAIntegerSetFixedEdges  (ctx, nfix, edge, PGA_TRUE);

    PGASetUp   (ctx);
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
