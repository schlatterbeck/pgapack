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

void edge_test (int argc, char **argv)
{
    int i, j;
    size_t sz;
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, 6, PGA_MINIMIZE);
    PGAInteger edges[][2] = {{1,3}, {3,4}};
    size_t elen = sizeof (edges) / (2 * sizeof (PGAInteger));

    PGASetCrossoverType (ctx, PGA_CROSSOVER_EDGE);
    PGASetRandomSeed    (ctx, 2);
    PGAIntegerSetFixedEdges (ctx, elen, edges, PGA_FALSE);
    printf ("Fixed edges:\n");
    for (sz=0; sz<elen; sz++) {
        PGAFixedEdge *p = ctx->ga.edges + sz;
        printf
            ( "%zu: %ld %ld %zd %zd\n"
            , sz, p->lhs, p->rhs
            , p->next ?  p->next - ctx->ga.edges : -1
            , p->prev ?  p->prev - ctx->ga.edges : -1
            );
    }
    PGASetUp (ctx);
    printf ("Initialized Population:\n");
    for (i=0; i<popsize; i++) {
        for (j=0; j<ctx->ga.StringLen; j++) {
            printf ("%d ", PGAGetIntegerAllele (ctx, i, PGA_OLDPOP, j));
        }
        printf ("\n");
    }
    ctx->ga.n_edges = 0;
    printf ("Edge-Crossover test:\n");
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

void pmx_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    PGAInteger *parent0, *parent1;
    printf ("PMX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_PMX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    /* Init random number generator and print first numbers */
    PGARandom01 (ctx, 1);
    for (i=0; i<14; i++) {
        printf ("rand: %d\n", PGARandomInterval (ctx, 0, l - 1));
    }
    /* Re-init, we'll get the same numbers as above */
    PGARandom01 (ctx, 1);
    /* Make crossovers */
    for (i=0; i<7; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerPartiallyMappedCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
}

void mx_test (int argc, char **argv)
{
    int l = 6;
    int i, j;
    PGAInteger *parent0, *parent1;
    printf ("MX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_MODIFIED);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    /* Init random number generator and print first 10 numbers */
    PGARandom01 (ctx, 3);
    for (i=0; i<5; i++) {
        printf ("rand: %d\n", PGARandomInterval (ctx, 0, l - 1));
    }
    /* Re-init, we'll get the same numbers as above */
    PGARandom01 (ctx, 3);
    /* do not do duplicate test: */
    PGARandomInterval (ctx, 0, l - 1);
    /* Make crossovers */
    for (i=0; i<4; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerModifiedCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
}

void ox_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    PGAInteger *parent0, *parent1;
    printf ("OX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_ORDER);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    PGARandom01 (ctx, 1);
    /* Make crossovers */
    for (i=0; i<7; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerOrderCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
}

void nox_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    PGAInteger *parent0, *parent1;
    printf ("NOX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_NOX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    PGARandom01 (ctx, 1);
    /* Make crossovers */
    for (i=0; i<7; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerNonWrappingOrderCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
}

int main (int argc, char **argv)
{
    edge_test (argc, argv);
    pmx_test  (argc, argv);
    mx_test   (argc, argv);
    ox_test   (argc, argv);
    nox_test  (argc, argv);
}
