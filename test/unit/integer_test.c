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
    PGADestroy (ctx);
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
    PGADestroy (ctx);
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
    PGADestroy (ctx);
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
    PGADestroy (ctx);
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
    PGADestroy (ctx);
}

void cx_test (int argc, char **argv)
{
    int l = 12;
    int i, j;
    PGAInteger *parent0, *parent1;
    static const PGAInteger p2 [] = {7, 10, 2, 4, 5, 3, 1, 11, 0, 8, 6, 9};
    printf ("CX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_CYCLE);
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
    /* Make crossovers */
    for (i=0; i<5; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [j] = p2 [j];
        }
        PGAIntegerCycleCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
    PGADestroy (ctx);
}

void pbx_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    static const int p [] = {8, 6, 4, 2, 7, 5, 3, 1};
    PGAInteger *parent0, *parent1;
    printf ("PBX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_PBX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    /* Init random number generator and print first 10 numbers */
    PGARandom01 (ctx, 3);
    for (i=0; i<3; i++) {
        printf ("rand: ");
        for (j=0; j<l; j++) {
            printf ("%d", PGARandomFlip (ctx, 0.5));
        }
        printf ("\n");
    }
    /* Re-init, we'll get the same numbers as above */
    PGARandom01 (ctx, 3);
    /* Make crossovers */
    for (i=0; i<3; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerPositionBasedCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
    PGADestroy (ctx);
    /* Verify that example in Dav91 p.81 (Uniform order based crossover)
     * is the same for the first child when using PBX
     */
    l = 8;
    ctx = PGACreate (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_PBX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    PGARandom01 (ctx, 5);
    /* The 841th gives exactly the bit pattern 01101100 from the example */
    for (i=0; i<840; i++) {
        PGARandomFlip (ctx, 0.5);
    }
    printf ("\n");
    for (j=0; j<l; j++) {
        parent0 [j] = j;
        parent1 [j] = p [j] - 1;
    }
    PGAIntegerPositionBasedCrossover
        (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
    PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
    PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    PGADestroy (ctx);
}

void uox_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    static const int p [] = {8, 6, 4, 2, 7, 5, 3, 1};
    PGAInteger *parent0, *parent1;
    printf ("UOX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_UOX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    /* Init random number generator and print first 10 numbers */
    PGARandom01 (ctx, 3);
    for (i=0; i<3; i++) {
        printf ("rand: ");
        for (j=0; j<l; j++) {
            printf ("%d", PGARandomFlip (ctx, 0.5));
        }
        printf ("\n");
    }
    /* Re-init, we'll get the same numbers as above */
    PGARandom01 (ctx, 3);
    /* Make crossovers */
    for (i=0; i<3; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerUniformOrderBasedCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
    PGADestroy (ctx);
    /* Verify that example in Dav91 p.81 (Uniform order based crossover)
     * is the same for both children
     */
    l = 8;
    ctx = PGACreate (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_PBX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    PGARandom01 (ctx, 5);
    /* The 841th gives exactly the bit pattern 01101100 from the example */
    for (i=0; i<840; i++) {
        PGARandomFlip (ctx, 0.5);
    }
    printf ("\n");
    for (j=0; j<l; j++) {
        parent0 [j] = j;
        parent1 [j] = p [j] - 1;
    }
    PGAIntegerUniformOrderBasedCrossover
        (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
    PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
    PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    PGADestroy (ctx);
}

void obx_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    static const char p [] = "eibdfajgch";
    PGAInteger *parent0, *parent1;
    printf ("OBX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_OBX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    /* Init random number generator and print first 10 numbers */
    PGARandom01 (ctx, 3);
    /* Make crossovers */
    for (i=0; i<3; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [l - 1 - j] = j;
        }
        PGAIntegerOrderBasedCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
    PGARandom01 (ctx, 1);
    /* The 31th gives exactly the bit pattern 0110100100 from the example */
    for (i=0; i<30; i++) {
        PGARandomFlip (ctx, 0.5);
    }
    printf ("\n");
    for (j=0; j<l; j++) {
        parent0 [j] = j;
        parent1 [j] = p [j] - 'a';
    }
    PGAIntegerOrderBasedCrossover
        (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
    PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
    PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    PGADestroy (ctx);
}

void aex_test (int argc, char **argv)
{
    int l = 6;
    int i, j;
    PGAInteger *parent0, *parent1;
    static const int p [] = {1, 2, 5, 4, 6, 3};
    printf ("AEX test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetCrossoverType (ctx, PGA_CROSSOVER_AEX);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent0 = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    parent1 = (PGAInteger *)PGAGetIndividual (ctx, 1, PGA_OLDPOP)->chrom;
    /* Init random number generator and print first 10 numbers */
    PGARandom01 (ctx, 2);
    for (i=0; i<4; i++) {
        printf ("rand: %d\n", PGARandomInterval (ctx, 0, l - 1));
    }
    PGARandom01 (ctx, 2);
    /* Make crossovers */
    for (i=0; i<3; i++) {
        /* Reset parents */
        for (j=0; j<l; j++) {
            parent0 [j] = j;
            parent1 [j] = p [j] - 1;
        }
        PGAIntegerAlternatingEdgeCrossover
            (ctx, 0, 1, PGA_OLDPOP, 0, 1, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 0, PGA_NEWPOP);
        PGAPrintString (ctx, stdout, 1, PGA_NEWPOP);
    }
}

void mutation_scramble_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    char *scramble [] = { "5 @ 4", "2 @ 2", "3 @ 2", "2 @ 5", "2 @ 3" };
    PGAInteger *parent;
    printf ("Scramble mutation test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetMutationType  (ctx, PGA_MUTATION_SCRAMBLE);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    PGARandom01 (ctx, 1);
    /* Make mutations */
    for (i=0; i<5; i++) {
        /* Reset parent */
        for (j=0; j<l; j++) {
            parent [j] = j;
        }
        printf ("scramble: %s\n", scramble [i]);
        PGAIntegerMutation (ctx, 0, PGA_OLDPOP, 0.5);
        PGAPrintString (ctx, stdout, 0, PGA_OLDPOP);
    }
}

void mutation_position_test (int argc, char **argv)
{
    int l = 10;
    int i, j;
    PGAInteger *parent;
    printf ("Position mutation test\n");
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, l, PGA_MINIMIZE);

    PGASetMutationType  (ctx, PGA_MUTATION_POSITION);
    PGASetRandomSeed    (ctx, 2);
    PGASetPopSize       (ctx, 8);
    PGASetUp (ctx);
    /* Now init two genes and cross them over */
    parent = (PGAInteger *)PGAGetIndividual (ctx, 0, PGA_OLDPOP)->chrom;
    PGARandom01 (ctx, 1);
    for (i=0; i<16; i++) {
        printf ("rand: %d\n", PGARandomInterval (ctx, 0, l - 1));
    }
    PGARandom01 (ctx, 1);
    /* Make mutations */
    for (i=0; i<8; i++) {
        /* Reset parent */
        for (j=0; j<l; j++) {
            parent [j] = j;
        }
        PGAIntegerMutation (ctx, 0, PGA_OLDPOP, 0.5);
        PGAPrintString (ctx, stdout, 0, PGA_OLDPOP);
    }
}

int main (int argc, char **argv)
{
    edge_test (argc, argv);
    pmx_test  (argc, argv);
    mx_test   (argc, argv);
    ox_test   (argc, argv);
    nox_test  (argc, argv);
    cx_test   (argc, argv);
    pbx_test  (argc, argv);
    uox_test  (argc, argv);
    obx_test  (argc, argv);
    aex_test  (argc, argv);
    mutation_scramble_test (argc, argv);
    mutation_position_test (argc, argv);
}
