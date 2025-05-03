/*  Magic Squares  */

#include <pgapack.h>
#include <unistd.h>

static PGAInteger square [5][5];
static double rowsum [5];
static double colsum [5];
static double diasum [2];

void phenotype (PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    PGAInteger *chrom = (PGAInteger *)ind->chrom;
    int len = PGAGetStringLength (ctx);
    int x, y;

    assert (len == 25);
    diasum [0] = diasum [1] = 0;
    for (y=0; y<5; y++) {
        rowsum [y] = 0;
        colsum [y] = 0;
        for (x=0; x<5; x++) {
            square [y][x] = chrom [5*y + x];
        }
    }
    for (y=0; y<5; y++) {
        for (x=0; x<5; x++) {
            rowsum [y] += square [y][x];
            colsum [x] += square [y][x];
            if (x == y) {
                diasum [0] += square [y][x];
            }
            if (y == 5 - x - 1) {
                diasum [1] += square [y][x];
            }
        }
    }
}

double evaluate (PGAContext *ctx, int p, int pop, double *aux)
{
    int i;
    int best = 5 * (5*5 + 1) / 2;
    double s = 0.0;
    phenotype (ctx, p, pop);
    for (i=0; i<2; i++) {
        s += abs (diasum [i] - best);
    }
    for (i=0; i<5; i++) {
        s += abs (colsum [i] - best);
        s += abs (rowsum [i] - best);
    }
    return s;
}

int stop_cond (PGAContext *ctx)
{
    double best = PGAGetBestReport (ctx, PGA_OLDPOP, 0);
    if (best == 0) {
        return PGA_TRUE;
    }
    return PGACheckStoppingConditions (ctx);
}

void print_str (PGAContext *ctx, FILE *fp, int p, int pop)
{
    int x, y;
    phenotype (ctx, p, pop);
    for (y=0; y<5; y++) {
        fprintf (fp, "    ");
        for (x=0; x<5; x++) {
            fprintf (fp, "%3ld ", square [y][x]);
        }
        fprintf (fp, "%3.0f\n", rowsum [y]);
    }
    fprintf (fp, "%3.0f ", diasum [1]);
    for (x=0; x<5; x++) {
        fprintf (fp, "%3.0f ", colsum [x]);
    }
    fprintf (fp, "%3.0f\n", diasum [0]);
    fflush (fp);
}

static const char *mut_types [] =
    {"PERMUTE", "SCRAMBLE", "POSITION"};
static int mut_values [] =
    {PGA_MUTATION_PERMUTE, PGA_MUTATION_SCRAMBLE, PGA_MUTATION_POSITION};
#define MUTATION_SIZE ((int)(sizeof (mut_values) / sizeof (int)))
static const char *cross_types [] =
    { "PMX"
    , "MODIFIED"
    , "ORDER"
    , "CYCLE"
    , "OBX"
    , "PBX"
    , "UOX"
    , "AEX"
    , "NOX"
    };
static int cross_values [] =
    { PGA_CROSSOVER_PMX
    , PGA_CROSSOVER_MODIFIED
    , PGA_CROSSOVER_ORDER
    , PGA_CROSSOVER_CYCLE
    , PGA_CROSSOVER_OBX
    , PGA_CROSSOVER_PBX
    , PGA_CROSSOVER_UOX
    , PGA_CROSSOVER_AEX
    , PGA_CROSSOVER_NOX
    };
#define CROSS_SIZE ((int)(sizeof (cross_values) / sizeof (int)))

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int popsize = 100;
    int seed = 1;
    int opt;
    int crossover = 0, mutation = 0;
    int i;

    while ((opt = getopt (argc, argv, "c:m:r:")) > 0) {
        switch (opt) {
        case 'c':
            crossover = atoi (optarg);
            if (crossover < 0 || crossover > CROSS_SIZE - 1) {
                fprintf (stderr, "Invalid crossover type: %d\n", crossover);
                fprintf (stderr, "Valid:\n");
                for (i=0; i<CROSS_SIZE; i++) {
                    fprintf (stderr, "%d: %s\n", i, cross_types [i]);
                }
                exit (EXIT_FAILURE);
            }
            break;
        case 'm':
            mutation = atoi (optarg);
            if (mutation < 0 || mutation > MUTATION_SIZE - 1) {
                fprintf (stderr, "Invalid mutation type: %d\n", mutation);
                fprintf (stderr, "Valid:\n");
                for (i=0; i<MUTATION_SIZE; i++) {
                    fprintf (stderr, "%d: %s\n", i, mut_types [i]);
                }
                exit (EXIT_FAILURE);
            }
            break;
        case 'r':
            seed = atoi (optarg);
            break;
        default:
            fprintf
                (stderr, "Usage: %s [-c cross] [-m mut] [-r seed]\n", argv [0]);
            exit (EXIT_FAILURE);
        }
    }
    crossover = cross_values [crossover];
    mutation  = mut_values   [mutation];

    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_INTEGER, 25, PGA_MINIMIZE);

    PGASetRandomSeed           (ctx, seed);
    PGASetPopSize              (ctx, popsize);
    PGASetNumReplaceValue      (ctx, (int)(popsize * 0.9));
    PGASetMaxGAIterValue       (ctx, 1000);
    PGASetMaxNoChangeValue     (ctx, 400);
    PGASetPrintOptions         (ctx, PGA_REPORT_STRING);
    /* PGASetPrintFrequencyValue  (ctx, 1); */
    PGASetSelectType           (ctx, PGA_SELECT_TRUNCATION);
    PGASetPopReplaceType       (ctx, PGA_POPREPL_RTR);
    /* PGASetRTRWindowSize       (ctx, 30); */
    PGASetIntegerInitPermute   (ctx, 1, 25);
    PGASetMutationType         (ctx, mutation);
    PGASetMutationScrambleMax  (ctx, 5);
    PGASetCrossoverType        (ctx, crossover);
    PGASetNoDuplicatesFlag     (ctx, PGA_TRUE);
    PGASetTruncationProportion (ctx, 0.5);
    PGASetUserFunction         (ctx, PGA_USERFUNCTION_STOPCOND, stop_cond);
    PGASetUserFunction         (ctx, PGA_USERFUNCTION_PRINTSTRING, print_str);
    /* PGASetTournamentSize       (ctx, 2); */
    /* PGASetMixingType           (ctx, PGA_MIX_MUTATE_OR_CROSS); */

    PGASetUp   (ctx);
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
