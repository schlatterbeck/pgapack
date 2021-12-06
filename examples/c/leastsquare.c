/* Least Squares
 * Test various selection methods and, if applicable, fitness scaling
 * algorithms of PGApack against a simple minimization problem of
 * minimizing a real gene with a least squares evaluation function
 * Note that we test maximization also by returning the negative of the
 * evaluation function.
 */
#include <pgapack.h>

double lower [5]      = { -5, -5, -5, -5, -5 };
double upper [5]      = {  5,  5,  5,  5,  5 };
int dimension = sizeof (lower) / sizeof(double);

typedef struct intconf
{ int   conf;
  char *name;
} intconf;

static const intconf selection_type [] =
{ { PGA_SELECT_PROPORTIONAL, "Proportional"}
, { PGA_SELECT_SUS,          "Stochastic Universal"}
, { PGA_SELECT_TOURNAMENT,   "Tournament"}
, { PGA_SELECT_PTOURNAMENT,  "Probabilistic Tournament"}
, { PGA_SELECT_TRUNCATION,   "Truncation"}
};
static const int nselect = sizeof (selection_type) / sizeof (intconf);

static const intconf fitness_type [] =
{ { PGA_FITNESS_RAW,        "Raw Fitness"}
, { PGA_FITNESS_NORMAL,     "Normal Fitness"}
, { PGA_FITNESS_RANKING,    "Ranking Fitness"}
};
static const int nfitness = sizeof (fitness_type) / sizeof (intconf);

static const intconf fitness_min [] =
{ { PGA_FITNESSMIN_RECIPROCAL, "Reciprocal Fitness Minimization"}
, { PGA_FITNESSMIN_CMAX,       "CMAX Fitness Minimization"}
};
static const int nfitmin = sizeof (fitness_min) / sizeof (intconf);

double evaluate_min (PGAContext *ctx, int p, int pop)
{
    int i;
    int dimension = PGAGetStringLength (ctx);
    double e = 0;

    for (i=0; i<dimension; i++) {
        double x = PGAGetRealAllele (ctx, p, pop, i);
        e += x * x;
    }
    return e;
}

double evaluate_max (PGAContext *ctx, int p, int pop)
{
    return -evaluate_min (ctx, p, pop);
}

static const struct evalconf
{
    int    maxormin;
    double (*f)(PGAContext *, int, int);
} evalconf [] =
{ { PGA_MINIMIZE, evaluate_min }
, { PGA_MAXIMIZE, evaluate_max }
};
int nevaltyp = sizeof (evalconf) / sizeof (evalconf [0]);

int main( int argc, char **argv )
{
    PGAContext *ctx;
    int i, sel, fit, fmin;
    int rank;

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    for (i=0; i<nevaltyp; i++) {
        int is_min = evalconf [i].maxormin == PGA_MINIMIZE;
        for (sel=0; sel<nselect; sel++) {
            int is_sus_or_prop =
                (  selection_type [sel].conf == PGA_SELECT_SUS
                || selection_type [sel].conf == PGA_SELECT_PROPORTIONAL
                );
            int nfit = is_sus_or_prop ? nfitness : 1;
            for (fit=0; fit<nfit; fit++) {
                int fminf = (is_min && is_sus_or_prop) ? nfitmin : 1;
                for (fmin=0; fmin<fminf; fmin++) {
                    if (rank == 0) {
                        printf ("* %s\n", is_min ? "Minimize" : "Maximize");
                        printf ("* %s Selection\n", selection_type [sel].name);
                        if (is_sus_or_prop) {
                            printf ("* %s\n", fitness_type [fit].name);
                            if (is_min) {
                                printf ("* %s\n", fitness_min [fmin].name);
                            }
                        }
                    }
                    ctx = PGACreate
                        ( &argc, argv
                        , PGA_DATATYPE_REAL
                        , dimension
                        , evalconf [i].maxormin
                        );

                    if (is_min && is_sus_or_prop) {
                        PGASetFitnessMinType (ctx, fitness_min [fmin].conf);
                    }
                    PGASetRandomSeed        (ctx, 1);
                    PGASetRealInitRange     (ctx, lower, upper);
                    PGASetMaxGAIterValue    (ctx, 100);
		    PGASetPopSize           (ctx, 30);
		    PGASetNumReplaceValue   (ctx, 10);
		    PGASetMutationType      (ctx, PGA_MUTATION_GAUSSIAN);
		    PGASetMutationOnlyFlag  (ctx, PGA_TRUE);
                    PGASetMutationRealValue (ctx,.5);
                    PGASetNoDuplicatesFlag  (ctx, PGA_TRUE);
                    PGASetSelectType        (ctx, selection_type [sel].conf);
                    PGASetFitnessType       (ctx, fitness_type [fit].conf);
                    PGASetUp                (ctx);
                    PGARun                  (ctx, evalconf [i].f);
                    PGADestroy              (ctx);
                }
            }
        }
    }
    
    MPI_Finalize();
    return 0;
}
