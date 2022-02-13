#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include "pgapack.h"

double compute_asf (PGAContext *ctx, double *point, int axis);
void compute_utopian (PGAContext *ctx, PGAIndividual **start, int n);
void compute_extreme (PGAContext *ctx, PGAIndividual **start, int n);
int compute_intersect (PGAContext *ctx, PGAIndividual **start, int n);
void compute_worst
    (PGAContext *ctx, PGAIndividual **start, int n, double *wpop, double *wof);
void compute_nadir (PGAContext *ctx, PGAIndividual **start, int n);
int ranking (PGAContext *ctx, PGAIndividual **start, int n, int goal);

double pop [][3] =
{{0.59200264, 1.25973595, 1.06928177}
,{0.33668049, 1.73701836, 0.58860729}
,{0.06632383, 0.39695692, 2.04663035}
,{0.25462559, 0.95805078, 1.84763014}
,{1.56691598, 0.32456604, 0.77699699}
,{1.39468614, 1.06112185, 0.28415417}
,{0.49165715, 1.57516065, 0.36687402}
,{0.25219645, 1.52832133, 0.69918447}
,{0.19350479, 0.24441518, 2.03823009}
,{0.00421519, 1.0023554 , 1.89960353}
,{1.9530772 , 0.0804578 , 0.06106157}
,{1.22995435, 0.03603701, 1.49230785}
,{1.94390398, 0.33138624, 0.13822646}
,{1.20384517, 1.55884737, 0.62435854}
,{0.59739374, 0.74291846, 1.33993911}
,{1.00879625, 0.34569624, 1.62184635}
,{0.32193694, 0.86087728, 2.07955597}
,{1.16748989, 1.11720521, 1.13571308}
,{0.16932101, 2.06714265, 0.07015315}
,{0.05331214, 0.35609524, 1.44694504}
,{0.30808162, 1.81378312, 0.99359522}
,{0.05091241, 1.23561754, 0.83076159}
,{0.25347899, 0.14074459, 1.34134605}
,{1.32704271, 1.05191765, 0.54601349}
,{0.46011744, 1.03260765, 1.42242991}
,{0.09326116, 0.46870793, 1.56971846}
,{1.51139799, 0.2659622 , 0.65585821}
,{0.87693924, 0.81612064, 1.25916094}
,{0.3334751 , 1.64835976, 0.93961617}
,{0.52438602, 0.02850446, 1.61295935}
,{1.76450275, 0.40316984, 0.92311876}
,{0.70738245, 1.67520806, 0.46809858}
,{1.54665274, 0.67557381, 0.69206195}
,{0.60142408, 1.45105265, 1.0711516 }
,{1.27369097, 1.01588981, 0.20965035}
,{1.8724852 , 0.2317311 , 0.43240781}
,{1.18106135, 1.08833986, 1.49905451}
,{0.02402091, 0.16046297, 2.11222513}
,{0.23147743, 0.16866407, 1.99651279}
,{1.19388453, 0.64289518, 1.50748858}
,{0.29738157, 1.32541357, 0.76349803}
,{1.55812855, 0.63794389, 1.33363254}
,{0.75230374, 0.57150114, 1.33130107}
,{1.53830223, 0.43311618, 1.6965542 }
,{1.28723505, 0.86137531, 1.33015377}
,{0.53359748, 1.18223901, 1.26516872}
,{0.55409408, 0.49449924, 1.8326513 }
,{0.40382719, 0.75353682, 1.51244444}
,{0.91001828, 1.1018896 , 0.65642448}
,{0.64773201, 1.69230677, 0.63216507}
,{0.28587971, 0.87879756, 1.64750156}
,{0.61210389, 0.57125203, 1.68256191}
,{0.59686865, 0.20512127, 1.55331987}
,{0.44216441, 0.13774006, 1.74968488}
,{1.37173124, 0.59756938, 1.53394827}
,{0.51538914, 1.19139214, 1.41718206}
,{0.6174288 , 1.55052357, 1.22549865}
,{0.55592539, 1.47347103, 1.05184097}
,{0.22494409, 1.21045483, 0.63243309}
,{0.63462193, 1.13308044, 0.78382333}
,{0.18563641, 1.29432024, 1.05764917}
,{0.49426461, 0.20622366, 2.19938632}
,{0.5185551 , 1.81387573, 0.02649165}
,{0.7309261 , 1.34686213, 1.18318878}
,{0.34320965, 0.35092957, 1.46805473}
,{0.05038524, 1.29390942, 0.96247625}
,{0.68962815, 1.3397539 , 0.86034788}
,{0.16956981, 1.37471672, 1.16848244}
,{1.28338372, 0.29661796, 0.78446261}
,{0.36740655, 0.74203689, 1.26766327}
,{1.33692614, 0.79758974, 0.40364995}
,{0.86546472, 0.17180019, 1.35935998}
,{0.20350703, 0.49376692, 2.22405968}
,{1.41688725, 0.13626354, 0.68248411}
,{0.7189143 , 0.44157166, 1.3491942 }
,{0.28365323, 2.01618638, 0.82661351}
,{0.02597332, 0.86652861, 1.32409954}
,{1.05212941, 0.89540367, 1.28523146}
,{0.85523304, 1.51616613, 0.30206971}
,{1.45396861, 1.051698  , 0.04712763}
,{0.00523844, 0.01969276, 1.98367522}
,{0.25657941, 1.31180652, 1.12945444}
,{0.98498806, 1.23558789, 1.09579641}
,{0.76633323, 1.33764359, 0.82152697}
,{0.82402416, 1.57936609, 0.10396752}
,{0.82595538, 1.20499078, 1.51732213}
,{1.66409283, 1.29775787, 0.21067013}
,{0.10414121, 0.15014165, 2.06789412}
,{0.4136045 , 0.92521967, 1.417854  }
,{0.28818922, 0.22475283, 1.80305745}
,{0.31148591, 1.68686981, 0.81425242}
,{0.97177884, 1.31839026, 0.98843984}
,{0.13679796, 0.18122792, 1.74460377}
,{1.96718331, 0.63979882, 0.11732475}
,{0.98873743, 0.4355161 , 1.42898936}
,{1.03291285, 0.97516195, 0.17959043}
,{0.23229372, 0.08935152, 1.63973078}
,{0.45227308, 1.56357525, 0.51152652}
,{0.37152299, 1.48067771, 0.87198931}
,{0.83534425, 1.53702138, 0.99724714}};

double pop2 [][3] =
{{1.75521090e+000, 1.56051315e-014, 2.85987040e-038}
,{1.86468364e+000, 6.63720331e-006, 3.35671457e-069}
,{2.08582584e+000, 4.77101511e-005, 6.09673412e-006}
,{2.09676718e+000, 4.64348421e-008, 1.51978772e-016}
,{1.77884601e+000, 7.08170085e-089, 2.24706475e-054}
,{1.77535140e+000, 1.41680851e-038, 2.80272957e-099}
,{1.69040065e+000, 1.35679883e-009, 6.48404988e-086}
,{1.69947876e+000, 4.48457393e-005, 3.55804871e-057}
,{2.06193229e+000, 2.38411042e-024, 1.25146404e-004}
,{7.75432984e-001, 2.00298062e+000, 2.94745937e-016}
,{1.95568722e+000, 2.16648944e-158, 2.13476299e-170}
,{1.93418435e+000, 3.50097313e-173, 2.40570964e-025}
,{1.97678681e+000, 4.27124281e-097, 2.39156288e-135}
,{2.06617325e+000, 9.03069541e-024, 4.07656624e-071}
,{1.64445856e+000, 8.16636464e-025, 4.81958816e-022}
,{1.94102069e+000, 5.52142357e-068, 2.50453317e-020}
,{2.27361081e+000, 2.11144643e-011, 1.53273598e-013}
,{1.97510110e+000, 1.42725682e-031, 3.97218008e-041}
,{2.07519324e+000, 1.55831633e-002, 6.41377605e-167}
,{1.49107208e+000, 1.13057423e-004, 1.10062374e-007}
,{2.09092203e+000, 3.94630645e-005, 2.40575446e-050}
,{1.48075671e+000, 1.63912075e-001, 9.01591074e-043}
,{1.37232279e+000, 1.62523118e-049, 1.02103949e-006}
,{1.77924249e+000, 2.87555736e-037, 1.72859876e-070}
,{1.81694620e+000, 9.42547087e-014, 1.70663229e-024}
,{1.64085368e+000, 4.07465922e-006, 2.28683332e-009}
,{1.66889478e+000, 8.09648668e-096, 2.70300333e-059}
,{1.73797629e+000, 1.99910768e-032, 4.88762601e-029}
,{1.92644078e+000, 3.78824598e-006, 3.83648985e-049}
,{1.69629922e+000, 1.98195031e-146, 5.17284901e-010}
,{2.03178841e+000, 1.09489412e-084, 1.78673352e-052}
,{1.87771889e+000, 5.27845376e-013, 9.74723976e-080}
,{1.82413936e+000, 2.07008479e-058, 7.17045741e-061}
,{1.90121815e+000, 9.39880727e-013, 3.71096885e-042}
,{1.64264240e+000, 4.15707075e-037, 3.26248146e-109}
,{1.93568506e+000, 8.07799257e-111, 1.39588474e-084}
,{2.19694197e+000, 1.30656913e-032, 3.07125118e-032}
,{2.11832985e+000, 1.60790369e-004, 2.23372643e-002}
,{2.01695133e+000, 6.33250086e-040, 2.35076840e-004}
,{2.02760852e+000, 1.82227234e-050, 1.69444540e-027}
,{1.55823174e+000, 6.49830992e-007, 5.11864066e-049}
,{2.14786237e+000, 7.36885653e-061, 3.29202083e-037}
,{1.63246348e+000, 1.16302718e-038, 5.42711479e-022}
,{2.33072082e+000, 6.27927563e-076, 1.20693374e-028}
,{2.04162938e+000, 9.11073096e-043, 9.86944793e-035}
,{1.81192362e+000, 6.19400859e-014, 4.54054083e-031}
,{1.97741259e+000, 1.35114529e-033, 1.90839482e-012}
,{1.73734922e+000, 1.31690927e-016, 1.59693934e-017}
,{1.57263697e+000, 1.78909647e-025, 1.53739354e-056}
,{1.91913825e+000, 9.43776419e-012, 2.87279108e-067}
,{1.88898750e+000, 5.87780083e-010, 2.36955342e-017}
,{1.87936549e+000, 2.60474509e-032, 2.25833230e-015}
,{1.67664234e+000, 6.21856252e-068, 1.49731676e-012}
,{1.80993891e+000, 6.92504248e-072, 4.33067356e-008}
,{2.14283295e+000, 1.91648848e-058, 1.27874796e-029}
,{1.92183407e+000, 2.55598551e-013, 5.43662354e-028}
,{2.07055273e+000, 3.32583667e-012, 1.16536924e-039}
,{1.89381608e+000, 1.38266643e-011, 7.28139954e-043}
,{1.38411428e+000, 8.60526653e-006, 2.24850805e-052}
,{1.51690319e+000, 2.02128728e-017, 1.76549952e-046}
,{1.68176916e+000, 1.96416498e-004, 1.17171240e-036}
,{2.26365322e+000, 4.25352908e-060, 2.44011859e-007}
,{1.88672901e+000, 9.94404451e-009, 3.99227539e-205}
,{1.93603364e+000, 9.00853923e-017, 4.57799492e-038}
,{1.54794351e+000, 7.82558080e-030, 2.50927239e-010}
,{1.60025772e+000, 2.05611317e-001, 2.26454930e-039}
,{1.73514436e+000, 6.04938926e-016, 2.10513648e-048}
,{1.81216735e+000, 8.34077865e-004, 2.50736091e-035}
,{1.53311369e+000, 2.49941861e-084, 6.00226642e-047}
,{1.51412565e+000, 2.17604075e-015, 2.65277197e-020}
,{1.60824568e+000, 7.26110004e-047, 1.66963356e-079}
,{1.62061841e+000, 1.02402054e-090, 3.79036447e-020}
,{2.28728274e+000, 1.33792722e-012, 3.13179823e-007}
,{1.57858221e+000, 8.98417405e-122, 6.63010952e-055}
,{1.59127254e+000, 7.69226283e-046, 2.00491099e-019}
,{2.19744316e+000, 3.09593915e-004, 3.51958627e-061}
,{1.54137142e+000, 3.59110264e-001, 2.48675463e-020}
,{1.88694037e+000, 4.83939502e-035, 2.11299779e-032}
,{1.76675675e+000, 1.76787784e-017, 2.17894169e-096}
,{1.79508061e+000, 3.23867721e-040, 5.78592045e-178}
,{1.36055481e+000, 2.96324918e-008, 1.44370122e+000}
,{1.74995333e+000, 5.50565026e-006, 2.71910577e-035}
,{1.92292711e+000, 1.53075176e-024, 1.37804834e-041}
,{1.74684388e+000, 9.28495950e-018, 6.49322907e-051}
,{1.78443893e+000, 3.75421843e-016, 2.51093298e-143}
,{2.10629331e+000, 3.80708668e-021, 2.83016047e-029}
,{2.12079286e+000, 1.04385660e-037, 4.92929710e-120}
,{2.07592663e+000, 2.11592659e-021, 1.01296568e-002}
,{1.74281671e+000, 8.14173036e-014, 4.05190612e-022}
,{1.83972362e+000, 9.07244132e-038, 3.52112593e-006}
,{1.89883128e+000, 1.28179808e-005, 3.31884942e-055}
,{1.91298727e+000, 9.44345491e-023, 2.21390229e-046}
,{1.75931785e+000, 2.54870391e-023, 5.09454165e-004}
,{2.07193576e+000, 4.51880720e-070, 1.67911730e-144}
,{1.79144818e+000, 4.28898118e-058, 2.38542190e-023}
,{1.43181791e+000, 4.25928969e-032, 4.94420389e-110}
,{1.65851171e+000, 1.97003171e-063, 1.09050780e-004}
,{1.70615881e+000, 7.06387810e-009, 1.49354978e-071}
,{1.75806745e+000, 1.12037051e-007, 2.20838283e-048}
,{2.01363765e+000, 8.84911811e-017, 2.14305904e-048}};

double ideal [] = {0.00421519, 0.01969276, 0.02649165};
double nadir [] = {1.96718331, 2.06714265, 1.98465785};
size_t npop = sizeof (pop) / 3 / sizeof (double);

void test_pop (PGAContext *ctx, void *p)
{
    int i, j;
    double (*pop) [3] = p;
    double (*extreme) [3];
    double wpop [3];
    double wof [3];
    PGAIndividual *start [npop];

    assert (sizeof (pop2) / 3 / sizeof (double) == npop);

    ctx->ga.extreme_valid = ctx->ga.utopian_valid = PGA_FALSE;
    ctx->ga.worst_valid = PGA_FALSE;
    for (i=0; i<npop; i++) {
        PGAIndividual *ind = ctx->ga.newpop + i;
        start [i] = ind;
        for (j=0; j<3; j++) {
            if (j==0) {
                ind->evalue = pop [i][j];
            } else {
                ind->auxeval [j-1] = pop [i][j];
            }
        }
    }
    ranking (ctx, start, npop, npop);
    compute_utopian   (ctx, start, npop);
    compute_extreme   (ctx, start, npop);
    compute_intersect (ctx, start, npop);
    printf ("Utopian: ");
    for (j=0; j<3; j++) {
        printf ("%e ", ctx->ga.utopian [j]);
    }
    printf ("\nIntersect: ");
    for (j=0; j<3; j++) {
        printf ("%e ", ctx->ga.nadir [j]);
    }
    compute_nadir   (ctx, start, npop);
    printf ("\nNadir (norm): ");
    for (j=0; j<3; j++) {
        printf ("%e ", ctx->ga.nadir [j] - ctx->ga.utopian [j]);
    }
    printf ("\n");
    printf ("Nadir: ");
    for (j=0; j<3; j++) {
        printf ("%e ", ctx->ga.nadir [j]);
    }
    printf ("\n");
    ctx->ga.worst_valid = PGA_FALSE;
    compute_worst (ctx, start, npop, wpop, wof);
    printf ("worst: ");
    for (j=0; j<3; j++) {
        printf ("%e ", ctx->ga.worst [j]);
    }
    printf ("\n");
    printf ("wpop: ");
    for (j=0; j<3; j++) {
        printf ("%e ", wpop [j]);
    }
    printf ("\n");
    printf ("Worst of front: ");
    for (j=0; j<3; j++) {
        printf ("%e ", wof [j]);
    }
    printf ("\n");
    printf ("Extreme:\n");
    extreme = ctx->ga.extreme;
    for (i=0; i<3; i++) {
        printf ("[ ");
        for (j=0; j<3; j++) {
            printf ("%e ", extreme [i][j]);
        }
        printf ("]\n");
    }
}

int main (int argc, char **argv)
{
    PGAContext *ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, npop, PGA_MINIMIZE);
    /* Example from slides EMO '19 */
    int i;
    double utop [] = {0.1, 0.1, 0.1};
    double f [][3] =
        { { 1.0, 0.1, 0.2 }
        , { 0.2, 1.0, 0.1 }
        , { 0.2, 0.5, 1.0 }
        , { 0.1, 0.9, 0.9 }
        , { 0.4, 0.4, 0.9 }
        , { 0.3, 0.3, 100 }
        };
    int dl = sizeof (f) / (3 * sizeof (double));
    printf ("%zu\n", npop);
    ctx->ga.NumAuxEval = 2;
    ctx->ga.NumConstraint = 0;
    ctx->ga.utopian = utop;
    for (i=0; i<dl; i++) {
        printf ("%e\n", compute_asf (ctx, f [i], 2));
    }
    PGASetPopReplaceType (ctx, PGA_POPREPL_NSGA_III);
    PGASetPopSize (ctx, 100);
    PGASetUp (ctx);
    test_pop (ctx, pop);
    test_pop (ctx, pop2);
    {
        double minasf = -1;
        double d [3];
        int asfidx = -1;
        for (i=0; i<npop; i++) {
            int j;
            double asf;
            memcpy (d, pop2 [i], 3 * sizeof (double));
            for (j=0; j<3; j++) {
                if (d [j] < 1e-3) {
                    d [j] = 0;
                }
            }
            asf = compute_asf (ctx, d, 0);
            if (asfidx < 0 || asf < minasf) {
                minasf = asf;
                asfidx = i;
            }
        }
        printf ("ASF22: %e\n", compute_asf (ctx, pop2 [22], 0));
        printf ("ASF: min: %e idx: %d\n", minasf, asfidx);
    }
}