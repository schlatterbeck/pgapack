#include <stdio.h>
#include "pgapack.h"


int main (int argc, char **argv)
{
    char *s1 = "abcdefghi";
    char *s2 = "abcdefgh";

    printf ("s1: %04x\n", PGAUtilHash (s1, strlen (s1)));
    printf ("s2: %04x\n", PGAUtilHash (s2, strlen (s2)));
}
