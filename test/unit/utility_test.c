#include <stdio.h>
#include "pgapack.h"


/* The algo does *not* produce the same hash when hashing two strings
 * after one another and the concatenated string in one go
 */
int main (int argc, char **argv)
{
    char *s1    = "abcdefghi";
    char *s2    = "abcdefgh";
    char *s3    = "abcdefghiabcdefgh";
    char *s12   = "abcdefghijkl";
    char *s13   = "abcdefghijklabcdefgh";
    char *s     = "This is the time for all good men to come to the aid "
                  "of their country";
    PGAHash h1  = PGAUtilHash (s1,  strlen (s1),  PGA_INITIAL_HASH);
    PGAHash h12 = PGAUtilHash (s12, strlen (s12), PGA_INITIAL_HASH);

    printf ("s1: %04x\n", h1);
    printf ("s2: %04x\n", PGAUtilHash (s2, strlen (s2), PGA_INITIAL_HASH));
    printf ("s3: %04x\n", PGAUtilHash (s3, strlen (s3), PGA_INITIAL_HASH));
    printf ("s1+s2: %04x\n", PGAUtilHash (s2, strlen (s2), h1));
    printf ("s12: %04x\n", h12);
    printf ("s12+s2: %04x\n", PGAUtilHash (s2,   strlen (s2),   h12));
    printf ("s13: %04x\n", PGAUtilHash (s13, strlen (s13), PGA_INITIAL_HASH));
    /* From the test uf lookup2: */
    printf ("s: %04x\n", PGAUtilHash (s, strlen (s), 0));
}
