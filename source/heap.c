/*
 *  
 *  ********************************************************************* 
 *  (C) COPYRIGHT 1995 UNIVERSITY OF CHICAGO 
 *  *********************************************************************
 *  
 *  This software was authored by
 *  
 *  D. Levine
 *  Mathematics and Computer Science Division Argonne National Laboratory
 *  Argonne IL 60439
 *  levine@mcs.anl.gov
 *  (708) 252-6735
 *  (708) 252-5986 (FAX)
 *  
 *  with programming assistance of participants in Argonne National 
 *  Laboratory's SERS program.
 *  
 *  This program contains material protectable under copyright laws of the 
 *  United States.  Permission is hereby granted to use it, reproduce it, 
 *  to translate it into another language, and to redistribute it to 
 *  others at no charge except a fee for transferring a copy, provided 
 *  that you conspicuously and appropriately publish on each copy the 
 *  University of Chicago's copyright notice, and the disclaimer of 
 *  warranty and Government license included below.  Further, permission 
 *  is hereby granted, subject to the same provisions, to modify a copy or 
 *  copies or any portion of it, and to distribute to others at no charge 
 *  materials containing or derived from the material.
 *  
 *  The developers of the software ask that you acknowledge its use in any 
 *  document referencing work based on the  program, such as published 
 *  research.  Also, they ask that you supply to Argonne National 
 *  Laboratory a copy of any published research referencing work based on 
 *  the software.
 *  
 *  Any entity desiring permission for further use must contact:
 *  
 *  J. Gleeson
 *  Industrial Technology Development Center Argonne National Laboratory
 *  Argonne IL 60439
 *  gleesonj@smtplink.eid.anl.gov
 *  (708) 252-6055
 *  
 *  ******************************************************************** 
 *  DISCLAIMER
 *  
 *  THIS PROGRAM WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY AN AGENCY 
 *  OF THE UNITED STATES GOVERNMENT.  NEITHER THE UNIVERSITY OF CHICAGO, 
 *  THE UNITED STATES GOVERNMENT NOR ANY OF THEIR EMPLOYEES MAKE ANY 
 *  WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR 
 *  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY 
 *  INFORMATION OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT 
 *  INFRINGE PRIVATELY OWNED RIGHTS.
 *  
 *  ********************************************************************** 
 *  GOVERNMENT LICENSE
 *  
 *  The Government is granted for itself and others acting on its behalf a 
 *  paid-up, non-exclusive, irrevocable worldwide license in this computer 
 *  software to reproduce, prepare derivative works, and perform publicly 
 *  and display publicly.
 */

/*****************************************************************************
*      FILE: heap.c: This file contains routines for sorting individuals for
*                     selection
*
*      Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*               Brian P. Walenz
*****************************************************************************/

#include <pgapack.h>

/******************************************************************************
   PGAAdjustHeap - Auxiliary routine called by PGA*HeapSort

   Category: Sorting

   Inputs:
       ctx      - context variable
       a        - array of values to be sorted
       idx      - array of integer indices corresponding to the array
                  a being sorted
       i        - point of combination -- combine the node at a[i] with
                  the two heaps at a[2i+1] and a[2i+2] to form a single
                  heap.  0 <= i <= n-1.
       n        - size of the arrays a and idx
       j        - temporary variable, integer
       item     - temporary variable, type must be same as a
       item_idx - temporary variable, integer

   Output:

   Example:

******************************************************************************/
#define PGAAdjustHeap(ctx, a, idx, i, n, j, item, item_idx) {       \
  item     = a[i];                                                  \
  item_idx = idx[i];                                                \
  j = 2*i+1;      /* let j be the left child */                     \
  while (j < n) {                                                   \
    if (j<n-1 && a[j] > a[j+1])                                     \
       j = j + 1;       /* j is the larger child */                 \
    if (item <= a[j])   /* a position for item has been found */    \
       break;                                                       \
    a[(j-1)/2]   = a[j];    /* move the larger child up a level */  \
    idx[(j-1)/2] = idx[j];                                          \
    j = j*2+1;                                                      \
  }                                                                 \
  a[(j-1)/2]   = item;                                              \
  idx[(j-1)/2] = item_idx;                                          \
}                                                                   \




/*I****************************************************************************
   PGADblHeapSort - Uses a heapsort algorithm to sort from largest to smallest
   element.  An integer array, intialized with the original indices of the
   elements of array a is sorted also so that the original locations are known

   Category: Sorting

   Inputs:
       ctx      - context variable
       a        - array of (double) values to be sorted
       idx      - array of integer indices corresponding to the array
                  a being sorted
       n        - size of the arrays a and idx

   Output:
       The sorted arrays a and idx

   Example:
      The following code sorts the population by fitness

      PGAContext *ctx;
      int i,j,n,idx[LARGE]
      double a[LARGE];
      :
      n = PGAGetPopsize(ctx);
      for(i=0;i<n;i++) {
        a[i]   = PGAGetFitness(ctx,p,PGA_OLDPOP);
        idx[i] = i;
      }
      PGADblHeapSort ( ctx, a, idx, n);

****************************************************************************I*/
void PGADblHeapSort ( PGAContext *ctx, double *a, int *idx, int n )
{
  int i;
  double temp_a;
  int temp_idx;
  int j, item_idx;
  double item;

    PGADebugEntered("PGADblHeapSort");

  /*  Create a heap from our array  */
  for (i=(n-2)/2; i>=0; i--)
    PGAAdjustHeap(ctx, a, idx, i, n, j, item, item_idx);

  for ( i=n-1; i>=1; i--)  /* interchange the new maximum with the   */
  {                        /* element at the end of the tree         */
    temp_a   = a[i];
    temp_idx = idx[i];
    a[i]     = a[0];
    idx[i]   = idx[0];
    a[0]     = temp_a;
    idx[0]     = temp_idx;
    PGAAdjustHeap(ctx, a, idx, 0, i, j, item, item_idx);
  }

    PGADebugExited("PGADblHeapSort");
}




/*I****************************************************************************
   PGAIntHeapSort - Uses a heapsort algorithm to sort from largest to smallest
   element.  An integer array, intialized with the original indices of the
   elements of array a is sorted also so that the original locations are known

   Category: Sorting

   Inputs:
       ctx      - context variable
       a        - array of (int) values to be sorted
       idx      - array of integer indices corresponding to the array
                  a being sorted
       n        - size of the arrays a and idx

   Output:
       The sorted arrays a and idx

   Example:
      The following code sorts the population by fitness

      PGAContext *ctx;
      int i,j,n,idx[LARGE],a[LARGE];
      :
      n = PGAGetPopsize(ctx);
      for(i=0;i<n;i++) {
        a[i]   = (int) PGAGetEvaluation(ctx,p,PGA_OLDPOP);
        idx[i] = i;
      }
      PGAIntHeapSort ( ctx, a, idx, n);

****************************************************************************I*/
void PGAIntHeapSort ( PGAContext *ctx, int *a, int *idx, int n )
{
  int i;                   /* index of for loops                      */
  int temp_a;
  int temp_idx;
  int j, item_idx;
  double item;

    PGADebugEntered("PGAIntHeapSort");

  /*  Create a heap from our elements.  */
  for (i=(n-2)/2; i>=0; i--)
    PGAAdjustHeap(ctx, a, idx, i, n, j, item, item_idx);

  for ( i=n-1; i>=1; i--)  /* interchange the new maximum with the   */
  {                        /* element at the end of the tree         */
    temp_a   = a[i];
    temp_idx = idx[i];
    a[i]     = a[0];
    idx[i]   = idx[0];
    a[0]     = temp_a;
    idx[0]     = temp_idx;
    PGAAdjustHeap(ctx, a, idx, 0, i, j, item, item_idx);
  }

    PGADebugExited("PGAIntHeapSort");
}


