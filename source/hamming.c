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
*     FILE: hamming.c: This file contains the routines that have to do with
*                      Hamming distances.
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle
*****************************************************************************/

#include "pgapack.h"

/*U****************************************************************************
  PGAHammingDistance - Calculates the mean Hamming distance for a population
  of binary strings.  For all other data types returns a value of 0.0 and
  prints a warning message.

  Category: Utility

  Inputs:
      ctx      - context variable
      popindex - symbolic constant of the population for which the
                 Hamming distance is to be calculated
  Output:
      The mean Hamming distance in the population

  Example:
      PGAContext *ctx;
      double hd;
      :
      hd = PGAHammingDistance(ctx, PGA_NEWPOP);

****************************************************************************U*/
double PGAHammingDistance( PGAContext *ctx, int popindex)
{
    int i, j, hd, count=0;
    double avg_hd = 0.;
    PGAIndividual *pop;      /* pointer to appropriate population          */

    PGADebugEntered("PGAHammingDistance");

    switch (popindex) {
    case PGA_OLDPOP:
        pop = ctx->ga.oldpop;
        break;
    case PGA_NEWPOP:
        pop = ctx->ga.newpop;
        break;
    default:
        PGAError( ctx, "PGAHammingDistance: Invalid value of popindex:",
                  PGA_FATAL, PGA_INT, (void *) &popindex );
        break;
    }

    switch (ctx->ga.datatype) {
    case PGA_DATATYPE_BINARY:
        for(i=0; i<ctx->ga.PopSize-1; ++i)
            for ( j = i+1; j<ctx->ga.PopSize; ++j ) {
                count++;
                hd = PGABinaryHammingDistance( ctx,
                                            (pop+i)->chrom, (pop+j)->chrom );
                avg_hd += (double) hd;
            }
        avg_hd /= (double) count;
        break;
    case PGA_DATATYPE_INTEGER:
        avg_hd = 0.0;
        PGAError( ctx,
        "PGAHammingDistance: No Hamming Distance for PGA_DATATYPE_INTEGER ",
                  PGA_WARNING,
                  PGA_DOUBLE,
                  (void *) &avg_hd );
        break;
    case PGA_DATATYPE_REAL:
        avg_hd = 0;
        PGAError( ctx,
        "PGAHammingDistance: No Hamming Distance for PGA_DATATYPE_REAL ",
                  PGA_WARNING,
                  PGA_DOUBLE,
                  (void *) &avg_hd );
        break;
    case PGA_DATATYPE_CHARACTER:
        avg_hd = 0;
        PGAError( ctx,
        "PGAHammingDistance: No Hamming Distance for PGA_DATATYPE_CHARACTER ",
                  PGA_WARNING,
                  PGA_DOUBLE,
                  (void *) &avg_hd );
        break;
    case PGA_DATATYPE_USER:
        avg_hd = 0;
        PGAError( ctx,
        "PGAHammingDistance: No Hamming Distance for PGA_DATATYPE_USER ",
                  PGA_WARNING,
                  PGA_DOUBLE,
                  (void *) &avg_hd );
        break;
    default:
        PGAError( ctx,
                 "PGAHammingDistance: Invalid value of datatype:",
                  PGA_FATAL,
                  PGA_INT,
                  (void *) &(ctx->ga.datatype) );
        break;
    }

    PGADebugExited("PGAHammingDistance");

    return(avg_hd);
}
