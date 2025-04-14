      include 'pgapackf.h'
      include 'mpif.h'

      double precision evaluate
      external         evaluate
      
      integer,parameter :: k18 = selected_int_kind(18)
      integer(kind=k18) ctx
      integer           rank
      integer           i, ierror
      double precision  lower(2), upper(2)

      data (lower(i), i=1,2)     / 0.1D0, 0.0D0 /
      data (upper(i), i=1,2)     / 1.0D0, 5.0D0 /

      call MPI_Init(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      if (rank .eq. 0) then
          print *, 'Example: CONSTR (Fortran)'
      endif

      ctx = PGACreate (PGA_DATATYPE_REAL, 2, PGA_MINIMIZE)
      call PGASetRandomSeed             (ctx, 1)
      call PGASetPopSize                (ctx, 100)
      call PGASetNumReplaceValue        (ctx, 100)
      call PGASetSelectType             (ctx, PGA_SELECT_LINEAR)
      call PGASetPopReplaceType         (ctx, PGA_POPREPL_NSGA_II)
      call PGASetMutationOnlyFlag       (ctx, PGA_TRUE)
      call PGASetMutationType           (ctx, PGA_MUTATION_DE)
      call PGASetDECrossoverProb        (ctx, 0.8D0)
      call PGASetDECrossoverType        (ctx, PGA_DE_CROSSOVER_BIN)
      call PGASetDEVariant              (ctx, PGA_DE_VARIANT_RAND)
      call PGASetDEScaleFactor          (ctx, 0.85D0)
      call PGASetRealInitRange          (ctx, lower, upper)
      call PGASetMaxGAIterValue         (ctx, 250)
      call PGASetNumAuxEval             (ctx, 3)
      call PGASetNumConstraint          (ctx, 2)
      call PGASetNoDuplicatesFlag       (ctx, PGA_TRUE)
      call PGASetMutationBounceBackFlag (ctx, PGA_TRUE)

      call PGASetUp                     (ctx)
      call PGARun                       (ctx, evaluate)
      call PGADestroy                   (ctx)

      call MPI_Finalize(ierror)

      stop
      end

      double precision function evaluate(ctx, p, pop, aux)
      include          'pgapackf.h'
      integer,parameter :: k18 = selected_int_kind(18)
      integer(kind=k18) ctx
      integer           p, pop
      double precision  aux(*)
      double precision  x1, x2

      x1 = PGAGetRealAllele (ctx, p, pop, 1)
      x2 = PGAGetRealAllele (ctx, p, pop, 2)
      aux(1) = (1 + x2) / x1
      aux(2) = 6 - (x2 + 9*x1)
      aux(3) = 1 + x2 - 9*x1
      evaluate = x1
      return
      end
