#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([PGAPACK], [1.1], [BUG-REPORT-ADDRESS])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([source/create.c])
AC_CONFIG_HEADER([include/config.h])

# Checks for programs.
AC_LANG([C])
AC_PROG_CC
AC_PROG_F77()
AM_CONDITIONAL([HAVE_FORTRAN],[test -n $ac_cv_prog_F77])
AC_ARG_ENABLE([fortran],
     [  --enable-fortran    Turn on fortran tests],
     [case "${enableval}" in
       yes) fortran=true ;;
       no)  fortran=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-fortran]) ;;
     esac],[fortran=false])
AM_CONDITIONAL([COMPILEFORTRAN], [test x$fortran = xtrue])
AC_F77_LIBRARY_LDFLAGS
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.

AC_SUBST(LIBPGAPACK)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(source/Makefile)
AC_CONFIG_FILES(test/Makefile)
AC_CONFIG_FILES(examples/Makefile)
AC_CONFIG_FILES(examples/c/Makefile)
AC_CONFIG_FILES(examples/fortran/Makefile)

AC_OUTPUT

