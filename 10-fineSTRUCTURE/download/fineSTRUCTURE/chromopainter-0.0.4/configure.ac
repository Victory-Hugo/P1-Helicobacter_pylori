AC_INIT([chromopainter C++], [0.0.1], [ghellenthal@gmail.com],
             [chromopainter], [http://www.paintmychromosomes.com])
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([1.10 -Wall no-define])
AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# find gsl, wherever it's hiding on the system
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
#AC_CHECK_LIB([zlib])
#AX_PATH_GSL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
