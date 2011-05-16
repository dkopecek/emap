#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.60])
AC_INIT([emap], [0.1], [xkopecek@fi.muni.cz])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign tar-pax])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

CFLAGS="$CFLAGS -pipe -std=c99 -W -Wall -Wnonnull -Wno-unknown-pragmas -D_GNU_SOURCE"

CFLAGS_DEBUGGING="-fno-inline-functions -O0 -g3 -fstack-protector -fstack-check"
CFLAGS_OPTIMIZED="-Ofast -ffast-math -fomit-frame-pointer -DNDEBUG -mtune=native"

# Some parts of the code pass small integer values thru a pointer
# argument, so all warnings muted by the options specified bellow
# are expected. These options aren't used if the sources are
# ./configured with `--enable-debug'.
CFLAGS_NONDEBUG="-Wno-pointer-to-int-cast -Wno-int-to-pointer-cast"

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_ARG_ENABLE([debug],
     [AC_HELP_STRING([--enable-debug], [enable debugging flags (default=no)])],
     [case "${enableval}" in
       yes) debug=yes ;;
       no)  debug=no ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
     esac], [debug=no])

if test "$debug" = "yes"; then
   CFLAGS="$CFLAGS $CFLAGS_DEBUGGING"
else
   CFLAGS="$CFLAGS $CFLAGS_NONDEBUG $CFLAGS_OPTIMIZED"
fi

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 src/Makefile
                 tests/Makefile])

AC_OUTPUT

echo "CFLAGS: $CFLAGS"
