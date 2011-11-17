

AC_INIT([fiter_fastq_by_bam], [1.0], [dcjones@cs.washington.edu])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])


AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AM_PROG_CC_C_O

opt_CFLAGS="-Wall -Wextra -g -O3"
dbg_CFLAGS="-Wall -Wextra -g -O0"


AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

AC_DEFINE(_FILE_OFFSET_BITS, 64)

AC_PROG_LIBTOOL

AC_CHECK_FUNC(fileno, ,
              AC_MSG_ERROR([The 'fileno' function is missing.]))

# check zlib
AX_CHECK_ZLIB

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/samtools/Makefile
                 src/hat-trie/Makefile])
AC_OUTPUT





