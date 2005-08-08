AC_DEFUN([AX_CC_OPTION], [
AC_REQUIRE([AC_PROG_CC])
AC_MSG_CHECKING([if ${CC-cc} accepts $2 option])
echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $2 -c conftest.c 2>&1`"; then
      $1=$3
      AC_MSG_RESULT([yes])
else
      $1=$4
      AC_MSG_RESULT([no])
fi
rm -f conftest*
])
