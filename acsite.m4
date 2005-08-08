AC_DEFUN([AX_GCC_OPTION], [
AC_REQUIRE([AC_PROG_CC])
if test "x$GCC" = "xyes"; then
        AC_MSG_CHECKING([if gcc accepts $2 option])
        if AC_TRY_COMMAND($CC $2) >/dev/null 2>&1; then
                $1=$3
                AC_MSG_RESULT([yes])
        else
                $1=$4
                AC_MSG_RESULT([no])
        fi
else
        unset $1
        AC_MSG_RESULT([sorry, no gcc available])
fi
])
