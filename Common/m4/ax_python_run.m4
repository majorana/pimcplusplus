# Copyright © 2005 Robert White <kranki@mac.com>
# Copyright © 2005 Dustin Mitchell <dustin@cs.uchicago.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. 

# AX_PYTHON_RUN( PYTHON_PROGRAM )
# -----------------
# Run a Python Test Program saving its output
# in ax_python_output and its condition code
# in ax_python_cc.

AC_DEFUN([AX_PYTHON_RUN],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        cat >conftest.py <<_ACEOF
$1
_ACEOF
        ax_python_output=`$PYTHON conftest.py`
        ax_python_cc=$?
        rm conftest.py
        if test -f "conftest.pyc"
        then
            rm conftest.pyc
        fi
    fi
])