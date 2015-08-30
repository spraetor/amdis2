#!/bin/bash

/opt/software/astyle/bin/astyle --style=break -s2 -Y -N -w -k1 -W1 -xy -r -n "../src/*.c*" "../src/*.h*"
