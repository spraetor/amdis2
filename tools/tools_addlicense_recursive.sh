#!/bin/bash

export LICENSELEN=`wc -l license.templ.txt | cut -f1 -d ' '`  
find .. -type f \( -name *.hpp -o -name *.h -o -name *.hh -o -name *.cc \) -print0 | xargs -0 ./addlicense.sh  
