#!/bin/bash

for x in $*; do  
  head -$LICENSELEN $x | diff license.templ.txt - || ( ( cat license.templ.txt; echo; cat $x) > /tmp/file;  
  mv /tmp/file $x )  
done
