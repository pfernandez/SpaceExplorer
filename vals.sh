#!/bin/bash
#Pulls a value from a file. Arguments are [valueName] [file],
#i.e. `./vals.sh /path/to/radstarMinus fort.50`.
#Printed number must be in 2nd column; columns may be
#delimited by any number of spaces. 

grep -m 1 $1 $2 | sed -e s%${1}%%g -e s%[:=]%%g | awk '{ printf "%g", $1}'
