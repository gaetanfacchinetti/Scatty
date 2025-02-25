#!/bin/bash
gcc -shared -o libscatty.so -fPIC integration.c cross_section.c -lgsl -lgslcblas -lm
echo "Shared library created: libscatty.so"