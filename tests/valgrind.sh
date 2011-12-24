#!/bin/sh
valgrind --log-file=valgrind.log \
         --leak-check=full       \
         --show-reachable=yes    \
         --read-var-info=yes     \
         --track-origins=yes     \
         -q \
         ${BUILDDIR}/src/emap -s 0.5 ${srcdir}/data/simple01.dat 2> vg-stderr.log 1> vg-stdout.log

if [ "$(wc -l valgrind.log | sed -n 's|^\([0-9]*\).*$|\1|p')" -gt 0 ]; then
    cat valgrind.log
    exit 1
else
    exit 0
fi
