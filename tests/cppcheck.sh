#!/bin/sh
cppcheck -j4 -f -q --enable=style \
    -I${srcdir}/.. -I${srcdir}/../src -I/usr/include -I/usr/local/include -I/usr/lib/gcc/*/*/include \
    ${srcdir}/../src 2> cppcheck.log

if [ "$(wc -l cppcheck.log | sed -n 's|^\([0-9]*\).*$|\1|p')" -gt 0 ]; then
    cat cppcheck.log
    exit 1
else
    exit 0
fi
