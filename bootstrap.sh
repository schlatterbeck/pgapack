#!/bin/sh

set -x
if [ ! -d config ]; then
    mkdir config
fi

aclocal -I config
autoheader
automake --add-missing --copy
autoconf

