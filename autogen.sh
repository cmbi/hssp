#!/bin/sh

touch NEWS README

aclocal \
&& autoheader \
&& automake --add-missing \
&& autoconf
