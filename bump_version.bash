#!/bin/bash

# Sets the version number for mkhssp and hsspconv

echo -e "#ifndef version_h\n#define version_h\n\n#define HSSP_VERSION \"$1\"\n\n#endif" > $(dirname $0)/src/version.h
