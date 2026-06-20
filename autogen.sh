#!/bin/sh
# Regenerate the Autotools build system from configure.ac / Makefile.am.
#
# Run this once after cloning from git, before ./configure.
# Requires autoconf (>= 2.57), automake, and the m4 macros under m4/.
# Release tarballs made with `make dist` already include configure, so
# tarball users do NOT need to run this.
#
# Note: no --force, so automake's --add-missing will NOT overwrite the
# hand-written, atlc-specific INSTALL (or COPYING) with generic versions.
set -e
autoreconf --install --warnings=all
echo "Bootstrap complete. Now run: ./configure && make"
