#!/bin/sh
ANAMAIN=$HOME/online-v9/src/
ln -s $PWD $ANAMAIN
cp -f Makefile.org $ANAMAIN/Makefile
cp -f common.mk $ANAMAIN/
cp -f $ANAMAIN/main/Makefile{.org,}
mkdir -p dict
