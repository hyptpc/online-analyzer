#!/bin/sh
bin=${HOME}/online-v9/bin/beam
RUN=$(printf %05d $1)
#conf=${HOME}/k18br_analyzer/e73/param/conf/analyzer_e73.conf
conf=${HOME}/k18br_analyzer/e73/param/conf/analyzer.online
input=/group/had/knucl/e15/e73_data/Run85/run$RUN.dat.gz
#$bin $conf $input
#output=/home/oper/data/dorami/test_$RUN.root
output=./e73_$RUN.root
$bin $conf $input $output
