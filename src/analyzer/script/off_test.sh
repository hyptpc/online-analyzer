#!/bin/sh
bin=${HOME}/online-v9/bin/test
RUN=$(printf %05d $1)
conf=${HOME}/k18br_analyzer/e73/param/conf/analyzer.raw
input=/group/had/knucl/e15/t98/test/run$RUN.dat.gz
#$bin $conf $input                                                                                                                                                                                                                          
#output=/home/oper/data/dorami/test_$RUN.root
output=./test_$RUN.root
$bin $conf $input $output
