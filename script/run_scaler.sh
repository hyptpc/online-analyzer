#!/bin/sh

. $(dirname $(readlink -f $0))/setebhost
program=scaler

#______________________________________________________________________________
top_dir=$(dirname $(readlink -f $0))/..
server=$top_dir/bin/$program

conf=/misc/software/param/conf/analyzer_e72_20250506.conf
data=$ebhost:8901

#______________________________________________________________________________
name=scaler
session=`tmux ls 2>/dev/null | grep $name`
if [ -z "$session" ]; then
    echo "create new session $name"
    tmux new-session -d -s $name \
	"while true; do $server $conf $data; sleep 1; done"
else
    echo "reattach session $name"
    tmux a -t $name
fi
