#!/bin/sh

name=scaler
session=`tmux ls 2>/dev/null | grep $name`
if [ -z "$session" ]; then
    echo "create new session $name"
    tmux new-session -d -s $name \
	"while true; do ./scaler_e72.py; sleep 1; done"
else
    echo "reattach session $name"
    tmux a -t $name
fi
