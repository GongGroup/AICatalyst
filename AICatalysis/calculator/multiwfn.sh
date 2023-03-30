#!/usr/bin/env bash

if [[ "$1" =~ "wfn" ]]; then
  echo -e "12\n 0\n -1\n -1\n q\n" | /mnt/d/Multiwfn_3.8_dev/Multiwfn $1 >$2
elif [[ "$1" =~ "polar" ]]; then
  echo -e "200\n 7\n 1\n 0\n 0\n q\n" | /mnt/d/Multiwfn_3.8_dev/Multiwfn $1 >$2
fi
