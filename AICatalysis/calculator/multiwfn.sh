#!/usr/bin/env bash

until [[ $# -eq 0 ]]; do
  if [[ "$1" =~ \. ]]; then
    file=$1
  elif [[ "$1" == "log" ]]; then
    log=$1
  elif [[ "$1" == "-j" ]]; then
    job=$2
    shift
  fi
  shift
done

if [[ "$file" =~ "wfn" ]]; then
  if [[ "$job" == "ALIE" ]]; then
    echo -e "12\n 2\n 2\n 0\n -1\n -1\n q\n" | /mnt/d/Multiwfn_3.8_dev/Multiwfn $file >$log
  elif [[ "$job" == "ESP" ]]; then
    echo -e "12\n 0\n -1\n -1\n q\n" | /mnt/d/Multiwfn_3.8_dev/Multiwfn $file >$log
  else
    echo -e "12\n 0\n -1\n -1\n q\n" | /mnt/d/Multiwfn_3.8_dev/Multiwfn $file >$log
  fi
elif [[ "$file" =~ "polar" ]]; then
  echo -e "200\n 7\n 1\n 0\n 0\n q\n" | /mnt/d/Multiwfn_3.8_dev/Multiwfn $file >$log
fi
