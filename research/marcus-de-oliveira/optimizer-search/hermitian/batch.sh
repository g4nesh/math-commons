#!/bin/sh
# Usage: sh batch.sh <python> <list file>; each line: tag cls n seed time extra-args...
PY="$1"; LIST="$2"
while read -r tag cls n seed t rest; do
  case "$tag" in ''|\#*) continue;; esac
  echo "== $tag"
  nice -n 10 "$PY" -B hsearch.py --cls "$cls" --n "$n" --seed "$seed" --time "$t" --restarts 100000 --out "runs/$tag.json" $rest
done < "$LIST"
