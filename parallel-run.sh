#!/bin/bash
# ./parallel-run.sh minK maxK nproc
# intended to be run on cycles, but can be adapted to run elsewhere.

minK="$1"
maxK="$2"
nproc="$3"

for K in $(seq $minK $maxK); do for mle in mle argmax; do for mincov in eps 0.001 0.1; do echo $K $mle $mincov; done; done; done | xargs -L 1 --max-procs=$nproc --verbose /bin/bash -c '/n/fs/cyc424/py/bin/python diag_tf_cv.py $0 $1 $2 data/cv-$0-$1-$2.p'
