#!/bin/bash
# xargs-compute-gd-mle.sh infile N nproc prefix
# Runs in parallel xargs-compute-gd-mle.sh for the posterior pickled
# file infile assuming that the number of posterior estimates is N.
#
# See compute-gd-mle.py documentation for details.
# Should only be run on cycles, but can be adapted to run elsewhere.
# Should only be run in the github root directory.

infile="$1"
N="$2"
nproc="$3"
prefix="$4"

seq 0 $(( $N / $nproc )) $N | \
  sed 'p' | head -n -1 | tail -n +2 | \
  sed '$!N;s/\n/ /' | \
  xargs --max-procs=$nproc --replace --verbose \
    /bin/bash -c "/n/fs/cyc424/py/bin/python3 compute-gd-mle.py $infile data/$prefix-mle-\$(echo {} | tr ' ' '-').p {}"

