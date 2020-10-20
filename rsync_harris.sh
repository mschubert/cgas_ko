#!/bin/bash

cd $(dirname $0)
PROJNAME=${PWD##*/}

rsync -auvr $@ --include='*.pdf' \
    --exclude '.snakemake' \
    --exclude '.bam' \
    --include '*.RData' \
    --include '*.rds' \
    --include 'data/**.txt' \
    --include 'data/**.tsv' \
    --include 'data/**.tsv.gz' \
    --include 'data/**.mtx' \
    --include 'data/**.mtx.gz' \
    --include='*/' \
    --include='docs/*' \
    --exclude='*' \
    harris:"$PROJNAME"/* .

cd -
