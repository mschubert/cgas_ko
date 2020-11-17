#!/bin/sh

cd $(dirname $0)
PROJNAME=${PWD##*/}

rsync -auvr $@ --prune-empty-dirs \
    --exclude "backup/" \
    --exclude "expr_stat1/" \
    --include "*/" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    . ~/Documents/Results/"$PROJNAME"

cd -
