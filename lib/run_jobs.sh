#!/bin/bash


WDIR=$(pwd)
mkdir $WDIR/$1
cd $WDIR/$1

cp -r $WDIR/input $WDIR/$1
cp -r $WDIR/lib $WDIR/$1

$WDIR/$1/lib/run_pacs.sh 0 5 10
cd $WDIR
