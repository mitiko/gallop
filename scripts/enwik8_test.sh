#!/usr/bin/sh

./gallop c /data/enwik8 enwik8.bin &&
./gallop d enwik8.bin enwik8.orig && 
cmp enwik8.orig /data/enwik8
