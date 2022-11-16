#!/usr/bin/sh

./gallop c /data/book1 book1.bin &&
./gallop d book1.bin book1.orig && 
cmp book1.orig /data/book1
