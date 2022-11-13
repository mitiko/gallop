# Gallop

A lightweight CM compressor for prototyping built with zig, roughly based on [lpaq](http://mattmahoney.net/dc/#lpaq).  
It *gallops*!

Currently: prediction is bitwise, model is order0-ish (12-bit context) with simple 12-bit counters.

## Usage

```bash
./gallop c /data/book1 book1.bin
./gallop c book1.bin book1.orig
cmp book1.orig /data/book1
```

## Goals

- Compress book1 in <2 s, enwik8 in <3 mins
- Under 1KLoC (but readable, lpaq style)
- Under 1GB memory usage
- Optimized for compression ratio

## TODO

- Add a magic number (I propose b"gllp" = 0x676c_6c70)
- Use (12-bit) states instead of counters
- Add state map
- Match model?
- Delayed counters?
- Entropy hashing?

### Build

```bash
zig build-exe gallop.zig
```

### References

Roughly based on [lpaq](http://mattmahoney.net/dc/#lpaq).
Check out my other compressor [weath3rb0i](https://github.com/Mitiko/weath3rb0i).
