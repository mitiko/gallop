const std = @import("std");
const print = std.debug.print;
const assert = std.debug.assert;
const File = std.fs.File;

/// ============================= Model =============================
/// A simple order0-ish model (with 12-bit context)
/// To be replaced with mixer + micromodels
/// MicroModels will share 2 hashtables, and hopefully 2 statetables - for big contexts, and for small ones
/// Mixer should use vectors (since zig is chill like that)
const Model = struct {
    ctx: u12,
    data: [1 << 12]Counter,

    const Self = @This();

    pub fn init() Self { return Self { .ctx = 0, .data = .{Counter.init()} ** (1 << 12) }; }
    pub fn p(self: *Self) u16 { return self.data[self.ctx].p(); }
    pub fn update(self: *Self, bit: u1) void {
        self.data[self.ctx].update(bit);
        self.ctx <<= 1; self.ctx |= bit;
    }
};

/// ============================= Counter =============================
/// A simple u12 counter (takes 3-bytes of memory)
/// To be replaced with a state table + state map
const Counter = struct {
    counts: [2]u12,

    const Self = @This();
    pub fn init() Self { return Self { .counts = [_]u12 { 0, 0 } }; }
    pub fn p(self: *Self) u16 {
        const n0 = @as(u64, self.counts[0]);
        const n1 = @as(u64, self.counts[1]);
        return @intCast(u16, (1 << 16) * (n1 + 1) / (n1 + n0 + 2));
    }
    pub fn update(self: *Self, bit: u1) void {
        const maxCount = (1 << 12) - 1;
        if (self.counts[0] == maxCount or self.counts[1] == maxCount) {
            self.counts[0] >>= 1;
            self.counts[1] >>= 1;
        }
        self.counts[bit] += 1;
    }
};

/// ============================= Arithmetic coder =============================
/// 32-bit (binary) arithmetic coder
/// Use `initAC(writer, Mode.c)` for encoding, and `initAC(reader, Mode.d)` for decoding
/// Initializing in wrong mode wouldn't compile because of the way zig emulates generics
/// To encode: `try ac.encode(bit, p1)`, To decode: `const bit = ac.decode(p1)`
/// Expected io is `std.io.BitReader` or `std.io.BitWriter`
/// `flush()` should be called exactly once
const Mode = enum { c, d }; // (compression, decompression) = (encode, decode)

fn initAC(io: anytype, comptime mode: Mode) !ArithmeticCoder(@TypeOf(io), mode) {
    return ArithmeticCoder(@TypeOf(io), mode).init(io);
}

fn ArithmeticCoder(comptime T: type, comptime mode: Mode) type { return struct {
    io:T, x: if (mode == Mode.d) u32 else void,
    revBits: if (mode == Mode.c) u64 else void,
    x1: u32 = 0, x2: u32 = (1 << 32) - 1,

    const Self = @This();
    const Q1: u32 = 1 << 30; const PREC_SHIFT: u32 = 31;
    const Q2: u32 = 2 << 30; const RLOW_MOD:   u32 = (1 << 31) - 1; // Modify x1 bits in E3 mapping, AND with
    const Q3: u32 = 3 << 30; const RHIGH_MOD:  u32 = (1 << 31) + 1; // Modify x2 bits in E3 mapping, OR with

    pub fn init(io: T) !Self { // initialize fields, read state in decode mode
        var self = if (mode == .c) Self { .io = io, .revBits = 0, .x = {} }
              else if (mode == .d) Self { .io = io, .x = 0, .revBits = {} };
        if (mode == .d) self.x = try self.io.read_u32();
        return self;
    }
    pub fn encode(self: *Self, bit: u1, p: u16) !void { return self.proc(bit, p); }
    pub fn decode(self: *Self, p: u16) !u1 { return self.proc({}, p); }
    pub fn flush(self: *Self) !void { // flush leading byte to stream
        comptime { assert(mode == .c); }
        try self.writeBit(self.x2 >> PREC_SHIFT);
        while (!self.io.bitQueue.isEmpty()) : (self.x2 <<= 1) {
            try self.writeBit(self.x2 >> PREC_SHIFT);
        }
    }

    fn incParity(self: *Self) void { self.revBits += 1; } // for E3 mapping
    fn writeBit(self: *Self, bit: u32) !void { // writes bit, conscious of any E3 mappings
        try self.io.writeBit(@intCast(u1, bit));
        while (self.revBits > 0) : (self.revBits -= 1) {
            try self.io.writeBit(@intCast(u1, bit ^ 1));
        }
    }

    // processes a single bit -> decompresses a bit in decode mode, compresses a bit in encode mode
    const returnType: type = anyerror!(if (mode == .d) u1 else void);
    fn proc(self: *Self, bit_: if (mode == .d) void else u1, prob: u16) returnType {
        const p = if (prob == 0) 1 else @as(u64, prob) << 16;
        const xmid = @intCast(u32, self.x1 + ((@as(u64, self.x2 - self.x1) * p) >> 32));

        const bit = if (mode == .c) bit_ else @boolToInt(self.x <= xmid);
        if (bit == 1) self.x2 = xmid else self.x1 = xmid + 1;

        while ((self.x1 ^ self.x2) >> PREC_SHIFT == 0) {
            if (mode == .c) try self.writeBit(self.x1 >> PREC_SHIFT)
            else self.x = (self.x << 1) | try self.io.readBit();
            self.x1 <<= 1;
            self.x2 = (self.x2 << 1) | 1;
        }

        while (self.x1 >= Q1 and self.x2 < Q3) {
            if (mode == .c) self.incParity()
            else self.x = ((self.x << 1) ^ Q2) | try self.io.readBit();
            self.x1 = (self.x1 << 1) & RLOW_MOD;
            self.x2 = (self.x2 << 1) | RHIGH_MOD;
        }

        if (mode == .d) return bit;
    }
};}

/// ============================== IO Types ===============================
/// Buffered Bit IO types with 8KB buffers tailored to AC's specific needs
/// 1. Reader returns 0 on EOF
/// 2. Read/Write methods for 32-bit and 64-bit unsigned integers
const BitBufReader = struct {
    inner: File.Reader, bitQueue: BitQueue,
    buf: [1 << 13]u8 = undefined, // 8KiB
    amt: usize = 0, idx: usize = 0, eof: bool = false,

    const Self = @This();
    pub fn init(inner: File.Reader) Self { return Self { .inner = inner, .bitQueue = BitQueue.init() }; }

    pub fn readBit(self: *Self) !u1 {
        const bit = self.bitQueue.pop();
        if (bit != null) return bit.?;

        const byte = try self.readByte();
        self.bitQueue.fill(byte);
        return self.bitQueue.pop().?;
    }

    pub fn readByte(self: *Self) !u8 {
        if (self.eof) return 0;
        if (self.idx == self.amt) {
            self.idx = 0;
            self.amt = try self.inner.read(&self.buf);
            if (self.amt == 0) { self.eof = true; return 0; }
        }

        const byte = self.buf[self.idx];
        self.idx += 1;
        return byte;
    }

    pub fn read_u32(self: *Self) !u32 {
        assert(self.bitQueue.isEmpty());
        var res: u32 = 0; comptime var i = 0;
        inline while (i < 4) : (i += 1) {
            res = (res << 8) | try self.readByte();
        }
        return res;
    }

    pub fn read_u64(self: *Self) !u64 {
        assert(self.bitQueue.isEmpty());
        var res: u64 = 0; comptime var i = 0;
        inline while (i < 8) : (i += 1) {
            res = (res << 8) | try self.readByte();
        }
        return res;
    }
};

const BitBufWriter = struct {
    inner: File.Writer, bitQueue: BitQueue,
    buf: [1 << 13]u8 = undefined, // 8KiB
    idx: usize = 0,

    const Self = @This();
    pub fn init(inner: File.Writer) Self { return Self { .inner = inner, .bitQueue = BitQueue.init() }; }

    pub fn writeBit(self: *Self, bit: u1) !void {
        self.bitQueue.push(bit);
        const byte = self.bitQueue.flush();
        if (byte == null) return;

        if (self.idx == self.buf.len) {
            try self.inner.writeAll(&self.buf);
            self.idx = 0;
        }

        self.buf[self.idx] = byte.?;
        self.idx += 1;
    }

    pub fn write_u64(self: *Self, valc: u64) !void {
        assert(self.bitQueue.isEmpty());
        assert(self.idx == 0);

        var val = valc; comptime var i = 0;
        inline while (i < 8) : (i += 1) {
            self.buf[self.idx] = @intCast(u8, val >> 56);
            self.idx += 1;
            val <<= 8;
        }
    }

    pub fn flush(self: *Self) !void {
        assert(self.bitQueue.isEmpty());
        try self.inner.writeAll(self.buf[0..self.idx]);
    }
};

/// ============================= Bit Queue ===============================
/// Handles bit io
const BitQueue = struct {
    t: u8, cnt: u4,

    const Self = @This();
    pub fn init() Self { return Self { .t = 0, .cnt = 0 }; }

    pub fn fill(self: *Self, byte: u8) void {
        assert(self.isEmpty());
        self.t = byte; self.cnt = 8;
    }

    pub fn flush(self: *Self) ?u8 {
        if (self.cnt != 8) return null;
        self.cnt = 0;
        return self.t;
    }

    pub fn pop(self: *Self) ?u1 {
        if (self.isEmpty()) return null;
        self.cnt -= 1;
        return @intCast(u1, (self.t >> @intCast(u3, self.cnt)) & 1);
    }

    pub fn push(self: *Self, bit: u1) void {
        assert(self.cnt < 8);
        self.t = (self.t << 1) | bit;
        self.cnt += 1;
    }

    pub inline fn isEmpty(self: Self) bool { return self.cnt == 0; }
};

/// ============================ User Interface =============================
pub fn main() !void {
    var args = std.process.args();
    _ = args.skip(); // skip program invokation
    const mode = parseMode(args.next());
    const inFile = try parseFile(args.next(), FileOptions.read);
    const outFile = try parseFile(args.next(), FileOptions.create);
    defer inFile.close(); defer outFile.close();

    var timer = try std.time.Timer.start();

    var reader = BitBufReader.init(inFile.reader());
    var writer = BitBufWriter.init(outFile.writer());
    var model = Model.init();

    const size = switch(mode) {
        .c => try getSize(inFile),
        .d => try reader.read_u64()
    };
    if (mode == .c) try writer.write_u64(size);

    var ac = switch(mode) {
        .c => try initAC(&writer, Mode.c),
        .d => try initAC(&reader, Mode.d)
    };

    var i: u64 = 0;
    while (i < size) : (i += 1) {
        comptime var k = 0;
        inline while (k < 8) : (k += 1) {
            const bit = switch(mode) {
                .c => try reader.readBit(),
                .d => try ac.decode(model.p())
            };
            switch(mode) {
                .c => try ac.encode(bit, model.p()),
                .d => try writer.writeBit(bit)
            }
            model.update(bit);
        }
    }
    try writer.flush();

    const ns = @intToFloat(f64, timer.lap());
    const inSize = try getSize(inFile);
    const outSize = try getSize(outFile);
    reportResult(mode, inSize, outSize, ns);
}

fn parseMode(arg: ?[:0]const u8) Mode {
    if (arg == null) exit(1);
    const mode = if (std.mem.eql(u8, arg.?, "c")) Mode.c
            else if (std.mem.eql(u8, arg.?, "d")) Mode.d
            else null;
    if (mode == null) exit(2);
    return mode.?;
}

const FileOptions = enum { create, read };
fn parseFile(arg: ?[:0]const u8, options: FileOptions) !File {
    if (arg == null) exit(3);
    if (options == .create) return std.fs.cwd().createFileZ(arg.?, .{});
    var pathBuf: [std.fs.MAX_PATH_BYTES]u8 = undefined;
    const path = try std.fs.realpathZ(arg.?, &pathBuf);
    return std.fs.openFileAbsolute(path, .{});
}

fn getSize(f: File) !u64 { return (try f.stat()).size; }

fn reportResult(mode: Mode, inSize: u64, outSize: u64, ns: f64) void {
    switch (mode) {
        .c => print("Compressed   {} -> {} in ", .{inSize, outSize}),
        .d => print("Decompressed {} -> {} in ", .{inSize, outSize})
    }

    if (ns < 1000) { print("{d:.0} ns\n", .{ns}); return; }
    const us = ns / 1000; if (us < 1000) { print("{d:.3} us\n", .{us}); return; }
    const ms = us / 1000; if (ms < 1000) { print("{d:.2} ms\n", .{ms}); return; }
    const s  = ms / 1000; if (s  < 300)  { print("{d:.2} sec\n", .{s}); return; }
    const m  = s / 60;    if (m  < 60)   { print("{d:.2} mins\n",.{m}); return; }
    const h = m / 60; print("{d:.2} hr\n", .{h});
}

fn exit(status: u8) void {
    print(
        \\gallop file compressor (C) 2022, Dimitar Rusev (mitiko)
        \\
        \\To compress:   ./gallop c input output
        \\To decompress: ./gallop d input output
        \\Example: (./gallop c /data/book1 book1.bin) && (./gallop d book1.bin book1.orig) && (cmp book1.orig /data/book1)
        \\
    ,.{});
    std.os.exit(status);
}
