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

    pub fn init() Self { return Self { .ctx = 0, .data = .{Counter.init()}**(1<<12) }; }
    pub fn p(self: Self) u16 { return self.data[self.ctx].p(); }
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
    pub fn p(self: Self) u16 {
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

fn initAC(writer: anytype, comptime mode: Mode) ArithmeticCoder(@TypeOf(writer), mode) {
    return ArithmeticCoder(@TypeOf(writer), mode).init(writer);
}

fn ArithmeticCoder(comptime T: type, comptime mode: Mode) type { return struct {
    io:T, x: if (mode == Mode.d) u32 else void,
    revBits: if (mode == Mode.c) u64 else void,
    x1: u32 = 0, x2: u32 = (1 << 32) - 1,

    const Self = @This();
    const Q1: u32 = 1 << 30; const PREC_SHIFT: u32 = 31;
    const Q2: u32 = 2 << 30; const RLOW_MOD:   u32 = (1 << 31) - 1; // Modify x1 bits in E3 mapping, AND with
    const Q3: u32 = 3 << 30; const RHIGH_MOD:  u32 = (1 << 31) + 1; // Modify x2 bits in E3 mapping, OR with

    pub fn init(io: T) Self { // initialize fields, read state in decode mode
        var self = if (mode == .c) Self { .io = io, .revBits = 0, .x = {} }
              else if (mode == .d) Self { .io = io, .x = 0, .revBits = {} };
        if (mode == .d) self.readState();
        return self;
    }
    pub fn encode(self: *Self, bit: u1, p: u16) !void { return self.proc(bit, p); }
    pub fn decode(self: *Self, p: u16) u1 { return self.proc({}, p); }
    pub fn flush(self: *Self) !void { // flush leading byte to stream
        comptime { assert(mode == .c); }
        try self.writeBit(self.x2 >> PREC_SHIFT);
        while (self.io.bit_count != 0) {
            self.x2 <<= 1; try self.writeBit(self.x2 >> 31);
        }
    }

    fn readBit(self: *Self) u1 { return self.io.readBitsNoEof(u1, 1) catch 0; } // TODO: return 0 only on EOF, otherwise return error
    fn incParity(self: *Self) void { self.revBits += 1; } // for E3 mapping
    fn writeBit(self: *Self, bit: u32) !void { // writes bit, conscious of any E3 mappings
        try self.io.writeBits(bit, 1);
        while (self.revBits > 0) {
            try self.io.writeBits(bit ^ 1, 1);
            self.revBits -= 1;
        }
    }
    fn readState(self: *Self) void { // reads 32-bits into state and pads with zeroes if necessary
        var bitsRead: usize = 0;
        var state = self.io.readBits(u32, 32, &bitsRead) catch 0;
        self.x = state << @intCast(u5, 32 - bitsRead);
    }

    // processes a single bit -> decompresses a bit in decode mode, compresses a bit in encode mode
    fn proc(self: *Self, bit_: if (mode == .d) void else u1, prob: u16) if (mode == .d) u1 else anyerror!void {
        const p = if (prob == 0) 1 else @as(u64, prob) << 16;
        const xmid = @intCast(u32, self.x1 + ((@as(u64, self.x2 - self.x1) * p) >> 32));

        const bit = if (mode == .c) bit_ else @boolToInt(self.x <= xmid);
        if (bit == 1) self.x2 = xmid else self.x1 = xmid + 1;

        while ((self.x1 ^ self.x2) >> PREC_SHIFT == 0) {
            if (mode == .c) try self.writeBit(self.x1 >> PREC_SHIFT)
            else self.x = (self.x << 1) | self.readBit();
            self.x1 <<= 1;
            self.x2 = (self.x2 << 1) | 1;
        }

        while (self.x1 >= Q1 and self.x2 < Q3) {
            if (mode == .c) self.incParity()
            else self.x = ((self.x << 1) ^ Q2) | self.readBit();
            self.x1 = (self.x1 << 1) & RLOW_MOD;
            self.x2 = (self.x2 << 1) | RHIGH_MOD;
        }

        if (mode == .d) return bit;
    }
};}

/// ============================ User Interface =============================
pub fn main() !void {
    var args = std.process.args();
    _ = args.skip(); // skip program invokation
    const mode = parseMode(args.next());
    const inFile = try parseFile(args.next(), FileOptions.read);
    const outFile = try parseFile(args.next(), FileOptions.create);
    defer inFile.close(); defer outFile.close();

    var timer = try std.time.Timer.start();

    var bufr = std.io.bufferedReader(inFile.reader());
    var bufw = std.io.bufferedWriter(outFile.writer());
    var reader = std.io.bitReader(.Big, bufr.reader());
    var writer = std.io.bitWriter(.Big, bufw.writer());
    var model = Model.init();

    if (mode == .c) { // Compression
        const size = try getSize(inFile);
        try writer.writeBits(size, 64);

        var ac = initAC(writer, Mode.c);
        while (true) {
            const bit = reader.readBitsNoEof(u1, 1) catch { break; };
            try ac.encode(bit, model.p());
            model.update(bit);
        }

        try ac.flush();
        try bufw.flush();
    } else { // Decompression
        const size = try reader.readBitsNoEof(u64, 64);
        var i: u64 = 0;

        var ac = initAC(reader, Mode.d);
        while (i / 8 < size) : (i += 1) {
            const bit = ac.decode(model.p());
            try writer.writeBits(bit, 1);
            model.update(bit);
        }

        try bufw.flush();
    }

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
    const s = ms / 1000; if (s < 300) { print("{d:.2} sec\n", .{s}); return; }
    const m = s / 60; if (m < 60) { print("{d:.2} mins\n", .{m}); return; }
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
