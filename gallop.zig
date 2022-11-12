const std = @import("std");
const print = std.debug.print;

const Mode = enum { c, d };

pub fn main() !void {
    var args = std.process.args();
    _ = args.skip(); const arg = args.next();
    if (arg == null) { print("Error: No args passed. Pass -c for compression , -d for decompression\n", .{}); std.os.exit(1); }
    
    const mode = if (std.mem.eql(u8, arg.?, "-d")) Mode.d
            else if (std.mem.eql(u8, arg.?, "-c")) Mode.c
            else null;
    if (mode == null) { print("Error: Invalid arg. Pass -c for compression , -d for decompression\n", .{}); std.os.exit(2); }

    var bufw = std.io.bufferedWriter(std.io.getStdOut().writer());
    var writer = std.io.bitWriter(.Big, bufw.writer());
    var model = Model.init();

    if (mode.? == .c) {
        const fileName = args.next();
        if (fileName == null) { print("Error: To compress use -c fileName\n", .{}); std.os.exit(3); }
        var path_buffer: [std.fs.MAX_PATH_BYTES]u8 = undefined;
        const path = try std.fs.realpathZ(fileName.?, &path_buffer);
        const file = try std.fs.openFileAbsolute(path, .{});
        defer file.close();
        const size = (try file.stat()).size;
        var bufr = std.io.bufferedReader(file.reader());
        var reader = std.io.bitReader(.Big, bufr.reader());

        try writer.writeBits(size, 64);
        var ac = initAC(writer, Mode.c);
        while (true) {
            const bit = reader.readBitsNoEof(u1, 1) catch { break; };
            try ac.encode(bit, model.p());
            model.update(bit);
        }
        try ac.flush();
        try bufw.flush();
    } else {
        var bufr = std.io.bufferedReader(std.io.getStdIn().reader());
        var reader = std.io.bitReader(.Big, bufr.reader());
        const size = try reader.readBitsNoEof(u64, 64);
        var ac = initAC(reader, Mode.d);

        var i: u64 = 0;
        while (i / 8 < size) : (i += 1) {
            const bit = ac.decode(model.p());
            try writer.writeBits(bit, 1);
            model.update(bit);
        }
        try bufw.flush();
    }
}

const Model = struct {
    ctx: u12,
    data: [1 << 12]Counter,

    const Self = @This();

    pub fn init() Self { return Self { .ctx = 0, .data = .{Counter.init()}**(1<<12) }; }
    pub fn p(self: Self) u16 { return self.data[self.ctx].p(); }
    pub fn update(self: *Self, bit: u1) void {
        self.data[self.ctx].update(bit);
        self.ctx <<= 1; self.ctx |= bit; self.ctx &= (1 << 12) - 1;
    }
};

const Counter = struct {
    // c0: u16, c1: u16,
    c0: u12, c1: u12,

    const Self = @This();
    pub fn init() Self { return Self { .c0 = 0, .c1 = 0 }; }
    pub fn p(self: Self) u16 {
        const n0 = @as(u64, self.c0);
        const n1 = @as(u64, self.c1);
        return @intCast(u16, (1 << 16) * (n1 + 1) / (n1 + n0 + 2));
    }
    pub fn update(self: *Self, bit: u1) void {
        // const maxCount = (1 << 16) - 1;
        const maxCount = (1 << 12) - 1;
        if (self.c0 == maxCount or self.c1 == maxCount) {
            self.c0 >>= 1;
            self.c1 >>= 1;
        }
        if (bit == 1) self.c1 += 1 else self.c0 += 1;
    }
};

fn initAC(writer: anytype, comptime mode: Mode) ArithmeticCoder(@TypeOf(writer), mode) { return ArithmeticCoder(@TypeOf(writer), mode).init(writer); }
fn ArithmeticCoder(comptime T: type, comptime mode: Mode) type { return struct {
    io:T, x: if (mode == Mode.d) u32 else void,
    revBits: if (mode == Mode.c) u64 else void,
    x1: u32 = 0, x2: u32 = (1 << 32) - 1,

    const Self = @This();
    pub fn init(io: T) Self {
        var self = if (mode == .c) Self { .io = io, .revBits = 0, .x = {} }
              else if (mode == .d) Self { .io = io, .x = 0, .revBits = {} };
        if (mode == .d) self.readState();
        return self;
    }
    pub fn encode(self: *Self, bit: u1, p: u16) !void { return self.code(bit, p); }
    pub fn decode(self: *Self, p: u16) u1 { return self.code({}, p); }
    pub fn flush(self: *Self) !void {
        try self.writeBit(self.x2 >> 31);
        while (self.io.bit_count != 0) {
            self.x2 <<= 1; try self.writeBit(self.x2 >> 31);
        }
    }

    fn readBit(self: *Self) u1 { return self.io.readBitsNoEof(u1, 1) catch 0; }
    fn incParity(self: *Self) void { self.revBits += 1; }
    fn writeBit(self: *Self, bit: u32) !void {
        try self.io.writeBits(bit, 1);
        while (self.revBits > 0) {
            try self.io.writeBits(bit ^ 1, 1);
            self.revBits -= 1;
        }
    }
    fn readState(self: *Self) void {
        var bitsRead: usize = 0;
        var state = self.io.readBits(u32, 32, &bitsRead) catch 0;
        self.x = state << @intCast(u5, 32 - bitsRead);
    }

    fn code(self: *Self, bit_: if (mode == .d) void else u1, prob: u16) if (mode == .d) u1 else anyerror!void {
        const p = if (prob == 0) 1 else @as(u64, prob) << 16;
        const xmid = @intCast(u32, self.x1 + ((@as(u64, self.x2 - self.x1) * p) >> 32));

        const bit = if (mode == .c) bit_ else @boolToInt(self.x <= xmid);
        if (bit == 1) self.x2 = xmid else self.x1 = xmid + 1;

        while ((self.x1 ^ self.x2) >> 31 == 0) {
            if (mode == .c) try self.writeBit(self.x1 >> 31)
            else self.x = (self.x << 1) | self.readBit();
            self.x1 <<= 1;
            self.x2 = (self.x2 << 1) | 1;
        }

        while (self.x1 >= (1 << 30) and self.x2 < (3 << 30)) {
            if (mode == .c) self.incParity()
            else self.x = ((self.x << 1) ^ (2 << 30)) | self.readBit();
            self.x1 = (self.x1 << 1) & ((1 << 31) - 1);
            self.x2 = (self.x2 << 1) | ((1 << 31) + 1);
        }

        if (mode == .d) return bit;
    }
};}
