const std = @import("std");
const print = std.debug.print;

fn ArithmeticCoder(comptime IOType: type) type {
    return struct {
        io: IOType,
        x: ?u32,
        x1: u32 = 0,
        x2: u32 = @as(u32, std.math.maxInt(u32)),

        var steps: u64 = 0;
        pub fn init(io: IOType) Self { return Self { .io = io, .x = null }; }
        pub fn init_dec(io: IOType) Self { return Self { .io = io, .x = 0 }; }

        const Self = @This();
        pub fn encode(self: *Self, bit: u1, p: u16) void { _ = self.code(bit, p); }
        pub fn decode(self: Self, p: u16) u1 { return self.code(null, p); }

        fn code(self: *Self, maybe_bit: ?u1, prob: u16) ?u1 {
            const p = if (prob == 0) 1 else @as(u64, prob) << 16;
            const xmid = @truncate(u32, self.x1 + ((@as(u64, self.x2 - self.x1) * p) >> 32));

            const bit = maybe_bit orelse @boolToInt(self.x.? <= xmid);
            switch (bit) {
                1 => self.x2 = xmid,
                0 => self.x1 = xmid + 1
            }

            // while ((self.x1 ^ self.x2) >> PREC_SHIFT) == 0 {
            //     self.x1 <<= 1;
            //     self.x2 = (self.x2 << 1) | 1;
            //     self.x = (self.x << 1) | u32::from(self.io.read_bit()?);
            // }

            // while ((self.x1 ^ self.x2) >> PREC_SHIFT) == 0 {
            //     self.io.write_bit(self.x1 >> PREC_SHIFT)?;
            //     self.x1 <<= 1;
            //     self.x2 = (self.x2 << 1) | 1;
            // }
            while ((self.x1 ^ self.x2) >> 31 == 0) {
                if (maybe_bit == null) {
                    self.x.? = (self.x.? << 1) | self.io.readBit();
                } else {
                    self.io.writeBit(self.x1 >> 31);
                }
                self.x1 <<= 1;
                self.x2 = (self.x2 << 1) | 1;
            }

            steps += 1;
            print("steps: {}\n", .{steps});
            if (self.x1 <= xmid and xmid < xmid) {
                return bit;
            } else {
                return null;
            }
        }
    };
}

pub fn arithmeticCoder(stream: anytype) ArithmeticCoder(@TypeOf(stream)) {
    return ArithmeticCoder(@TypeOf(stream)).init(stream);
}

pub fn arithmeticDecoder(stream: anytype) ArithmeticCoder(@TypeOf(stream)) {
    return ArithmeticCoder(@TypeOf(stream)).init_dec(stream);
}

fn BitIO(comptime ReaderType: type, comptime WriterType: type) type {
    return struct {
        reader: ReaderType,
        writer: WriterType,
        const Self = @This();

        pub fn init(reader: ReaderType, writer: WriterType) Self { return Self { .reader = reader, .writer = writer }; }
        pub fn readBit(self: *Self) u1 { return self.reader.readBitsNoEof(u1, 1) catch 0; }
        pub fn writeBit(self: *Self, bit: u32) void { return self.writer.writeBits(bit, 1) catch @panic("write_err"); }
    };
}

pub fn bit_io(reader: anytype, writer: anytype) BitIO(@TypeOf(reader), @TypeOf(writer)) {
    return BitIO(@TypeOf(reader), @TypeOf(writer)).init(reader, writer);
}
