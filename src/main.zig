const std = @import("std");
const print = std.debug.print;
const io = std.io;
const ac32 = @import("ac.zig");

pub fn main() void {
    const bit_io = ac32.bit_io(
        io.bitReader(.Big, io.getStdIn().reader()),
        io.bitWriter(.Big, io.getStdOut().writer())
    );
    var ac = ac32.arithmeticCoder(bit_io);
    // ac.init(writer);

    print("std-type: {}\n", .{@TypeOf(io.bitWriter)});
    ac.encode(0, 1 << 15);
    ac.encode(1, 1 << 15);
    ac.encode(0, 1 << 15);
    ac.encode(0, 1 << 15);
    // var ac2 = ArithmeticCoder { };
    // ac2.encode(1, 1 << 15);
    // ac2.encode(0, 1 << 15);
}
