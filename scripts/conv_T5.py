#!/usr/bin/python3
# Copyright Dimitar Rusev 2022

import sys

if len(sys.argv) < 3:
    print("Usage: ./scripts/conv_T5.py fsm.txt out.txt")
    print("Converts FSMs from T5's format (15-bit probabilities of bit = 0) to gallop's 16-bit bit = 1")
    sys.exit(1)

input = sys.argv[1]
output = sys.argv[2]

with open(input, "r") as file:
    with open(output, "w") as out:
        for line in file:
            if line == "\n":
                continue

            data = line.split(',')
            if len(data) != 3:
                print(f"Error parsing FSM for {input}")
                sys.exit(2)

            next0 = int(data[0])
            next1 = int(data[1])
            prob = int(data[2].split('\n')[0])
            # make sure 
            if prob < 0 or prob > (1 << 15):
                print("Make sure this FSM is in T5's format.")
                print("Probabilities are not in range [0, 1 << 15]")
                sys.exit(3)

            if next0 < 0  or next1 < 0 or next0 > (1 << 15) or next1 > (1 << 15):
                print("Make sure this FSM is in T5's format.")
                print("States are not in range [0, 1 << 15]")
                sys.exit(4)

            prob = ((1 << 15) - prob) << 1
            out.write(f"{next0},{next1},{prob}\n")
