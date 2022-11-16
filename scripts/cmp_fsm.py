#!/usr/bin/python3
# Copyright Dimitar Rusev 2022

import sys

if len(sys.argv) < 3:
    print("Usage: ./scripts/cmp_fsm.py fsm_file1.txt fsm_file2.txt [diff|unr]")
    print("Optional parameters: diff, unr")
    print("\tdiff: Print the differences in probabilities")
    print("\tunr: Print unreachable states")
    print("Exits quietly if FSMs are equivalent (like unix's cmp)")
    sys.exit(1)

# Format of FSM files:
# Like T5 from GDCC but the first line may optionally contain a single integer
# to represent initial state if nonzero, otherwise initial state is assumed to be zero
# Each line i (0-index based and excluding the optional inital)
# contains a triplet s0, s1, p and represents state i
#
# s0, s1 are states if next bit is 0 or 1
# p is probability that next bit = 1, scaled by 65536 = 2^16
# For example 32768 = 2^15 represents 1/2, as 32768/65536 = 1/2
# Some probabilities are not expressable this way, but we can find close enough approximations:
# 1/3 ~ 21845/65536, so we can represent it with 21845
# Note this is different than T5's 2^15 scaling and p evaluated for bit=0
# It is however possible to convert T5 statetables to this format using conv_T5.py
#
# --- example1:
# 1, 2, 32768
# 3, 4, 16384
# 4, 5, 49152
# 0, 0, 8192
# 0, 0, 32768
# 0, 0, 1337
# ---
# 
# Notice the current parsing of the FSM allows for single line comments after a third comma:
# 
# --- example2:
# 1
# 1, 1, 32768, This is state 0
# 0, 0, 32768, This is the inital state, as pointed by the first line
# ---
#
# This format allows for unreachable states:
#
# --- example3:
# 1, 2, 32768
# 0, 0, 32768
# 0, 0, 32768
# 9, 9, 64000, This state is unreachable. Pass the `unr` parameter to print unreachable states
# 

# TODO: Seperate unreachable states detection into a seperate script
# TODO: Detect equivalent but non-isomorphic FSMs (hint: need two way mapping) (naive way would be to run DFS on fsm1 then on fsm2)
# TODO: Detect if all probabilities are (1 - p) (predct bit = 0, instead of bit = 1)
# TODO: Detect if all probabilities are p << 1 (if 15-bit probabilities, instead of 16)
# TODO: Add a bunch of asserts that all states, probabilities fit in 16 bits
# TODO: Abstract the parseFSM code so it's reusable by conv_T5
# TODO: Could pass some option to indicate an fsm is of T5's format

# read filenames from args
filename1 = sys.argv[1]
filename2 = sys.argv[2]
opt_print_prob_diffs = False
opt_print_unreachable = False

if len(sys.argv) == 4:
    if sys.argv[3] == "diff":
        opt_print_prob_diffs = True
    elif sys.argv[3] == "unr":
        opt_print_unreachable = True

def parseFSM(filename):
    i = 0
    initial_state = 0
    states = []
    with open(filename, "r") as file:
        for line in file:
            if line == "\n":
                continue

            data = line.split(',')
            if len(data) == 1 and i == 0:
                initial_state = int(data[0])
                continue

            if len(data) < 3:
                print(f"Error parsing FSM for {filename}")
                sys.exit(2)

            i += 1
            next0 = int(data[0])
            next1 = int(data[1])
            prob = int(data[2].split('\n')[0])
            states.append((next0, next1, prob))

    assert(len(states) == i)
    return (initial_state, states)

(init_state1, fsm1) = parseFSM(filename1)
(init_state2, fsm2) = parseFSM(filename2)

# DFS on fsm1 and check there're equivalent nodes in fsm2
# Note it is only possible to check for isomorphism of the two graphs
# because they're deterministic and have a defined start which must match
visited = []
stack = []
stack.append(init_state1)

# Keep track of mapping between fsm1 and fsm2
mapping = { init_state1: init_state2 }
diff_probs = [] # keep track of differences in probabilities

while len(stack) > 0:
    node = stack.pop()
    if node in visited:
        continue
    else:
        visited.append(node)

    (next0_a, next1_a, p_a) = fsm1[node]
    (next0_b, next1_b, p_b) = fsm2[mapping[node]] # node must be in mapping already
    stack.append(next1_a)
    stack.append(next0_a)

    if p_a != p_b:
        diff_probs.append(((node, mapping[node]), (p_a, p_b)))

    # if no mapping has been created yet, we can create it
    if not next0_a in mapping:
        mapping[next0_a] = next0_b

    if not next1_a in mapping:
        mapping[next1_a] = next1_b

    # these do not match
    if mapping[next0_a] != next0_b or mapping[next1_a] != next1_b:
        print("FSMs are not isomorphic.")
        sys.exit(3)

# TODO: Maybe print if FSMs have the exact same mapping, or if they're isomorphic but have a different mapping

# If we've counted some states twice, then they're not truly isomorphic
# but they are equivalent in terms of regular grammars they match
# For now we don't report this equivalence because it's only detected if
# DFS-ing on one of the FSMs - aka order of FSMs matters for this detection (with current alg)
if len(visited) > len(fsm1) or len(visited) > len(fsm2):
    print("FSMs are not isomorphic.")
    exit(4)

if len(diff_probs) > 0:
    print("FSMs are isomorphic but differ in probabilities")
    end = ":\n" if opt_print_prob_diffs else "\n"
    print(f"There are {len(diff_probs)} (out of {len(visited)} total) states with differing probabilities", end=end)
    if opt_print_prob_diffs:
        for ((nodeA, nodeB), (pA, pB)) in diff_probs:
            print(f"({nodeA: <4}, {nodeB: <4}) -- ({pA: <4}, {pB: <4})")

def print_unr(fsm, filename):
    count = len(fsm) - len(visited)
    if count <= 0:
        return
    end = ":\n" if opt_print_unreachable else "\n"
    print(f"There are {count} unreachable states in {filename}", end=end)
    if opt_print_unreachable:
        unreachable_states = []
        for s in range(len(fsm)):
            if not s in visited:
                unreachable_states.append(s)
        for i in range(len(unreachable_states)):
            end = ", "
            if (i + 1) % 16 == 0:
                end = "\n"
            elif i == len(unreachable_states) - 1:
                end = "\n"
            print(f"{unreachable_states[i]: >4}", end=end)
        print()

print_unr(fsm1, filename1)
print_unr(fsm2, filename2)
