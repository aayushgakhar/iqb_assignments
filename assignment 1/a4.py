from a1 import global_alignment
from a2 import local_alignment


s = 'ATCAGAGTA'
t = 'TTCAGTA'
match = 2
mismatch = -1

# original result of 1st and 2nd part
gap = -1
print('previous alignment:')
global_alignment(s, t, match, mismatch, gap)
local_alignment(s, t, match, mismatch, gap)

# new result of 1st and 2nd part after making gap = -2
gap = -2
print('\nnew alignment:')
global_alignment(s, t, match, mismatch, gap)
local_alignment(s, t, match, mismatch, gap)

