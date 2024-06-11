from math import pi

from rs_paths import Line, combine_segments, Arc

l = Line((0, 0), (1, 1))
k = Line((0, 1), (1, 0))
m = Line((1, 0), (0, -1))
n = Arc((0, -1), (0, 0), 1, 1, 0, False, False)

o = Arc.from_center_params((0, 0), 1, 1, 0, 0, pi / 2)
combined = combine_segments([l, k, n, m, o])
print(combined)
