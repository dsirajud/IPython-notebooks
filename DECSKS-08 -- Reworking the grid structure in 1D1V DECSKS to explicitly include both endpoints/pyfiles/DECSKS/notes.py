

# everywhere x.N occurs or is affecting the size of vectors

# write x.N be the number of gridpoints in the closed domain x in [a,b]

# then need to define x.width = L / (x.N - 1)

# construct grid to cover x in [a,b], x[0] = a, x[x.N-1] = b

#     for all i in range(Num): # Num instead of self.gridpoints to accomodate for time grid which has t.N = time steps, with t.N + 1 total cells
#         x.gridvalues[i] = a + i*x.width




# need: len(f) = x.N for i = 0, 1, 2, ... , x.N

# x[i] = a + i*x.width

# x.width = L / (x.N - 1)
