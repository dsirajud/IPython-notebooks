Periodic plasma situation

The output files here correspond to doing the regular diagnostics and calculcations, and including all grid points, 0:z.Ngridpoints

the output files in the parent directory only involve points from 0:z.N

For some reason, the data in this folder is not conserving I1 and I2 and a couple others satisfyingly, only to within 10^{-3} or so. Since the periodic boudary is incorporated by simply duplicating the left-end point, there is a difference here between the active grid and the total grid, but it is not clear why there is such a difference in I1 and I2, since the density does not leave the grid, and is only reallocated at each step.
