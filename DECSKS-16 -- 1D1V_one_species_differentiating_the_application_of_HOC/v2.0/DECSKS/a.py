    # could just do

    for var in phasespace_vars:
        if HOC[var] == 'FD':
            W[var] = assemble_finite_difference_weight_matrix(
                eval('N' + var + '_active'),
                N,
                FD_schemes
                )
        elif HOC[var] == 'FOURIER':
            # ensure the correct number of grid points
            # is passed for the generalized velocity Nvz_active
            # for x,y,z, 'vz' = vx, vy, vz
            # for vx, vy, vz, 'vz' = ax, ay, az, which have
            # the same number of dims as x, y, z, respectively

            if var[0] == 'v':
                Nvz_active = eval('N' + var[1] + '_active')
            else:
                Nvz_active = eval('Nv' + var + '_active')

            Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                          var,
                                                          eval('a' + var),
                                                          eval('b' + var),
                                                          eval('N' + var),
                                                          eval('N' + var + '_active'),
                                                          Nvz_active,
                                                          N)

    print Xi['x'].shape
    print Xi['vx'].shape
    print xi.viewkeys()
