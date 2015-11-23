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

    if HOC['x'] == None:
        pass
    elif HOC['x'] == 'FD' and 'x' in phasespace_vars:
        W['x'] = assemble_finite_difference_weight_matrix(Nx_active,
                                                          N,
                                                          FD_schemes
                                                          )
    elif HOC['x'] == 'FOURIER' and 'x' in phasespace_vars:
        Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                      'x',
                                                      ax, bx,
                                                      Nx, Nx_active,
                                                      Nvx_active,
                                                      N)


    if HOC['y'] == None:
        pass
    elif HOC['y'] == 'FD' and 'y' in phasespace_vars:
        W['y'] = assemble_finite_difference_weight_matrix(Ny_active,
                                                          N,
                                                          FD_schemes
                                                          )
    elif HOC['y'] == 'FOURIER' and 'y' in phasespace_vars:
        Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                      'y',
                                                      ay, by,
                                                      Ny, Ny_active,
                                                      Nvy_active,
                                                      N)

    if HOC['z'] == None:
        pass
    elif HOC['z'] == 'FD' and 'z' in phasespace_vars:
        W['z'] = assemble_finite_difference_weight_matrix(Nz_active,
                                                          N,
                                                          FD_schemes
                                                          )
    elif HOC['z'] == 'FOURIER' and 'z' in phasespace_vars:
        Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                       'z',
                                                       az, bz,
                                                       Nz, Nz_active,
                                                       Nvz_active,
                                                       N)
    if HOC['vx'] == None:
        pass
    elif HOC['vx'] == 'FD' and 'vx' in phasespace_vars:
        W['vx'] = assemble_finite_difference_weight_matrix(Nvx_active,
                                                          N,
                                                          FD_schemes
                                                          )
    elif HOC['vx'] == 'FOURIER' and 'vx' in phasespace_vars:
        Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                       'vx',
                                                       avx, bvx,
                                                       Nvx, Nvx_active,
                                                       Nx_active, # ax will have same number of grid points as x
                                                       N)
    if HOC['vy'] == None:
        pass
    elif HOC['vy'] == 'FD' and 'vy' in phasespace_vars:
        W['vy'] = assemble_finite_difference_weight_matrix(Nvy_active,
                                                          N,
                                                          FD_schemes
                                                          )
    elif HOC['vy'] == 'FOURIER' and 'vy' in phasespace_vars:
        Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                       'vy',
                                                       avy, bvy,
                                                       Nvy, Nvy_active,
                                                       Ny_active, # ay will have same number of grid points as y
                                                       N)
    if HOC['vz'] == None:
        pass
    elif HOC['vz'] == 'FD' and 'vz' in phasespace_vars:
        W['vz'] = assemble_finite_difference_weight_matrix(Nvz_active,
                                                          N,
                                                          FD_schemes
                                                          )
    elif HOC['vz'] == 'FOURIER' and 'vz' in phasespace_vars:
        Xi, xi = assemble_spectral_derivative_operator(Xi, xi,
                                                       'vz',
                                                       avz, bvz,
                                                       Nvz, Nvz_active,
                                                       Nz_active, # az will have same number of grid points as z
                                                       N)


