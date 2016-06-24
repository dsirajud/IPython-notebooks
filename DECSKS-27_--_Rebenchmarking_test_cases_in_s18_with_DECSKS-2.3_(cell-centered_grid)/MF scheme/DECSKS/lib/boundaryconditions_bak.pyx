@cython.boundscheck(False)
def collector_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_k,
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        float dt,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j
    cdef DTYPE_t sigma_nk = 0
    cdef DTYPE_t vzwidth = vz.width

    if k == 0:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] <= 0 and Uf_old[i,j] <= 0:
                    sigma_nk += -vzprepointvaluemesh[i,j] * (f_old[i,j] + Uf_old[i,j]) * charge * vzwidth * dt
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber

    if k == 1:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] <= 0 and Uf_old[i,j] <= 0:
                    sigma_nk += vzprepointvaluemesh[i,j] * Uf_old[i,j] * charge * vzwidth * dt
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber

    # update cumulative charge density
    sim_params['sigma'][z.str]['lower'] += sigma_nk
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def collector_upper_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        float dt,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j
    cdef DTYPE_t sigma_nk = 0
    cdef DTYPE_t vzwidth = vz.width

    if k == 0:
        for i in range(Nz):
            for j in range(Nvz): # rightward velocities (incoming fluxes) at the right boundary
                if zpostpointmesh[i,j] >=  Nz - 1 and Uf_old[i,j] >= 0:
                    sigma_nk += vzprepointvaluemesh[i,j] * (f_old[i,j] + Uf_old[i,j])  * charge * vzwidth * dt
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber

    if k == 1:
        for i in range(Nz):
            for j in range(Nvz): # rightward velocities (incoming fluxes) at the right boundary
                if zpostpointmesh[i,j] >=  Nz - 1 and Uf_old[i,j] >= 0:
                    sigma_nk += vzprepointvaluemesh[i,j] * Uf_old[i,j]  * charge * vzwidth * dt
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber

    # update cumulative charge density sigma = sum_k sum_j sum_n sigma_nk * dvx * dt
    # k = k1, k2, j = 0, ..., Nvx -1, n = 0, 1, 2, ... Nt
    sim_params['sigma'][z.str]['upper'] += sigma_nk
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old
