import DECSKS
import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t

# PYTHON METHODS
def periodic(f_k,
             f_old,
             Uf,
             t,
             z,
             vz,
             dt,
             sim_params,
             charge,
             k = 0
             ):
    """Applies periodic boundary conditions to
    postpointmesh

    inputs:
    f_old -- (ndarray, ndim=2) density array
    z -- (instance) phase space variable being evolved

        z.postpointmesh -- (ndarray, ndim=3),
                           shape = (2, x.N, vx.N)

    outputs:
    f_old -- (ndarray, ndim=2) Array with both periodic
             BCs being enforce
    z    -- (instance) phase sapce variable being evolved with
             updated attribute z.postpointmesh

    f_old, Uf returned for symmetry with nonperiodic routine below
    """
    z.postpointmesh[k,:,:] = np.mod(z.postpointmesh[k,:,:], z.N)
    vz.postpointmesh[k,:,:] = vz.prepointmesh.copy()     # assign to a copy so that changes to postpointmesh do not affect prepointmesh

    return f_k, f_old, Uf, z

def nonperiodic(f_k,
                f_old,
                Uf,
                t,
                z,
                vz,
                dt,
                sim_params,
                charge,
                k = 0
                ):
    """orchestrates applying nonperiodic boundary conditions
    to the array w with total active grid points Nw. Nonperiodic
    boundary conditions require handling both left and right
    boundaries

    inputs:
    f_old -- (ndarray, ndim=2) density array
    z -- (instance) phase space variable being evolved

    outputs:
    f_old -- (ndarray, ndim=2) density with both left and right
             nonperiodic BCs enforced
    Uf -- (ndarray, ndim=2) high order fluxes with both left and right
             nonperiodic BCs enforced


    z returned (no changes) for symmetry with periodic routine above
    """

    # initialize vz.postpointmesh to be prepointmesh values. Some boundary routines may require modifying these, most will not
    vz.postpointmesh[k,:,:] = vz.prepointmesh.copy() # assign to a copy so that changes to postpointmesh do not affect prepointmesh

    # lower boundary
    f_k, f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(
                                              f_k,     # only lower symmetry boundary would need access to this
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              dt,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)

    # upper boundary
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              dt,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)


    #    vz.postpointmesh[k,:,:] = vz.prepointmesh
    return f_k, f_old, Uf, z

# CYTHON METHODS
@cython.boundscheck(False)
def absorbing_lower_boundary(np.ndarray[DTYPE_t, ndim=2] f_k,
                              np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    if k == 0:
        #        print "k = 0 output LABC:"
        for j in range(Nvz):
            for i in range(Nz):
                if vzprepointvaluemesh[i,j] < 0:
                    # there is  a case where z.CFL.frac < 0, so k < 0 for near enough i, but if computed Uf > 0, the lim[Uf] = 0
                    # so we are not counting the cases with negative postpoints suddenly because we are skipping over the limited Uf = 0 cases
                    # thus we need to actually check sign on z.CFL.numbers = sign on vz.prepointvaluemesh[i,j] * dt = sign on vz.prepointvaluemesh
                    # if we are correct and only use a boundary routine with a split scheme with positive time steps since V-P system is not
                    # reversible

                    # note that the case vzprepointmesh[i,j] = 0 (z.CFL.numbers[i,j] = 0), in an edge case produces a
                    # postpoint z.postpointmesh[1,i,j] = 0 + np.sign(z.CFL.numbers[stage,i,j]) = 0 + 0 = 0
                    # so that the remap to z.postpointmesh[1,i,j] assigns Uf = 0 to i = 0
                    # wheras the remap to z.postpointmesh[0,i,j] assigns f + Uf = f to i = 0.
                    # Since the remap routine is summative, there is no overwriting of information and this proceeds consistently.

                    if zpostpointmesh[i,j] == 0 and -Uf_old[i,j] <= 1/2. * f_old[i,j]: # then $k_{true}\in [-1/2,0]$
                        Uf_old[i,j] = 2*Uf_old[i,j]
                    elif zpostpointmesh[i,j] <= 0: # then $k_{true}\in (-\infty , -1/2)$
                        Uf_old[i,j] = 0
                        f_old[i,j] = 0
                        zpostpointmesh[i,j] = 0 # dummy index to avoid off-grid index reference error


    if k == 1:
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] <= -1:
                        # then $k_{true}\in (-infty , 0]$, the k2 proportion is always allocated to i <= -1 (beyond the wall)
                        # in an ABC, the proportion being absorbed is unimportant as this proportion is zeroed out
                        # see notebook s24 for derivation of the ABC prescription
                        #                        print "chi_2 f[x=%g,vx=%g] -> 0" % (z.prepointvaluemesh[i,j], vzprepointvaluemesh[i,j])

                        Uf_old[i,j] = 0
                        zpostpointmesh[i,j] = 0 # dummy index to avoid off-grid index reference error

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def absorbing_upper_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    if k == 0:
        for j in range(Nvz):
            for i in range(Nz):
                if vzprepointvaluemesh[i,j] > 0: # the case for Uf = 0 indicates any such prepoint density remains where it is, zero is absorbed at a wall nearby
                    if zpostpointmesh[i,j] == Nz-1 and Uf_old[i,j] <= 1/2. * f_old[i,j]: # then $k_{true}\in [Nz-1, Nz-1/2]$
                        Uf_old[i,j] = 2*Uf_old[i,j]

                    elif zpostpointmesh[i,j] >= Nz-1: # then $k_{true}\in (Nz-1/2, \infty )$
                        Uf_old[i,j] = 0
                        f_old[i,j] = 0
                        zpostpointmesh[i,j] = Nz-1 # dummy index to avoid off-grid index reference error

    if k == 1:
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] > 0:
                    if zpostpointmesh[i,j] >=  Nz:
                        # then $k_{true}\in (Nz-1, \infty )$, the k2 proportion is always allocated to i >= Nz (beyond the wall)
                        # in an ABC, the proportion being absorbed is unimportant as this proportion is zeroed out
                        # see notebook s24 for derivation of the ABC prescription

                        Uf_old[i,j] = 0
                        zpostpointmesh[i,j] = Nz-1

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old


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
    cdef DTYPE_t sigma_nks = 0    # n = time step, k = 0 (nearest gridpoint) or 1 (contiguous gridpoint), s = species with charge q_s = int charge
                                  # total charge sigma is sum over k, integration over n, and sum over s
                                  # $\sigma_{nks} = q_s N_{ks}
                                  #
                                  # where N_{ks} is the total number allocated from proportion k of species s that reach the wall, i.e.
                                  #
                                  #        N_{ks}  \simeq \sum_{i,j} \Delta x\Delta v_x chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
    cdef DTYPE_t vzwidth = vz.width
    cdef DTYPE_t  zwidth =  z.width

    if k == 0: # often referred to as "k1"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] == 0 and -Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [-1/2, 0]$

                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = -1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = 0 (on-grid)

                    elif zpostpointmesh[i,j] <= 0: # $k_{true}\in (-\infty , -1/2)$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k1, $\chi_{k_1} f_{i,j}^n$ where $\chi_{k_1} = (1 + U_{i,j}^n)$, $U_{i,j}^n <= 0$

                        # charge collection:
                        #
                        sigma_nks += (f_old[i,j] + Uf_old[i,j])

                        # remove this same proportion ("f + Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        f_old[i,j] = 0
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = 0

    if k == 1: # often referred to as "k2"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] == -1 and -Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [-1/2, 0]$ (k1 = 0, k2 = -1/2)

                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = -1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = -1/2 -> -1 (wall)

                        #collect reproportioned amount  allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = -U_{i,j}^n$, $U_{i,j}^n <= 0$
                        # charge collection:
                        sigma_nks += -Uf_old[i,j] # Uf < 0 since vzprepointvaluemesh[i,j] < 0 and split scheme must have dt > 0 for bdry cases

                        # remove this same reproportioned amount ("-Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = 0

                    elif zpostpointmesh[i,j] <= -1: # $k_{true}\in (-\infty , -1/2)$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = -U_{i,j}^n$, $U_{i,j}^n <= 0$

                        # charge collection:
                        #
                        #    $$\sigma_{m\ell}^\alpha} \equiv
                        sigma_nks += -Uf_old[i,j] # Uf < 0 since vzprepointvaluemesh[i,j] < 0 and split scheme requires dt > 0 for bdry cases
                        # note: minus signs have been cancelled on on (-vzprepointvaluemesh[i,j])*(-Uf_old[i,j])

                        # remove this same proportion ("-Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        Uf_old[i,j] = 0
                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = 0

    # at this point sigma_nks = sum (chi * f[i,j]) over i,j, need sigma_nks = charge * sum( vzwidth * zwidth * chi * f[i,j])
    # where chi = 1 + U or -U
    # include multiplicative factors from integration (vzwidth, zwidth) and charge
    sigma_nks *= charge
    sigma_nks *= vzwidth # for integration in velocity space
    sigma_nks *= zwidth # for integration in configuration space

    # update cumulative charge density, this sum below accomplishes the following required sums and integrals:
    #
    #    o  integrates over time with width dt = split_coeff*t.width
    #    o  sums over k after k1 and k2 both go through this routine
    #    o  sums over species s when fe and fi both go through this routine
    sim_params['sigma'][z.str]['lower'] += sigma_nks
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
    cdef DTYPE_t sigma_nks = 0    # n = time step, k = 0 (nearest gridpoint) or 1 (contiguous gridpoint), s = species with charge q_s = int charge
                                  # total charge sigma is sum over k, integration over n, and sum over s
                                  # $\sigma_{nks} = q_s \int dv_x |v_{x,j}| chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
                                  #          \simeq q_s\sum_{i,j} \Delta v_x |v_{x,j}| chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
    cdef DTYPE_t vzwidth = vz.width
    cdef DTYPE_t  zwidth =  z.width

    if k == 0: # often referred to as "k1"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] > 0:
                    if zpostpointmesh[i,j] == Nz-1 and Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [Nz-1,Nz-1/2]$

                        # maps on-grid k1 = Nz-1
                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = Nz-1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = Nz-1 (on-grid)

                    elif zpostpointmesh[i,j] >= Nz-1: # $k_{true}\in (Nz-1/2,\infty )$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k1, $\chi_{k_1} f_{i,j}^n$ where $\chi_{k_1} = (1 - U_{i,j}^n)$, $U_{i,j}^n >= 0$

                        # charge collection:
                        #
                        sigma_nks += (f_old[i,j] - Uf_old[i,j])

                        # remove this same proportion ("f + Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        f_old[i,j] = 0
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = Nz-1

    if k == 1: # often referred to as "k2"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] == Nz and Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [-1/2, 0]$ (k1 = 0, k2 = -1/2)

                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = Nz-1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = Nz-1/2 -> Nz (wall)

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = -U_{i,j}^n$, $U_{i,j}^n <= 0$
                        # charge collection:
                        sigma_nks += Uf_old[i,j]

                        # remove this same proportion ("Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = Nz-1

                    elif zpostpointmesh[i,j] >= Nz: # $k_{true}\in (Nz-1/2, \infty )$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = U_{i,j}^n$, $U_{i,j}^n <= 0$

                        # charge collection:
                        #
                        #    $$\sigma_{m\ell}^\alpha} \equiv
                        sigma_nks += Uf_old[i,j]
                        # note: minus signs have been cancelled on on (-vzprepointvaluemesh[i,j])*(-Uf_old[i,j])

                        # remove this same proportion ("-Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        Uf_old[i,j] = 0
                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = Nz-1

    # at this point sigma_nks = sum (chi[k,i,j] * f[i,j]) over i,j, need sigma_nks = charge * sum( vzwidth * zwidth * chi * f[i,j])
    # where chi[0,i,j] = 1 + U[i,j] and chi[1,i,j] = -U
    # include multiplicative factors from integration (vzwidth, zwidth) and charge
    sigma_nks *= charge
    sigma_nks *= vzwidth # for integration in velocity space
    sigma_nks *= zwidth # for integration in configuration space

    # update cumulative charge density, this sum below accomplishes the following required sums and integrals:
    #
    #    o  integrates over time with width dt = split_coeff*t.width
    #    o  sums over k after k1 and k2 both go through this routine
    #    o  sums over species s when fe and fi both go through this routine
    sim_params['sigma'][z.str]['upper'] += sigma_nks
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def symmetric_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_k,
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        float dt,
        int Nz, int Nvz, int k,
        int charge,
        sim_params,
        z, vz
        ):

    """
    This routine assumes a time substep is positive as the Vlasov-Poisson system is irreversible
    when boundaries are present (i.e. the upper boundary cannot be periodic/open in this setup)

    this boundary condition must not inherently operate on a symmetric velocity grid,
    but it is more natural and consistent this way. Some rationale is given below:

    For example, if an exiting particle at the symmetry plane has
    a velocity V1 < 0, it is the physical case that there exists an entering partner particle from the
    complement domain with oppositely directly velocity -V1 > 0. If a symmetric velocity
    grid, i.e. [avx, bvx] = [-V, V], is used then it is always possible to factor in both
    exiting (V1 < 0) density with entering partners (-V1 > 0) provided that V1 is inside [-V,V].

    If we have an asymmetric case where [avx, bvx] where avx != -bvx, then in the handling
    of the symmetry plane, there would be cases where such an exiting particle having V1 < 0 would have
    its entering partner with velocity -V1 not be representable (i.e. if avx < V1 is in the velocity grid,
    but bvx < -V1). We permit options so that the user may elect to take the following alternative to
    a symmetric velocity grid, i.e.

    We could still argue that it would be appropriate to handle such situations where the entering partner
    velocity is off-grid (at too largely positive of a value for a given grid [avx, bvx])
    then to take the partner density to be zero and interpret the  situation as the "window" or
    control volume in phase space. Thus, any entering particle with off-grid velocities will not re-enter
    the domain as their phase space coordinates are not in the simulation domain.

    a CUTOFF upper boundary condition should be chosen for the distribution function on the associated
    velocity variable for any case. CAUTION: periodic BCs on velocity for the distribution function are not
    appropriate, as the reflection in velocity exiting grid point j -> vz.N - 1 - j will not successfully mirror
    to the exact grid point needed if vz.N != vz.Ngridpoints (as would be the case for periodic boundary conditions)
    """
    vzpostpointmesh = vzprepointmesh.copy() # by this point, vzpostpointmesh has not been assigned (is an array of zeroes), init to prepointmesh values
                                            # and replace exiting particles with their partner velocities if encountered when scanning of all [i,j]

    cdef int i, j
    if k == 0:
        for j in range(Nvz/2):
            for i in range(Nz):
                if zpostpointmesh[i,j] <= 0:
                    # prepare for partner remap in lib.convect_configuration.remap_step
                    zpostpointmesh[i,j] = -zpostpointmesh[i,j]
                    vzpostpointmesh[i,j] = Nvz - 1 - vzprepointmesh[i,j]
                    Uf_old[i,j] = -Uf_old[i,j]

    if k == 1:
        for j in range(Nvz / 2):
            for i in range(Nz):
                if zpostpointmesh[i,j] <= 0:
                    # prepare for partner remap in lib.convect_configuration.remap_step
                    zpostpointmesh[i,j] = -zpostpointmesh[i,j]
                    vzpostpointmesh[i,j] = Nvz - 1 - vzpostpointmesh[i,j]
                    Uf_old[i,j] = -Uf_old[i,j]

    z.postpointmesh[k,:,:] = zpostpointmesh
    vz.postpointmesh[k,:,:] = vzpostpointmesh

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def cutoff_lower_boundary(np.ndarray[DTYPE_t, ndim=2] f_k,
                          np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] < 0:
                f_old[i,j] = 0
                Uf_old[i,j] = 0
                zpostpointmesh[i,j] = 0 # set postpoint whose density/flux is zero
                                        # to a dummy gridpoint, e.g. lower boundary

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def cutoff_upper_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] > Nz - 1:
                f_old[i,j] = 0
                Uf_old[i,j] = 0

                zpostpointmesh[i,j] = Nz - 1 # set postpoint whose density/flux is zero
                                             # to a dummy gridpoint, e.g. upper boundary

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

#==============================================================================#
# DEPRECATED METHODS
#==============================================================================#

@cython.boundscheck(False)
def collector_lower_boundary_flux_integrated(
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
    cdef DTYPE_t sigma_nks = 0    # n = time step, k = 0 (nearest gridpoint) or 1 (contiguous gridpoint), s = species with charge q_s = int charge
                                  # total charge sigma is sum over k, integration over n, and sum over s
                                  # $\sigma_{nks} = q_s \int dv_x |v_{x,j}| chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
                                  #          \simeq q_s\sum_{i,j} \Delta v_x |v_{x,j}| chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
    cdef DTYPE_t vzwidth = vz.width

    if k == 0: # often referred to as "k1"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] == 0 and -Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [-1/2, 0]$

                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = -1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = 0 (on-grid)

                    elif zpostpointmesh[i,j] <= 0: # $k_{true}\in (-\infty , -1/2)$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k1, $\chi_{k_1} f_{i,j}^n$ where $\chi_{k_1} = (1 + U_{i,j}^n)$, $U_{i,j}^n <= 0$

                        # charge collection:
                        #
                        sigma_nks += (-vzprepointvaluemesh[i,j])*(f_old[i,j] + Uf_old[i,j])

                        # remove this same proportion ("f + Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        f_old[i,j] = 0
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = 0

    if k == 1: # often referred to as "k2"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] == -1 and -Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [-1/2, 0]$ (k1 = 0, k2 = -1/2)

                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = -1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = -1/2 -> -1 (wall)

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = -U_{i,j}^n$, $U_{i,j}^n <= 0$
                        # charge collection:
                        sigma_nks += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                        # note: minus signs have been cancelled on on (-vzprepointvaluemesh[i,j])*(-Uf_old[i,j])

                        # remove this same proportion ("-Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = 0

                    elif zpostpointmesh[i,j] <= -1: # $k_{true}\in (-\infty , -1/2)$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = -U_{i,j}^n$, $U_{i,j}^n <= 0$

                        # charge collection:
                        #
                        #    $$\sigma_{m\ell}^\alpha} \equiv
                        sigma_nks += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                        # note: minus signs have been cancelled on on (-vzprepointvaluemesh[i,j])*(-Uf_old[i,j])

                        # remove this same proportion ("-Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        Uf_old[i,j] = 0
                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = 0

    # include multiplicative factors from integration (vzwidth) and charge

    sigma_nks *= charge
    sigma_nks *= vzwidth # for integration in velocity space

    # update cumulative charge density, this sum below accomplishes the following required sums and integrals:
    #
    #    o  integrates over time with width dt = split_coeff*t.width
    #    o  sums over k after k1 and k2 both go through this routine
    #    o  sums over species s when fe and fi both go through this routine
    sim_params['sigma'][z.str]['lower'] += sigma_nks * dt
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def collector_upper_boundary_flux_integrated(
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
    cdef DTYPE_t sigma_nks = 0    # n = time step, k = 0 (nearest gridpoint) or 1 (contiguous gridpoint), s = species with charge q_s = int charge
                                  # total charge sigma is sum over k, integration over n, and sum over s
                                  # $\sigma_{nks} = q_s \int dv_x |v_{x,j}| chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
                                  #          \simeq q_s\sum_{i,j} \Delta v_x |v_{x,j}| chi_{k,i,j} f(x^n_i -> wall at time n+1,v^n_{x,j})
                                  #
    cdef DTYPE_t vzwidth = vz.width

    if k == 0: # often referred to as "k1"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] > 0:
                    if zpostpointmesh[i,j] == Nz-1 and Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [Nz-1,Nz-1/2]$

                        # maps on-grid k1 = Nz-1
                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = Nz-1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = Nz-1 (on-grid)

                    elif zpostpointmesh[i,j] >= Nz-1: # $k_{true}\in (Nz-1/2,\infty )$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k1, $\chi_{k_1} f_{i,j}^n$ where $\chi_{k_1} = (1 - U_{i,j}^n)$, $U_{i,j}^n >= 0$

                        # charge collection:
                        #
                        sigma_nks += vzprepointvaluemesh[i,j]*(f_old[i,j] - Uf_old[i,j])

                        # remove this same proportion ("f + Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        f_old[i,j] = 0
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = Nz-1

    if k == 1: # often referred to as "k2"
        for i in range(Nz):
            for j in range(Nvz):
                if vzprepointvaluemesh[i,j] < 0:
                    if zpostpointmesh[i,j] == Nz and Uf_old[i,j] <= 1/2.*f_old[i,j]: # then $k_{true}\in [-1/2, 0]$ (k1 = 0, k2 = -1/2)

                        # reproportion based on half-cell width to acknowledge the boundary at half width, i = Nz-1/2, see notebook s24
                        Uf_old[i,j] = 2*Uf_old[i,j] # maps to i = Nz-1/2 -> Nz (wall)

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = -U_{i,j}^n$, $U_{i,j}^n <= 0$
                        # charge collection:
                        sigma_nks += vzprepointvaluemesh[i,j] * Uf_old[i,j]

                        # remove this same proportion ("Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks
                        Uf_old[i,j] = 0

                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = Nz-1

                    elif zpostpointmesh[i,j] >= Nz: # $k_{true}\in (Nz-1/2, \infty )$, density packet hits wall or beyond --> collect charge

                        #collect proportion allocated to k2, $\chi_{k_2} f_{i,j}^n$ where $\chi_{k_2} = U_{i,j}^n$, $U_{i,j}^n <= 0$

                        # charge collection:
                        #
                        #    $$\sigma_{m\ell}^\alpha} \equiv
                        sigma_nks += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                        # note: minus signs have been cancelled on on (-vzprepointvaluemesh[i,j])*(-Uf_old[i,j])

                        # remove this same proportion ("-Uf") from the distribution function since it's efflux to the wall
                        # has been tallied in sigma_nks

                        Uf_old[i,j] = 0
                        # assign a dummy index for the zero contribution to remap to (evades an off-grid index reference error)
                        zpostpointmesh[i,j] = Nz-1

    # include multiplicative factors from integration (vzwidth) and charge

    sigma_nks *= charge
    sigma_nks *= vzwidth # for integration in velocity space

    # update cumulative charge density, this sum below accomplishes the following required sums and integrals:
    #
    #    o  integrates over time with width dt = split_coeff*t.width
    #    o  sums over k after k1 and k2 both go through this routine
    #    o  sums over species s when fe and fi both go through this routine
    sim_params['sigma'][z.str]['upper'] += sigma_nks * dt
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old
