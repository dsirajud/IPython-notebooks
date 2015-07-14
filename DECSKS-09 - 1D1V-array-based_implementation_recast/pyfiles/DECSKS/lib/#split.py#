import DECSKS
import time

def scheme(
        f,
        t,x,v,
        n,
        sim_params
        ):
    """Steps through 1D-1V Vlasov with chosen splitting scheme.

    inputs:
    f -- (ndarray, dim=3) f(t,x,v) with f(0,x,v) = f_0 initialized
    x -- (instance) x, space
    v -- (instance) v, velocity
    t -- (instance) time
    n -- (int) current time step index, t^n
    sim_params -- (dict) simulation parameters

    outputs:
    f -- (ndarray, dim=3) f(n+1,x,v)
    """

    # retrieve sub-dictionary containing splitting coefficients and composition order
    splitting = sim_params['splitting']
    coeff = splitting['order']['coeffs']
    stage = splitting['order']['stages']
    tic = time.time()
    for s in range(len(stage)):
        split_coeff = splitting[coeff[s]][int(stage[s])]
        if s == 0: # on first pass, the previous time step (n - 1) needs to be
            if coeff[s] == 'a': # convect x
                for j in v.prepoints:
                    x.MCs   = x.generate_Lagrangian_mesh(x.prepointvalues, v.prepointvalues[j], split_coeff*t.width) 
                    f[n,:,j] = DECSKS.lib.convect.scheme(
                        f[n-1,:,j],
                        x,n,
                        sim_params)

            elif coeff[s] == 'b': # convect v
                E = DECSKS.lib.fieldsolvers.Gauss(sim_params['ni'], f, x, v, n-1) # calculate accelerations at time zero (n-1)
                a = -E
                for i in x.prepoints:
                    v.MCs   = v.generate_Lagrangian_mesh(v.prepointvalues, a[i], split_coeff*t.width)
                    f[n,i,:] = DECSKS.lib.convect.scheme(
                        f[n-1,i,:],
                        v,n,
                        sim_params)

        else: # each subsequent steps overwrites the previous step, all at time n until all split steps complete
            if coeff[s] == 'a': # convect x
                for j in v.prepoints:
                    x.MCs   = x.generate_Lagrangian_mesh(x.prepointvalues, v.prepointvalues[j], split_coeff*t.width)
                    f[n,:,j] = DECSKS.lib.convect.scheme(
                        f[n,:,j],
                        x,n,
                        sim_params)

            elif coeff[s] == 'b': # convect v
                E = DECSKS.lib.fieldsolvers.Gauss(sim_params['ni'], f, x, v, n)
                a = -E
                for i in x.prepoints:
                    v.MCs   = v.generate_Lagrangian_mesh(v.prepointvalues, a[i], split_coeff*t.width)
                    f[n,i,:] = DECSKS.lib.convect.scheme(
                        f[n,i,:],
                        v,n,
                        sim_params)

    toc = time.time()
    print "time step %d of %d completed in %g seconds" % (n,t.N, toc - tic)

    return f
