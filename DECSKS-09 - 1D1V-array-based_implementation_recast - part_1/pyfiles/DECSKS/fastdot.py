def _force_forder(x):
    """
    Converts arrays x to fortran order. Returns
    a tuple in the form (x, is_transposed).
    """
    if x.flags.c_contiguous:
        return (x.T, True)
    else:
        return (x, False)
        
def dot(A,B):
    """
    Uses blas libraries directly to perform dot product
    """
    from scipy import linalg
    
    A, trans_a = _force_forder(A)
    B, trans_b = _force_forder(B)
    gemm_dot = linalg.get_blas_funcs("gemm", arrays=(A,B))
    
    # gemm is implemented to compute: C = alpha*AB  + beta*C
    return gemm_dot(alpha=1.0, a=A, b=B, trans_a=trans_a, trans_b=trans_b)
