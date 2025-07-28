import numpy as np

def auto_cutoff(mat, q=0.25):
    """
    Filter out genes with expression below the q-th percentile across all cells. (Default is 0.25)
    """
    mat = mat.X
    v = [arr[arr > 0] for x in mat for arr in [x.toarray()] if not np.all(arr <= 0)]
    v = np.concatenate(v)
    c = np.percentile(v, q * 100) 
    print("c:", c)
    return c