#!/usr/bin/env python3
import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import anndata as ad
from scipy import sparse

name = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
print("GPU", name)
x = ad.AnnData(
    sparse.random(2000, 500, density=0.03, format="csr", random_state=0, dtype=np.float32)
)
rsc.get.anndata_to_GPU(x)
rsc.pp.scale(x, max_value=10, zero_center=False)
rsc.pp.pca(x, n_comps=20, zero_center=False)
rsc.pp.neighbors(x, n_neighbors=15, n_pcs=20, algorithm="cagra")
print("RSC_SMOKE=OK")
