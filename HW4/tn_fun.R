### function to sample from truncated normal distribution

tn.fun = function(n, x, mu, sig, a, b, rng_c, GPU) {
  if(GPU) {
    m = loadModule('tnC.ptx')
    
    # compute grid/block dims
    source('utility.R')
    dims = compute_grid(n)
    
    # call the kernel
    this = .cuda(m$tn_kernel, 'x' = x, as.integer(n), 
                 mu, sig, a, b, 
                 as.integer(rng_c), 2000L, 
                 gridDim = as.integer(dims$grid_dims), 
                 blockDim = as.integer(dims$block_dims), 
                 outputs = 'x')
    return(this)
  } else return(rtnorm(n, mu, sig, a, b)) # if not GPU
}