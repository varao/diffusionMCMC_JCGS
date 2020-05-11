The file wrapper.R runs both the proposed HMC sampler as well as the Euler-Maruyama and random-weight particle MCMC samplers.

The HMC sampler is implemented by the function sde_hmc() in the file sder.R. The function returns a list, which consists of:
1) theta: a vector of posterior samples of the diffusion parameter
2) xObsMat: a matrix of the posterior samples of the diffusion at the observation times.
3) x0 and xT: vectors of posterior samples of the diffusion at the end-points.
4) dfm: a dataframe of posterior samples of the diffusion skeleton. The dataframe has 3 columns:
        a) time: the Poisson times 
        b) omega: the diffusion evaluated on the Poisson times 
        c) chainInd: indicates which MCMC sample the row of the dataframe belongs to (since the skeleton has random length)
