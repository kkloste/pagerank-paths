---
title: "Simultaneous computation of graph diffusion for multiple values of error tolerance"
layout: project
---

ppr-path: seeded PageRank solution paths for refined local clustering
==========================================

### Kyle Kloster
### David F. Gleich


_These are research codes and may not work for you._


Demo
--------

	compile ;
	load ./data/usps_3nn.mat;
	n = size(G,1);
	seed = randi(n);
	ppr_pathplot_rho(G,seed,'rho',0.1);
	
	
    
Reusable codes
--------------
* `ppr_paths_rho_mex.cpp` C++ MEX code for computing the full solution paths for inputs "A", "seeds", "epsmin", "alpha", and "rho".
* `ppr_path_rho.m` Matlab wrapper for the mex code `ppr_paths_rho_mex.cpp`.
* `ppr_fast_grid_mex.cpp` C++ MEX code for computing PPR on a grid of eps parameters. Takes inputs "A", "seeds", "epsmin", "alpha", and "theta" (determines the fineness of the grid, i.e. how many eps values there are between 1e-1 and epsmin).
* `pprgrow_path_comp.m` modifies pprgrow.m to compute 10,000 instances of pprgrow, for use in comparing our ppr_paths code with the naive method.
* `ppr_pathplot_rho.m`  Given a graph and a seed, this computes the solution paths and also outputs images of the pathplot and conductance curve for the solution paths.
* `compile.m` will compile the necessary mex files.
* `/util/` contains standard codes for computing conducance and sweep cuts, etc.

Our sparse data structure codes used in the ppr paths code above; these are hash-table based 

* `sparseheap.hpp` a heap implemented using three hashtables
* `sparselist.hpp` our plain wrapper for Google's sparsehash
* `sparserank.hpp` for bubble sorting
* `sparsevec.hpp`  similar to sparselist.
* `sparse_maxshelf.hpp` used for linear-time approximate sorting
* `sparse_and_dense_container.hpp` used to dynamically switch between sparse and full arrays depending on problem size.


Codes from others
-----------------

taken from pprgrow project from [Gleich & Seshadhri 2012]:

* `pprgrow_mex.cc` C++ MEX code for computing a set of best conductance via seeded personalized pagerank with input graph "A" and input parameters "seeds, eps, alpha".
*  `pprgrow.m` Matlab script that calls pprgrow_mex, for more convenient selection of parameters. Makes ~30 calls to pprgrow_mex for ~30 different error tolerances, eps, and returns best-conductance set found via sweep-cuts over each of the ~30 solution vectors.

* `sparsehash/`  Google's sparse hashtable package


---------

# Experiment notes


* datasets must be made symmetric and loopless.
* we include the netscience dataset and USPS dataset that we used. We also include the code we used for generating the adjacency matrices for the USPS-digits dataset, in `/data/making_data/`

Reproducing Paper Results
--------
1. All experiments require you to specify the directory in which the relevant datasets are saved.
2. Datasets should be converted to `.mat` format.
* Before executing any experiment scripts, `compile.m` must be run to compile all mex codes.
* Furthermore, the relevant datasets must be collected, and directories in the experiment scripts must be changed to point to the appropriate directories containing those datasets.

### Experiment scripts:
To reproduce figures from the paper, run the following codes:

* Figure 1: run `/experiments/netscience/netscience_eps.m`
* Figure 2: run `/experiments/netscience/exact_paths_ACL.m` (Note this took about 2 minutes to run on a 2014 macbook air.)
* Figure 3: run `/experiments/prpaths/prpath_netsci.m` (Note that with rho=0.9 this will take 30-70 minutes, even though computing the PageRank info itself takes less than a second. This is because we haven't yet optimized the plotting features for runtime.)
* Figure 4: run `/experiments/prpaths/prpath_fbA.m` (This could take almost 15 minutes, even though computing the PageRank info itself takes less than a second. This is because we haven't yet optimized the plotting features for runtime.)
* Figure 5: From `/experiments/senate/` execute:
	* Top left and right: `senate_paths_setup.m` (this calls `senate_paths.m`)
	* (a) through (f): `senate_layout.m`
* Figure 6: run `/experiments/prpaths/prpath_usps.m`

* Figure 7:
	1. Generate results by executing `/experiments/rho_scaling/rho-scaling-exps.sh` (this shell script will execute `rho_scaling.m` for all datasets.)
	2. Generate images by executing `/plotting/rho_scaling_plot.m`

* Figure 8: run `/experiments/rho_scaling/lj_anomaly.m`

* Table 2:
	1. Generate data by executing (from `/experiments/timing_experiments/`)
		* `path_fullgrow_exps.sh` to generate data for algorithm `multi diff`
		* `path_time_exps.sh` to generate data for algorithms `Single diff` and `ppr-path`
	2. Generate table data by executing `/plotting/path_grow_compare_ejam.m`

* Table 3 and 4:
	1. Generate data by executing `/experiments/timing_experiments/grid_experiments.sh` (this calls `grid_v_grow_new.m` for all datasets)
	2. Generate tables by executing `/plotting/grid_bi_column_ejam.m`	
	
