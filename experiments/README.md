
ppr-path: diffusion based clustering
=====================================

### Kyle Kloster
### David F. Gleich


Experiment notes
--------

* datasets must be made symmetric and loopless.
* we include code for generating the adjacency matrices for the USPS-digits dataset.

Reproducing Paper Results
--------
1. All experiments require you to specify the directory in which the relevant datasets are saved.
2. Datasets should be converted to `.mat` format.
* Before executing any experiment scripts, `compile.m` must be run to compile all mex codes.
* Furthermore, the relevant datasets must be collected, and directories in the experiment scripts must be changed to point to the appropriate directories containing those datasets.

### DONE:
* Figure 1: run `/experiments/netscience/netscience_eps.m`
* Figure 2: run `/experiments/netscience/exact_paths_ACL.m` (Note this took about 2 minutes to run on a 2014 macbook air.)
* Figure 3: run `/experiments/prpaths/prpath_netsci.m` (Note that with rho=0.9 this will take 30-70 minutes, even though computing the PageRank info itself takes less than a second. This is because we haven't yet optimized the plotting features for runtime.)
* Figure 4: run `/experiments/prpaths/prpath_fbA.m` (This could take almost 15 minutes, even though computing the PageRank info itself takes less than a second. This is because we haven't yet optimized the plotting features for runtime.)
* Figure 5: From `/experiments/senate/` execute:
	* Top left and right: `senate_paths_setup.m` (this calls `senate_paths.m`)
	* (a) through (f): `senate_layout.m`
* Figure 6: run `/experiments/prpaths/prpath_usps.m`
* Figure 7: run `/experiments/rho_scaling/lj_anomaly.m`

### CHECK

* Table 2:
	1. Generate data by executing (from /ppr-all/experiments/timing_experiments)
		* `path_fullgrow_exps.sh` to generate data for algorithm `mult`
		* `path_time_exps.sh` to generate data for algorithms `single` and `path`
	2. Generate table data by executing (from /ppr-all/results)
		* `path_grow_compare_ejam.m`
* Table 3 and 4:
	1. Generate data by executing (from /ppr-all/experiments/timing_experiments)
		* `new_grid_experiments.sh`
	2. Generate table data by executing (from /ppr-all/results)
		* `grid_bi_column_ejam.m`
		
		
* Exact paths Figure

* Rho scaling Figure

* new senate images?
* lj anomaly script?
* rho scaling scripts