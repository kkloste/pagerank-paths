---
title: "Simultaneous computation of graph diffusion for multiple values of error tolerance"
layout: project
---

ppr-path: seeded PageRank solution paths for refined local clustering
==========================================

### Kyle Kloster
### David F. Gleich


_These are research codes and may not work for you._

Download
--------

(nothing yet)

Demo
--------
	(nothing yet)
    
Reusable codes
--------------
* `ppr_path_rho_mex.cpp` C++ MEX code for computing the full solution paths for inputs "A", "seeds", "epsmin", "alpha", and "rho".
* `ppr_path_rho.m` Matlab wrapper for the mex code ppr_path_rho_mex.cpp.
* `ppr_grid_fast_mex.cpp` C++ MEX code for computing PPR on a grid of eps parameters. Takes inputs "A", "seeds", "epsmin", "alpha", and "theta" (determines the fineness of the grid, i.e. how many eps values there are between 1e-1 and epsmin).
* `pprgrow_path_comp.m` modifies pprgrow.m to compute 10,000 instances of pprgrow, for use in comparing our ppr_paths code with the naive method.

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
