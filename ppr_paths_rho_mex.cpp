/**
 * @file ppr_paths_rho_mex.cpp
 * Implement a seeded ppr clustering scheme that finds
 * the best cluster for all tolerances eps in an interval.
 * Returns information about solution vector, residual,
 * and best cluster at every push step, and every
 * new value of epsilon reached.
 *
 *  Call with debugflag = 1 to display parameter values
 * before/after each call to a major function
 *
 *
 * USAGE:
 * [step_stats,eps_stats,bestset] = ppr_paths_rho_mex(A,set,alpha,eps_min,rho,debugflag)
 *
 *
 * TO COMPILE:
 *
 * if ismac
 *      mex -O -largeArrayDims ppr_paths_rho_mex.cpp
 * else
 * mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims ppr_paths_rho_mex.cpp
 *
 *
 */

#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#include "sparsehash/dense_hash_map.h"
#include <vector>
#include <utility> // for pair sorting
#include <assert.h>
#include <algorithm>
#include <math.h>

#include "sparselist.hpp" // include our sparse hashtable functions
#include "sparseheap.hpp" // include our heap functions
#include "sparserank.hpp" // include our sorted-list functions
#include "sparsevec.hpp" // include our sparse hashtable functions

#include "mex.h"


#define DEBUGPRINT(x) do { if (debugflag>=1) { \
mexPrintf x; mexEvalString("drawnow"); } \
} while (0)

int debugflag = 0;


struct sparserow {
    mwSize n, m;
    double volume;
    mwIndex *ai;
    mwIndex *aj;
    double *a;
};


/**
 * Returns the degree of node u in sparse graph s
 */
mwIndex sr_degree(sparserow *s, mwIndex u) {
    return (s->ai[u+1] - s->ai[u]);
}


struct sweep_info {
    mwIndex num_sweeps;
    std::vector<mwIndex> rank_of_best_cond;
    std::vector<double> cond;
    std::vector<double> vol;
    std::vector<double> cut;
    double best_cond_global;
    mwIndex rank_of_bcond_global;
    double vol_of_bcond_global;
    double cut_of_bcond_global;

    double best_cond_this_sweep;
    mwIndex back_ind;
    
    sweep_info(size_t initial_size)     // constructor
    : num_sweeps(0), rank_of_best_cond(initial_size,0), cond(initial_size,0.0), vol(initial_size,0.0),
        cut(initial_size,0.0), best_cond_global(1.0), rank_of_bcond_global(0),
        vol_of_bcond_global(0.0), cut_of_bcond_global(0.0), best_cond_this_sweep(1.0),
        back_ind(0)
    { 

    }
};

//  ep_stats.update(nsteps, cur_eps, loc_bcond,cur_cut, cur_vol, (rank_of_bcond+1) );
struct eps_info {
    mwIndex num_epsilons;
    
    std::vector<double> epsilons;
    std::vector<double> conds;
    std::vector<double> cuts;
    std::vector<double> vols;
    std::vector<mwIndex> setsizes;
    std::vector<mwIndex> stepnums;
        
    void update(mwIndex stepn, double eps, double cond, double cut, double vol, mwIndex sets){
        epsilons.push_back(eps);
        conds.push_back(cond);
        cuts.push_back(cut);
        vols.push_back(vol);
        setsizes.push_back(sets);
        stepnums.push_back(stepn);
        num_epsilons++;
        return;
    }
    
    eps_info(size_t init_size)
    : num_epsilons(0)
    {
    }
};


struct rank_record {
    mwIndex lastval; //sentinal
    std::vector<int> starts;
    std::vector<int> ends;
    std::vector<mwIndex> nodes;
    std::vector<mwIndex> deg_of_pushed;
    std::vector<size_t> size_of_solvec;
    std::vector<size_t> size_of_r;
    std::vector<double> val_of_push;
    std::vector<double> global_bcond;

    mwIndex nrank_changes;
    mwIndex nrank_inserts;
    mwIndex nsteps; // this will be the size of the data in these vectors
    mwIndex size_for_best_cond; // this gives size of the vectors when
    // the global best conductance was attained
    
    void update_record(mwIndex start_val, mwIndex end_val, mwIndex node_val, mwIndex deg,
                    size_t r_sz, double val_pushed, double gbcond){
        starts.push_back(start_val);
        ends.push_back(end_val);
        nodes.push_back(node_val);
        deg_of_pushed.push_back(deg);
        size_of_solvec.push_back(nrank_inserts);
        size_of_r.push_back(r_sz);
        val_of_push.push_back(val_pushed);
        global_bcond.push_back(gbcond);
        nsteps++;
        nrank_changes++;
        if (start_val == nrank_inserts){ // a new entry was added
            nrank_inserts++;
        }
    }
    
    rank_record()
    : nrank_changes(0), nrank_inserts(0), nsteps(0), size_for_best_cond(0)
    {
        lastval = std::numeric_limits<mwIndex>::max();
    }
};


typedef google::dense_hash_map<mwIndex, mwIndex> rank_map;
    
/*
 * returns index of the best conductance in the vector swinfo.conductances
 */
bool resweep(mwIndex r_end, mwIndex r_start, sparserow* G,
            sparse_max_rank<mwIndex,double,size_t>& rankinfo, sweep_info& swinfo){

    // ensure sweep_info vectors are big enough
    if ( r_start >= swinfo.cut.size() ){
        swinfo.cut.resize((r_start+1)*2);
        swinfo.vol.resize((r_start+1)*2);
        swinfo.cond.resize((r_start+1)*2);
        swinfo.rank_of_best_cond.resize((r_start+1)*2);
    }
    double old_bcond = swinfo.best_cond_global;
    (swinfo.num_sweeps) += (r_start-r_end+1);
    double total_degree = G->volume;
    bool was_there_new_bcond = 0;
    mwIndex rank_of_best_cond = 0;

//  FAST WAY/*
    mwIndex gindex = rankinfo.rank_to_index(r_end);
    double deg = (double)sr_degree(G,gindex);
    std::vector<double> neighbors_ranked_less_than(r_start-r_end+1,0.0);
    std::vector<double> oldcut(r_start-r_end+1,0.0);
    std::vector<double> oldvol(r_start-r_end+1,0.0);

    // get rankings of neighbors of the shifted node
    for (mwIndex nzi = G->ai[gindex]; nzi < G->ai[gindex+1]; nzi++){
        mwIndex temp_gindex = G->aj[nzi];
        mwIndex temp_rank = rankinfo.index_to_rank(temp_gindex);
        if ( (temp_rank < r_end) && (temp_rank >= 0) ){ neighbors_ranked_less_than[0] += 1.0; }
        if ( (temp_rank > r_end) && (temp_rank <= r_start) ){
            neighbors_ranked_less_than[temp_rank-r_end] = 1.0;
        }
    }
    for (mwIndex j = 1; j <= (r_start-r_end); j++){
        neighbors_ranked_less_than[j] += neighbors_ranked_less_than[j-1];
    }
    
    // get old cut/vol information
    if (r_end == 0){
        oldcut[0] = 0.0;
        oldvol[0] = 0.0;
    }
    else{
        oldcut[0] = swinfo.cut[r_end-1];
        oldvol[0] = swinfo.vol[r_end-1];
        rank_of_best_cond = swinfo.rank_of_best_cond[r_end-1];
    }    
    for (mwIndex j = 1; j <= (r_start-r_end); j++){
        oldcut[j] = swinfo.cut[r_end-1+j];
        oldvol[j] = swinfo.vol[r_end-1+j];
    }

    // update volumes and cuts from r_end to r_start
    double cur_cond = 1.0;
    for (mwIndex j = 0; j <= (r_start-r_end); j++){
        double cur_vol = oldvol[j] + deg;
        swinfo.vol[r_end+j] = cur_vol;
        double cur_cut = oldcut[j] + deg - 2.0*neighbors_ranked_less_than[j];
        swinfo.cut[r_end+j] = cur_cut;
        
        if (cur_vol == 0.0 || cur_vol == total_degree) { cur_cond = 1.0; }
        else { cur_cond = cur_cut/std::min(cur_vol,total_degree-cur_vol); }
    }

    // finally, compute conductance values from r_end to r_start
    for (mwIndex cur_rank = r_end; cur_rank <= r_start; cur_rank++){
        double cur_cut = swinfo.cut[cur_rank];
        double cur_vol = swinfo.vol[cur_rank];
        if (cur_vol == 0.0 || cur_vol == total_degree) { cur_cond = 1.0; }
        else { cur_cond = cur_cut/std::min(cur_vol,total_degree-cur_vol); }
        swinfo.cond[cur_rank] = cur_cond;
    }

    // ... and update 'rank_of_best_cond' for all indices r_end to r_start
    if ( r_start > swinfo.back_ind ) { swinfo.back_ind = r_start; }
    for (mwIndex cur_rank = r_end; cur_rank <= swinfo.back_ind; cur_rank++){
        if ( swinfo.cond[cur_rank] < swinfo.cond[rank_of_best_cond] ){ rank_of_best_cond = cur_rank; }
        swinfo.rank_of_best_cond[cur_rank] = rank_of_best_cond;
    }
    rank_of_best_cond = swinfo.rank_of_best_cond[swinfo.back_ind];
    cur_cond = swinfo.cond[rank_of_best_cond];
    swinfo.best_cond_this_sweep = cur_cond;


    if (cur_cond < old_bcond){ // if current best_cond_this_sweep improves...
        swinfo.rank_of_bcond_global = rank_of_best_cond;
        swinfo.best_cond_global = cur_cond; // update best_cond_global
        swinfo.vol_of_bcond_global = swinfo.vol[rank_of_best_cond];
        swinfo.cut_of_bcond_global = swinfo.cut[rank_of_best_cond];
        was_there_new_bcond = 1; // signal that a new best_cond_global was found
    }
    return was_there_new_bcond;
} // END resweep()

           

/*
 * reorders an ordered list to reflect an update to the rankings
 */
mwIndex rank_permute(std::vector<mwIndex> &cluster, mwIndex r_end, mwIndex r_start)
{
    mwIndex temp_val = cluster[r_start];
    for (mwIndex ind = r_start; ind > r_end; ind--){ cluster[ind] = cluster[ind-1]; }
    cluster[r_end] = temp_val;
    return (r_start-r_end);
}



/**
 *  graphdiffseed inputs:
 *      G   -   adjacency matrix of an undirected graph
 *      set -   seed vector: the indices of a seed set of vertices
 *              around which cluster forms; normalized so
 *                  set[i] = 1/set.size(); )
 *  output:
 *      p = f(tP) * set
 *              with infinity-norm accuracy of eps * f(t)
 *              in the degree weighted norm
 *  parameters:
 *      t   - the value of t
 *      eps - the accuracy
 *      max_push_count - the total number of steps to run
 */

void graphdiffseed(sparserow* G, sparsevec& set, const double t, const double eps_min,
        const double rho, const mwIndex max_push_count, eps_info& ep_stats, rank_record& rkrecord,
        std::vector<mwIndex>& cluster )
{
    DEBUGPRINT(("ppr_all_mex::graphdiffseed()  BEGIN \n"));
    G->volume = (double)(G->ai[G->m]);
    mwIndex npush = 0;
    mwIndex nsteps = 0;
    double best_eps = 1.0;
    double cur_eps = 1.0;
    std::vector<double> epsilons;
    std::vector<double> conds;
//    mwIndex stagnant = 0;
    
    // ***** initialize residual, solution, and bookkeeping vectors
    sparse_max_heap<mwIndex,double,size_t> r(1000);
    sparse_max_rank<mwIndex,double,size_t> solvec(1000);
    for (sparsevec::map_type::iterator it=set.map.begin(),itend=set.map.end();
         it!=itend;++it) {
        r.update(it->first,it->second); // "update" handles the heap internally
    }
    sweep_info spstats(1000);
    
    cur_eps = r.look_max();

    DEBUGPRINT(("ppr_all_mex::graphdiffseed()  variables declared, begin WHILE loop \n"));

    while ( (npush < max_push_count) && (cur_eps > eps_min) ) {
        // STEP 1: pop top element off of heap
        double rij, rij_temp, rij_res;
        mwIndex ri = r.extractmax(rij_temp); // heap handles sorting internally
        double degofi = (double)sr_degree(G,ri);
        rij_res = cur_eps*rho;
        r.update(ri, rij_res ); // handles the heap internally
        rij = rij_temp - rij_res;


        // STEP 2: update soln vector
        bool new_bcond = 0;            
        size_t rank_start;
        size_t old_size = solvec.hsize;
        size_t rank_end = solvec.update(ri, rij, rank_start, debugflag); // handles sorting internally.
                // Sets rank_start to the rank ri had before it was updated.
                // Sets rank_end to the rank ri has after it was updated.

        // STEP 3: update sweeps for new solution vector
        if ( rank_start == old_size ){ // CASE (1): new entry
            new_bcond = resweep(rank_end, old_size, G, solvec, spstats);
            rkrecord.update_record(old_size, rank_end, ri, degofi, r.hsize, rij, spstats.best_cond_global);            
        }
        else if( (rank_start < old_size) && (rank_start > rank_end) ){ // CASE (2): existing entry changes rank
            new_bcond = resweep(rank_end, rank_start, G, solvec, spstats);
            rkrecord.update_record(rank_start, rank_end, ri, degofi, r.hsize, rij, spstats.best_cond_global);
        } 
        else {
            // CASE (3): no changes to sweep info, just resweep anyway.
            new_bcond = resweep(rank_end, rank_start, G, solvec, spstats);
            rkrecord.update_record(rank_start, rank_end, ri, degofi, r.hsize, rij, spstats.best_cond_global);
        }

        // STEP 4: update residual
        double update = t*rij;
        for (mwIndex nzi=G->ai[ri]; nzi < G->ai[ri+1]; ++nzi) {
            mwIndex v = G->aj[nzi];
            r.update(v,update/(double)sr_degree(G,v)); // handles the heap internally            
        }

        // STEP 5: update cut-set stats, check for convergence
        double cur_max = r.look_max();
        if (cur_max < cur_eps){ // we've reached a lower value of || ||_{inf},
            cur_eps = cur_max;  // so update cut stats for new val of cur_eps
            mwIndex rank_of_bcond = spstats.rank_of_best_cond[spstats.back_ind];            
            double loc_bcond = spstats.cond[rank_of_bcond];
            double cur_vol = spstats.vol[rank_of_bcond];
            double cur_cut = spstats.cut[rank_of_bcond];
            ep_stats.update(nsteps, cur_eps, loc_bcond,cur_cut, cur_vol, (rank_of_bcond+1) );
        }         
        if ( new_bcond == 1 ){ // new best_cond_global, so update
            best_eps = cur_eps;
            rkrecord.size_for_best_cond = rkrecord.nrank_changes;
        }
        nsteps++;        
        npush+=degofi;
    }//END 'while'
    DEBUGPRINT(("WHILE done \n"));

    
    //reconstruct bestcluster from the record of rank changes, rkrecord
    cluster.resize(rkrecord.nrank_inserts);
    mwIndex cluster_length = 0;
    mwIndex num_rank_swaps = 0;
    for (mwIndex j = 0; j < rkrecord.size_for_best_cond ; j++){
        mwIndex rs = rkrecord.starts[j];
        mwIndex re = rkrecord.ends[j];
        mwIndex rn = rkrecord.nodes[j];
        if (rs == rkrecord.size_of_solvec[j]){ // node rn added for the first time
            cluster[cluster_length] = rn;
            rs = cluster_length;
            cluster_length++;
        }
        num_rank_swaps += rank_permute(cluster, re, rs);   
    }    
    cluster.resize(spstats.rank_of_bcond_global+1); // delete nodes outside best cluster

}  // END graphdiffseed()



/** Cluster will contain a list of all the vertices in the cluster
 * @param set - the set of starting vertices to use
 * @param t - scaling factor in f(t*A)
 * @param eps - the solution tolerance eps
 * @param p - the solution vector
 * @param r - the residual vector
 * @param a - vector which supports .push_back to add vertices for the cluster
 * @param stats - a structure for statistics of the computation
 */

 void hypercluster_graphdiff_multiple(sparserow* G, const std::vector<mwIndex>& set,
                        double t, double eps, double rho, eps_info& ep_stats, rank_record& rkrecord,
                         std::vector<mwIndex>& cluster)
{
    // reset data
    sparsevec r; r.map.clear();
    DEBUGPRINT(("beginning of hypercluster_graphdiff_multiple() \n"));
    size_t maxdeg = 0;
    for (size_t i=0; i<set.size(); ++i) { //populate r with indices of "set"
        assert(set[i] >= 0); assert(set[i] < G->n); // assert that "set" contains indices i: 1<=i<=n
        size_t setideg = sr_degree(G,set[i]);
        r.map[set[i]] += 1.0/(double)(set.size()*(double)setideg);
        // r is normalized to be stochastic, then degree-normalized
    DEBUGPRINT(("i = %i \t set[i] = %i \t setideg = %i \n", i, set[i], setideg));
        maxdeg = std::max(maxdeg, setideg);
    }
    DEBUGPRINT(("at last, graphdiffseed: t=%f eps=%f \n", t, eps));
    
    const mwIndex max_npush = std::min( 
        (mwIndex)std::numeric_limits<int>::max() , (mwIndex)(1/((1-t)*eps)) );
    graphdiffseed(G, r, t, eps, rho, max_npush, ep_stats, rkrecord, cluster);

}  // END hyper_cluster_graphdiff_multiple()
            

void copy_array_to_index_vector(const mxArray* v, std::vector<mwIndex>& vec)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);
    
    vec.resize(n);
    
    for (size_t i=0; i<n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1, "Only positive integer elements allowed");
        vec[i] = (mwIndex)elem - 1;
    }
}  // END copy_array_to_index_vector()



// USAGE
// [step_stats,eps_stats,bestset] = ppr_paths_rho_mex(A,set,alpha,eps_min,rho,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 2 || nrhs > 6) {
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongNumberArguments",
                          "ppr_stats_mex needs two to six arguments, not %i", nrhs);
    }
    if (nrhs == 6) {
        debugflag = (int)mxGetScalar(prhs[5]);
    }
    DEBUGPRINT(("ppr_stats_mex: preprocessing start: \n"));
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongInputMatrix",
                          "ppr_stats_mex needs sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongInputMatrixDimensions",
                          "ppr_stats_mex needs square input matrix");
    }

    if ( nlhs > 3 ){
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongNumberOutputs",
                          "ppr_stats_mex needs 0 to 3 outputs, not %i", nlhs);
    }
    
    double eps = pow(10,-4);
    double alpha = 0.99;
    double rho = 0.0;
    
    if (nrhs >= 5) { rho = mxGetScalar(prhs[4]); }
    if (nrhs >= 4) { eps = mxGetScalar(prhs[3]); }
    if (nrhs >= 3) { alpha = mxGetScalar(prhs[2]); }

    if ( rho < 0 || rho >= 1.0 ){
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongArgumentsRho",
                          "ppr_stats_mex needs 0 <= rho < 1, not %f", rho);
    }
    if ( eps <= 0 || eps > 1.0 ){
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongArgumentsEps",
                          "ppr_stats_mex needs 0 < eps <= 1, not %f", eps);
    }
    if ( alpha <= 0 || alpha >= 1.0 ){
        mexErrMsgIdAndTxt("ppr_stats_mex:wrongArgumentsAlpha",
                          "ppr_stats_mex needs 0 < alpha < 1, not %f", alpha);
    }
    
    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);

    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    DEBUGPRINT(("ppr_stats_mex: preprocessing end: \n"));
    
    eps_info ep_stats(1000);
    rank_record rkrecord;
    std::vector<mwIndex> bestclus;

    
    DEBUGPRINT(("ppr_stats_mex: call to hypercluster_graphdiff() start\n"));    
    hypercluster_graphdiff_multiple(&r, seeds, alpha, eps, rho, ep_stats, rkrecord, bestclus);
    DEBUGPRINT(("ppr_stats_mex: call to hypercluster_graphdiff() DONE\n"));
    
    if (nlhs > 2) { // sets output "bestset" to the set of best conductance
        mxArray* cassign = mxCreateDoubleMatrix(bestclus.size(),1,mxREAL);
        plhs[2] = cassign;
        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<bestclus.size(); ++i) {
            ci[i] = (double)(bestclus[i] + 1);
        }
    }
    if (nlhs > 0) { // step_stats output
        mwIndex numsteps = rkrecord.nsteps;
        mxArray* cassign = mxCreateDoubleMatrix(numsteps,8,mxREAL);
        plhs[0] = cassign;
        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<numsteps; ++i) {
            ci[i] = (double)rkrecord.starts[i] + 1;
            ci[i + numsteps] = (double)rkrecord.ends[i] + 1;
            ci[i + 2*numsteps] = (double)rkrecord.nodes[i]+1; // adjust so index is correct in matlab
            ci[i + 3*numsteps] = (double)rkrecord.deg_of_pushed[i];
            ci[i + 4*numsteps] = (double)rkrecord.size_of_solvec[i];
            ci[i + 5*numsteps] = (double)rkrecord.size_of_r[i];
            ci[i + 6*numsteps] = rkrecord.val_of_push[i];
            ci[i + 7*numsteps] = rkrecord.global_bcond[i];
        }
    }
    if (nlhs > 1) { // outputs info about successive eps values
        mwIndex m = ep_stats.num_epsilons;
        mxArray* solvecptr = mxCreateDoubleMatrix(m,6,mxREAL);
        plhs[1] = solvecptr;
        double *ci = mxGetPr(solvecptr);
        for (mwIndex j = 0; j < m; j++) {
            ci[j] = ep_stats.epsilons[j];
            ci[j + m] = ep_stats.conds[j];
            ci[j + 2*m] = ep_stats.cuts[j];
            ci[j + 3*m] = ep_stats.vols[j];
            ci[j + 4*m] = (double)ep_stats.setsizes[j];
            ci[j + 5*m] = (double)ep_stats.stepnums[j];
        }
    }
}