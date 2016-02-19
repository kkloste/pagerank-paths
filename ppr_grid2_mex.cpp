/**
 * @file ppr_grid2_mex.cpp
 * Implement a seeded ppr clustering scheme that finds
 * the best cluster for all tolerances eps in an interval.
 *
 * This version of ppr_grid outputs the actual pagerank solution vector,
 * unlike the other versions. The pagerank vector is scaled bny degrees,
 * i.e. it outputs D^(-1)*x, where x is the actual solution.
 *
 * Call with debugflag = 1 to display parameter values
 * before/after each call to a major function
 *
 * USAGE:
 * [bestset,cond,cut,vol,eps_stats,sol_track,soln_vec] = ppr_grid2_mex(A,set,alpha,eps_min,theta,deg_weight_flag,debugflag)
 *
 *
 */

#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 


#include "mex.h"
#include <vector>
#include <assert.h>
#include <math.h>
#include <utility> // for pair sorting
#include <algorithm>

#include "maxshelf.hpp"
#include "sparsehash/dense_hash_map.h"
#include "sparsevec.hpp"


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


struct bestset_stats {
    double conductance;
    double volume;
    double cut;
    std::vector<double> cuts_eps;
    std::vector<double> vols_eps;
    std::vector<double> setsizes_eps;

    bestset_stats()    // constructor
    : conductance(0.0), volume(0.0), cut(0.0)
    {
    }
};

typedef google::dense_hash_map<mwIndex, mwIndex> rank_map;

struct greater2nd {
  template <typename P> bool operator() (const P& p1, const P& p2) {
    return p1.second > p2.second;
  }
};

void cluster_from_sweep(sparserow* G, double* p, std::vector<mwIndex>& p_nnzs,
      std::vector<mwIndex>& cluster, double *outcond, double* outvolume,
      double *outcut)
{
  typedef std::vector< std::pair<mwIndex, double> > vertex_prob_type;
  vertex_prob_type prpairs(p_nnzs.size());
  for (size_t nzi=0; nzi<p_nnzs.size(); ++nzi) {
    mwIndex node_ind = p_nnzs[nzi];
    prpairs[nzi] = std::make_pair( node_ind, p[node_ind] );
  }
  // now we have to do the sweep over p in sorted order by value

  std::sort(prpairs.begin(), prpairs.end(), greater2nd());
  
  // compute cutsize, volume, and conductance
  std::vector<double> conductance(prpairs.size());
  std::vector<mwIndex> volume(prpairs.size());
  std::vector<mwIndex> cutsize(prpairs.size());
  size_t i=0;
  rank_map rank;
  rank.set_empty_key((mwIndex)-1);
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    rank[it->first] = i;
  }
  mwIndex total_degree = G->ai[G->m];
  mwIndex curcutsize = 0;
  mwIndex curvolume = 0;
  i=0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    mwIndex v = it->first;
    mwIndex deg = G->ai[v+1]-G->ai[v];
    mwIndex change = deg;
    for (mwIndex nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
      mwIndex nbr = G->aj[nzi];
      if (rank.count(nbr) > 0) {
        if (rank[nbr] < rank[v]) {
          change -= 2;
        }
      }
    }
    curcutsize += change;
    curvolume += deg;
    volume[i] = curvolume;
    cutsize[i] = curcutsize;
    if (curvolume == 0 || total_degree-curvolume==0) {
      conductance[i] = 1;
    } else {
      conductance[i] = (double)curcutsize/
                        (double)std::min(curvolume,total_degree-curvolume);
    }
    if (curvolume >= total_degree-curvolume) {
        break;
    }
  }
  
  // we stopped the iteration when it finished, or when it hit target_vol
  size_t lastind = i;
  double mincond = std::numeric_limits<double>::max();
  size_t mincondind = 0; // set to zero so that we only add one vertex 
  for (i=0; i<lastind; i++) {
    if (conductance[i] < mincond) {
      mincond = conductance[i];
      mincondind = i;
    }
  }

  if (lastind == 0) {
    // add a case 
    mincond = 0.0;
  }
  i = 0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend && i<mincondind+1; ++it, ++i) {
    cluster.push_back(it->first);
  }
  if (outcond) { *outcond = mincond; }
  if (outvolume) { *outvolume = volume[mincondind]; }
  if (outcut) { *outcut = cutsize[mincondind]; }
}


void compute_pagerank(sparserow* G, sparsevec& set, const double alpha,
                    const double eps_min, const double theta, std::vector<mwIndex>& cluster,
                    bestset_stats* stats, double* soln_vec, std::vector<mwIndex>& tracker)
{
    // declare variables and set constants
    const mwIndex max_npush = std::min( (mwIndex)std::numeric_limits<int>::max() , (mwIndex)(1/((1-alpha)*eps_min)) );
    mwIndex npush = 0;
    mwIndex nstep = 0;
    const mwIndex n = G->n;
    std::vector<mwIndex> temp_cluster;
    double best_cond = 1.0;
    double best_vol = 0.0;
    double best_cut = 0.0;
    stats->volume = 0.0;
    stats->cut = 0.0;
    stats->conductance = 1.0;

    // initialize error-tracking structures
    const mwIndex epsvec_size =  1 + ceil(log2(eps_min/0.1)/log2(theta));
    std::vector<double> epsvec_vals;
    epsvec_vals.resize( epsvec_size, 0.0 );
    epsvec_vals[0] = 0.1;
    for (mwIndex j = 1; j < epsvec_size; j++){ epsvec_vals[j] = epsvec_vals[j-1]*theta; }
    double cur_eps = epsvec_vals[0];
    mwIndex eps_num = 0;
    stats->cuts_eps.resize(epsvec_size,0.0);
    stats->vols_eps.resize(epsvec_size,0.0);
    stats->setsizes_eps.resize(epsvec_size,0.0);
       
    // initialize residual and solution
    const double eps_max = 0.1;
    max_shelf residual(theta, eps_min, eps_max, n, debugflag);
    std::vector<mwIndex> soln_nnzs;
    for (sparsevec::map_type::iterator it=set.map.begin(),itend=set.map.end();
         it!=itend;++it) {
        residual.update(it->first,it->second); // "update" handles the heap internally
    } 
    
    mwIndex top_shelf = 0;
    double current_rmax = residual.look_max(top_shelf);
    DEBUGPRINT(( "initial residual = %f\n", current_rmax ));
    while ( current_rmax > eps_min && npush < max_npush) {
        if ( residual.find_top_shelf(top_shelf) == 0 ){ break; } // find current top_shelf
        mwIndex new_top_shelf = top_shelf; 
        double rij = 0.0; 
        mwIndex ri = 0;
        
        // STEP 1: pop element from top_shelf
        if( residual.pop_from_shelf(top_shelf,ri,rij) == 0 ){ break; } // top shelf emptied!

        // STEP 2: update soln vector
        if (soln_vec[ri] == 0.0 ){
            soln_nnzs.push_back(ri);
            tracker.push_back(ri);
            tracker.push_back(nstep);
        } // this is first time ri is touched
        soln_vec[ri] += rij;
        

        // STEP 3: update residual
        double update = alpha*rij;
        for (mwIndex nzi=G->ai[ri]; nzi < G->ai[ri+1]; ++nzi) {
            mwIndex v = G->aj[nzi];
            mwIndex dummy_shelf = residual.update(v,update/(double)sr_degree(G,v)); // handles the heap internally            
            if ( dummy_shelf < new_top_shelf ){ new_top_shelf = dummy_shelf; }
        }
        npush+=sr_degree(G,ri);

        // STEP 4: check if new epsilon value is satisfied, if so sweep
        residual.find_top_shelf(new_top_shelf);
        current_rmax = epsvec_vals[new_top_shelf];
        top_shelf = new_top_shelf;
        if (current_rmax <= cur_eps ){ // if error reaches new eps, sweep
            // get new value of cur_eps
            for (mwIndex j = eps_num; j < epsvec_size; j++){
                if ( epsvec_vals[j] < cur_eps){
                    eps_num = j;
                    cur_eps = epsvec_vals[j];
                    break;
                }
            }
            temp_cluster.clear();
            cluster_from_sweep(G, soln_vec, soln_nnzs, temp_cluster,  &stats->conductance, &stats->volume, &stats->cut );
            // update the per epsilon bestset stats
            stats->vols_eps[eps_num] = stats->volume;
            stats->cuts_eps[eps_num] = stats->cut;
            stats->setsizes_eps[eps_num] = cluster.size();
            if ( stats->conductance < best_cond ){ // new best_cond_global, so update
                best_cond = stats->conductance;
                best_vol = stats->volume;
                best_cut = stats->cut;
                cluster = temp_cluster;
            }
        }
        nstep++;
    }//END 'while'

    stats->volume = best_vol;
    stats->cut = best_cut;
    stats->conductance = best_cond;
    return;
}  // END compute_pagerank()



void transform_array_to_residual(const mxArray* v, sparsevec& vec, sparserow* G, bool neighb_flag, double alpha_scaling)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);
    double seeds_volume = 0.0;
    
    for (size_t i=0; i<n; ++i) {
        double temp_elem = p[i];
        mxAssert(temp_elem >= 1, "Only positive integer elements allowed");
        mwIndex elem = (mwIndex)temp_elem - 1;
        size_t setideg = sr_degree(G,elem);
        seeds_volume += (double) setideg;
    }
    for (size_t i=0; i<n; ++i) {
        double temp_elem = p[i];
        mxAssert(temp_elem >= 1, "Only positive integer elements allowed");
        mwIndex elem = (mwIndex)temp_elem - 1;
        size_t setideg = sr_degree(G,elem);
        if (neighb_flag == 1){ vec.map[elem] += alpha_scaling/seeds_volume; }
        if (neighb_flag != 1){ vec.map[elem] += alpha_scaling/((double)n*(double)setideg); }
        /* Note this code computes the solution to (I-alpha*P^T)y = D^-1*v
         * so "degree weighted" is implicitly performed by setting
         */
    }
    
}  // END transform_array_to_residual()


// USAGE
// [bestset,cond,cut,vol,eps_stats,sol_track,soln_vec] = ppr_grid2_mex(A,set,alpha,eps_min,theta,deg_weight_flag,debugflag)
void mexFunction(int num_outputs, mxArray* ptr_outputs[], int num_inputs, const mxArray* ptr_inputs[])
{
    if (num_inputs == 7){ debugflag = (int)mxGetScalar(ptr_inputs[6]); }
    DEBUGPRINT(("ppr_grid2_mex: preprocessing start: \n"));

    if (num_inputs < 2 || num_inputs > 7) {
        mexErrMsgIdAndTxt("ppr_grid2_mex:wrongNumberArguments",
                          "ppr_grid2_mex needs two to seven arguments, not %i", num_inputs);
    }

    const mxArray* mat = ptr_inputs[0];
    const mxArray* set = ptr_inputs[1];

    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("ppr_grid2_mex:wrongInputMatrix",
                          "ppr_grid2_mex needs input 1 to be a sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("ppr_grid2_mex:wrongInputMatrixDimensions",
                          "ppr_grid2_mex needs input 1 to be a square input matrix");
    }

    mxArray* cond = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(1,1,mxREAL);
    if (num_outputs > 1) { ptr_outputs[1] = cond; }
    if (num_outputs > 2) { ptr_outputs[2] = cut; }
    if (num_outputs > 3) { ptr_outputs[3] = vol; }
    double* condptr = mxGetPr(cond);
    double* cutptr = mxGetPr(cut);
    double* volptr = mxGetPr(vol);

    if ( num_outputs > 7 ){
        mexErrMsgIdAndTxt("ppr_grid2_mex:wrongNumberOutputs",
                          "ppr_grid2_mex needs 0 to 7 outputs, not %i", num_outputs);
    }
    
    double eps_min = pow(10,-4);
    double alpha = .99;
    double theta = .66;
    bool deg_weight_flag = 1; // degree weighting the seed vector is standard
    if (num_inputs >= 5) { theta = mxGetScalar(ptr_inputs[4]); }
    if (num_inputs >= 4) { eps_min = mxGetScalar(ptr_inputs[3]); }
    if (num_inputs >= 3) { alpha = mxGetScalar(ptr_inputs[2]); }
    
    sparserow G;
    G.m = mxGetM(mat);
    G.n = mxGetN(mat);
    G.ai = mxGetJc(mat);
    G.aj = mxGetIr(mat);
    G.a = mxGetPr(mat);
    G.volume = (double)(G.ai[G.m]);

    sparsevec seeds; seeds.map.clear();
    if (num_inputs >= 6){ deg_weight_flag = (int)mxGetScalar(ptr_inputs[5]); }
    // transform_array_to_residual( set, seeds, &G, deg_weight_flag, 1.0-alpha);
    // eps_min = eps_min*(1.0-alpha);
    transform_array_to_residual( set, seeds, &G, deg_weight_flag, 1.0);
    
    bestset_stats stats;    
    std::vector<mwIndex> bestclus;

    mxArray* soln = mxCreateDoubleMatrix(G.n,1,mxREAL); // initialized to 0
    double* solnptr = mxGetPr(soln);

    std::vector<mwIndex> soltracker;
    
    DEBUGPRINT(("ppr_grid2_mex: inputs handled, call compute_pagerank \n"));

    compute_pagerank(&G, seeds, alpha, eps_min, theta, bestclus, &stats, solnptr, soltracker);

    DEBUGPRINT(("ppr_grid2_mex: compute_pagerank DONE \n"));

    if (num_outputs > 6){
        ptr_outputs[6] = soln;
    }

    
    if (num_outputs > 5){
        mwIndex num_insert = (soltracker.size()/2);
        mxArray* stptr = mxCreateDoubleMatrix(num_insert,2,mxREAL); // initialized to 0
        ptr_outputs[5] = stptr;
        double* slptr = mxGetPr(stptr);
        for (mwIndex j = 0; j < num_insert; j++){
            slptr[j] = (double)soltracker[2*j];
            slptr[j+num_insert] = (double)soltracker[2*j+1];
        }
    }

    if (num_outputs > 4){
        mwIndex eps_n = stats.cuts_eps.size();
        mxArray* epassign = mxCreateDoubleMatrix(eps_n,3,mxREAL);
        ptr_outputs[4] = epassign;
        double* ep = mxGetPr(epassign);
        for (size_t j = 0; j < eps_n; j++){
            ep[j] = stats.setsizes_eps[j];
            ep[j+eps_n] = stats.cuts_eps[j];
            ep[j+2*eps_n] = stats.vols_eps[j];
        }
    }
    
    DEBUGPRINT(("ppr_grid2_mex::  done with compute_pagerank \n" ));
    *condptr = stats.conductance;
    *cutptr = stats.cut;
    *volptr = stats.volume;     
    if (num_outputs > 0) { // sets output "bestset" to the set of best conductance
        mxArray* cassign = mxCreateDoubleMatrix(bestclus.size(),1,mxREAL);
        ptr_outputs[0] = cassign;
        double *ci = mxGetPr(cassign);
        for (size_t i = 0; i < bestclus.size(); i++) {
            ci[i] = (double)(bestclus[i] + 1);
        }
    }   
}
