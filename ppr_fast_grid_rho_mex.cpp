/**
 * @file ppr_fast_grid_rho_mex.cpp
 * Implement a seeded ppr clustering scheme that finds
 * the best cluster for all tolerances eps in an interval
 *
 * Call with debugflag = 1 to display parameter values
 * before/after each call to a major function.
 * 
 * This function automatically switches between dense and sparse computation
 * depending on the set.
 *
 * USAGE:
 * [bestset,cond,cut,vol] = ppr_fast_grid_rho_mex(A,set,alpha,eps_min,theta,rho,debugflag)
 *
 */

#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif 

#define MINSPARSE 10000


#include "mex.h"
#include <vector>
#include <assert.h>
#include <math.h>
#include <utility> // for pair sorting
#include <algorithm>

#include "sparse_maxshelf.hpp"
#include "maxshelf.hpp"
#include "sparse_and_dense_container.hpp"

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
mwIndex sr_degree(const sparserow *s, mwIndex u) {
    return (s->ai[u+1] - s->ai[u]);
}


struct greater2nd {
  template <typename P> bool operator() (const P& p1, const P& p2) {
    return p1.second > p2.second;
  }
};

struct bestset_stats {
  double conductance;
  mwIndex volume;
  mwIndex cut;
  std::vector<double> cuts_eps;
  std::vector<double> vols_eps;
  std::vector<double> setsizes_eps;
  std::vector<mwIndex> cluster;
  const sparserow* G;
  
  typedef std::vector< std::pair<mwIndex, double> > vertex_prob_type;
  vertex_prob_type prpairs;
  std::vector< unsigned int > denserank;
  
  bestset_stats(const sparserow* _G) 
  : conductance(std::numeric_limits<double>::max()), volume(0), cut(0),
          G(_G), denserank(0)
  {
    assert(G->n < (mwIndex)(std::numeric_limits<double>::max()) );
    prpairs.resize(std::min(G->n, (mwSize)MINSPARSE));
  }
  
  typedef google::dense_hash_map<mwIndex, unsigned int> rank_map;
   
  template <typename RankType>
  void _internal_sweep(const RankType &rank) {

    double curcond = 0.;
    mwIndex curvol = 0, curcut = 0;
    size_t i=0;
    double mincond = std::numeric_limits<double>::max();
    size_t mincondind=0;
    mwIndex mincut=0, minvol=0;

    // do the sweep
    mwIndex total_degree = G->ai[G->m];
    for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
      it!=itend; ++it, ++i) {
      mwIndex v = it->first;
      mwIndex deg = G->ai[v+1]-G->ai[v];
      mwIndex change = deg;
      unsigned int rankv = rank.get(v);
      for (mwIndex nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
        mwIndex nbr = G->aj[nzi];
        // need the function here to handle the two dense/sparse cases
        // and the fact that nbr might not be in the sparse rank
        if (rank.get(nbr) < rankv) {
          change -= 2;
        }    
      }
      curcut += change;
      curvol += deg;
      curcond = (double)curcut/
                            (double)std::min(curvol,total_degree-curvol);
      if (curvol >= total_degree-curvol) {
        break;
      }
      if (curcond < mincond) {
        mincondind = i;
        mincond = curcond;
        minvol = curvol;
        mincut = curcut;
      }
    }
       
    if (i>0) {
        cuts_eps.push_back(mincut);
        vols_eps.push_back(minvol);
        setsizes_eps.push_back(mincondind+1);
    } else {
        cuts_eps.push_back(0);
        vols_eps.push_back(0);
        setsizes_eps.push_back(0);
    }
    
    if (mincond < conductance) {
      conductance = mincond;
      volume = minvol;
      cut = mincut;
      cluster.clear();        
      i = 0;
      for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
        it!=itend && i<(mincondind+1); ++it, ++i) {
        cluster.push_back(it->first);
      }
        DEBUGPRINT(( "* new sweep: mincond=%f size=%i solnlen=%i\n", 
                            mincond, mincondind+1, prpairs.size()));
    } else {
        DEBUGPRINT(( "  new sweep: mincond=%f size=%i solnlen=%i\n", 
                            mincond, mincondind+1, prpairs.size()));
    }
  }
    
  template <typename SolutionType>
  void update_cluster_from_sweep(SolutionType p, std::vector<mwIndex>& p_nnzs) {  
    if (prpairs.capacity() < p_nnzs.size()) {
        prpairs.reserve( std::min( G->n, 
                                   std::max(10*prpairs.capacity(), p_nnzs.size())));
    }
    prpairs.resize(p_nnzs.size());
    for (size_t nzi=0; nzi<p_nnzs.size(); ++nzi) {
      mwIndex node_ind = p_nnzs[nzi];
      prpairs[nzi] = std::make_pair( node_ind, p[node_ind] );
    }

    // now we have to do the sweep over p in sorted order by value
    std::sort(prpairs.begin(), prpairs.end(), greater2nd());
    
    // TODO adjust this breakpoint
    // used 200000 as of 2/18/2014 based on timing with europe_osm
    // to match the original ppr_grid_rho_mex.cpp
    if (p_nnzs.size() >= (size_t)((double)(G->n)*0.05)) {
    //if (p_nnzs.size() >= 200000) {
        
        denserank.resize(G->n, (unsigned int)(G->n + 1));
       
        size_t i=0;
        for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
            it!=itend; ++it, ++i) {
          denserank[it->first] = i;
        }
        
        _internal_sweep(make_dense_array(&denserank[0], (size_t)G->n, (unsigned int)(G->n+1)));
    } else {
      rank_map rank;
      rank.set_empty_key((mwIndex)-1);
      
      
      size_t i=0;
      for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
        it!=itend; ++it, ++i) {
        rank[it->first] = i;
      }  

      _internal_sweep(make_const_sparse_array(rank, (unsigned int)(G->n+1)));
    }
  }
};

struct ppr_grid {
    
    const sparserow *G;
    const double alpha;
    const double theta;
    const double rho;
    std::vector<mwIndex> tracker;
    std::vector<mwIndex> soln_nnzs;
  
    const double eps_min;
    size_t npush;
    size_t nstep;
    size_t max_npush;
    double current_rmax;
    mwIndex top_shelf;
    std::vector<double> epsvec_vals;
    
    double *dense_soln_vec;
    
    ppr_grid(const sparserow* _G, const double _alpha, const double _eps_min,
             const double _theta, const double _rho)
    : G(_G), alpha(_alpha), theta(_theta), eps_min(_eps_min), rho(_rho)
    {}
            
  
    /** Clear the current value of epsilon. 
    */
    template <typename ShelfType, typename SolutionType>
    void clear_eps(ShelfType& residual, SolutionType& soln, 
              double cur_eps) 
    {    
        while ( current_rmax > cur_eps && npush < max_npush) {
            
            if ( residual.find_top_shelf(top_shelf) == 0 ){ break; } 
            
            mwIndex new_top_shelf = top_shelf; 
        
            double rij = 0.0; 
            mwIndex ri = 0;
            // STEP 1: pop element from top_shelf
            if( residual.pop_from_shelf(top_shelf,ri,rij) == 0 ){ break; }
            double rij_remainder = cur_eps*rho;
            rij = rij - rij_remainder;
            mwIndex dummy_shelf2 = residual.update(ri,rij_remainder); // handles the heap internally
//            if ( dummy_shelf2 < new_top_shelf ){ new_top_shelf = dummy_shelf2; }
            
            rij_res = cur_eps*rho;
            r.update(ri, rij_res ); // handles the heap internally
            rij = rij_temp - rij_res;

            
            
            // STEP 2: update soln vector
            if (soln[ri] == 0.0 ){
                soln_nnzs.push_back(ri);
                tracker.push_back(ri);
                tracker.push_back(nstep);
            } // this is first time ri is touched
            soln[ri] += rij;
        
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
            nstep++;
        }
    }
    
    template <typename ResidualType>
    void _init_residual(ResidualType& residual, std::vector<mwIndex>& set) {
        //const mwIndex n = G->n;
        for (size_t i = 0; i < set.size(); i++){
            // residual normalized to be stochastic, then degree-normalized
            mwIndex ind = set[i];
            assert(ind >= 0); assert(ind < G->n); // assert that "set" contains indices i: 1<=i<=n
            size_t setideg = sr_degree(G,ind);
            double val = 1.0/(double)(set.size()*(double)setideg);
            residual.update(ind,val);
        }
    }
    
    template <typename SparseSolution, typename SparseResidual, 
                typename DenseSolution, typename DenseResidual>
	void convert_to_dense(SparseSolution& cursoln, SparseResidual& curresid,
                           DenseSolution& newsoln, DenseResidual& newresid)
    {
        for (size_t i=0; i<soln_nnzs.size(); ++i) {
            newsoln[soln_nnzs[i]] = cursoln[soln_nnzs[i]];
        }
        
        typedef typename google::dense_hash_map<mwIndex, double>::const_iterator SparseIter;
        for (SparseIter it = curresid.values.begin(), itend = curresid.values.end(); 
                it != itend; ++it) {
            newresid.update(it->first, it->second);
        }
    }
      
    void compute(std::vector<mwIndex>& set, bestset_stats& stats) {
        
        // initial setup
        npush = 0;
        nstep = 0;
        top_shelf = 0;
        soln_nnzs.clear();
        tracker.clear();
        current_rmax = 0; // TODO avoid making this a class-wide variable
        
        max_npush = std::min( 
                (mwIndex)std::numeric_limits<int>::max() , 
                (mwIndex)(1/((1-alpha)*eps_min)) );
        
        const mwIndex epsvec_size =  1 + ceil(log2(eps_min/0.1)/log2(theta));
        epsvec_vals.resize( epsvec_size, 0.0 );
        epsvec_vals[0] = 0.1;
        for (mwIndex j = 1; j < epsvec_size; j++) { 
            epsvec_vals[j] = epsvec_vals[j-1]*theta; 
        }
        // end initial setup
        
        double cur_eps = epsvec_vals[0];
        mwIndex eps_num = 0;
        
        const size_t switch_threshold = (size_t)((double)(G->n)*0.05);
        const size_t minsparse = MINSPARSE;
        
        size_t last_sweep_step = 0; // track the step number every time we
        
        // start with the sparse computation
        sparse_max_shelf sparse_residual(theta, eps_min, 0.1, G->n, debugflag);
        google::dense_hash_map<mwIndex, double> sparse_soln_vec;
        sparse_soln_vec.set_empty_key(G->n);
        
        _init_residual(sparse_residual, set);
        
        current_rmax = sparse_residual.look_max(top_shelf);

        DEBUGPRINT(( "initial residual = %f\n", current_rmax ));
        
        if (G->n >= minsparse) {
            while ( current_rmax > eps_min && npush < max_npush) {
                DEBUGPRINT(( "nstep = %6i / %7i -- clearing eps %i = %e, current_rmax %e\n",
                    nstep, npush, eps_num, cur_eps, current_rmax ));
                clear_eps(sparse_residual, sparse_soln_vec, cur_eps); 
                eps_num += 1;
                cur_eps = epsvec_vals[eps_num];
                if (last_sweep_step < nstep) {
                    stats.update_cluster_from_sweep(sparse_soln_vec, soln_nnzs);
                    last_sweep_step = nstep;
                }
                // TODO tune this value
                if (sparse_residual.values.size() > switch_threshold) {
                    break;
                }
            }
        }
        
        if (current_rmax > eps_min) {
            //std::vector<double> dense_soln_vec(G->n, 0.0);
            mxArray* matlab_soln = mxCreateDoubleMatrix(G->n,1,mxREAL); // initialized to 0
            dense_soln_vec = mxGetPr(matlab_soln);
            
            max_shelf dense_residual(theta, eps_min, 0.1, G->n, debugflag);
            convert_to_dense(sparse_soln_vec, sparse_residual, 
                             dense_soln_vec, dense_residual);
        
            while ( current_rmax > eps_min && npush < max_npush) {
                clear_eps(dense_residual, dense_soln_vec, cur_eps); 
                // get new value of cur_eps
                for (size_t j = eps_num; j < epsvec_vals.size(); j++){
                    if ( epsvec_vals[j] < cur_eps){
                        eps_num = (mwIndex)j;
                        cur_eps = epsvec_vals[j];
                        break;
                    }
                }
                stats.update_cluster_from_sweep(dense_soln_vec, soln_nnzs);
            }
        }
    }
    
};
    
void compute_pagerank(sparserow* G, std::vector<mwIndex>& set, const double alpha,
                    const double eps_min, const double theta, 
                    bestset_stats& stats, 
                    double* soln_vec, 
                    std::vector<mwIndex>& tracker)
{
    // declare variables and set constants
    const mwIndex max_npush = std::min( 
                (mwIndex)std::numeric_limits<int>::max() , 
                (mwIndex)(1/((1-alpha)*eps_min)) );
    mwIndex npush = 0;
    mwIndex nstep = 0;
    const mwIndex n = G->n;

    // initialize error-tracking structures
    const mwIndex epsvec_size =  1 + ceil(log2(eps_min/0.1)/log2(theta));
    std::vector<double> epsvec_vals;
    epsvec_vals.resize( epsvec_size, 0.0 );
    epsvec_vals[0] = 0.1;
    for (mwIndex j = 1; j < epsvec_size; j++) { 
        epsvec_vals[j] = epsvec_vals[j-1]*theta; 
    }
    double cur_eps = epsvec_vals[0];
    mwIndex eps_num = 0;
       
    // initialize residual and solution
    //max_shelf residual(theta, eps_min, 0.1, n, debugflag);
    sparse_max_shelf residual(theta, eps_min, 0.1, n, debugflag);
    std::vector<mwIndex> soln_nnzs;

    for (mwIndex i = 0; i < set.size(); i++){
        // residual normalized to be stochastic, then degree-normalized
        mwIndex ind = set[i];
        assert(ind >= 0); assert(ind < n); // assert that "set" contains indices i: 1<=i<=n
        size_t setideg = sr_degree(G,ind);
        double val = 1.0/(double)(set.size()*(double)setideg);
        residual.update(ind,val);
    }
    
    mwIndex top_shelf = 0;
    double current_rmax = residual.look_max(top_shelf);
    DEBUGPRINT(( "initial residual = %f\n", current_rmax ));
    while ( current_rmax > eps_min && npush < max_npush) {
        // find current top_shelf
        if ( residual.find_top_shelf(top_shelf) == 0 ){ break; } 
        mwIndex new_top_shelf = top_shelf; 
        
        double rij = 0.0; 
        mwIndex ri = 0;
        // STEP 1: pop element from top_shelf
        if( residual.pop_from_shelf(top_shelf,ri,rij) == 0 ){ break; } 

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
            stats.update_cluster_from_sweep(soln_vec, soln_nnzs);
        }
        nstep++;
    }//END 'while'
}  // END compute_pagerank()



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
// [bestset,cond,cut,vol,eps_stats,sol_track] = ppr_grid_rho_mex(A,set,alpha,eps_min,theta,rho,debugflag)
void mexFunction(int num_outputs, mxArray* ptr_outputs[], int num_inputs, const mxArray* ptr_inputs[])
{
    if (num_inputs == 7){ debugflag = (int)mxGetScalar(ptr_inputs[6]); }
    DEBUGPRINT(("ppr_grid_rho_mex: preprocessing start: \n"));

    if (num_inputs < 2 || num_inputs > 7) {
        mexErrMsgIdAndTxt("ppr_grid_rho_mex:wrongNumberArguments",
                          "ppr_grid_rho_mex needs two to seven arguments, not %i", num_inputs);
    }

    const mxArray* mat = ptr_inputs[0];
    const mxArray* set = ptr_inputs[1];

    if ( mxIsSparse(mat) == false ){
        mexErrMsgIdAndTxt("ppr_grid_rho_mex:wrongInputMatrix",
                          "ppr_grid_rho_mex needs input 1 to be a sparse input matrix");
    }
    if ( mxGetM(mat) != mxGetN(mat) ){
        mexErrMsgIdAndTxt("ppr_grid_rho_mex:wrongInputMatrixDimensions",
                          "ppr_grid_rho_mex needs input 1 to be a square input matrix");
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
        mexErrMsgIdAndTxt("ppr_grid_rho_mex:wrongNumberOutputs",
                          "ppr_grid_rho_mex needs 0 to 7 outputs, not %i", num_outputs);
    }
    
    double eps_min = pow(10,-4);
    double alpha = .99;
    double theta = .8;
    double rho = 0.0;
    
    if (num_inputs >= 6) { rho = mxGetScalar(prhs[5]); }
    if (num_inputs >= 5) { theta = mxGetScalar(ptr_inputs[4]); }
    if (num_inputs >= 4) { eps_min = mxGetScalar(ptr_inputs[3]); }
    if (num_inputs >= 3) { alpha = mxGetScalar(ptr_inputs[2]); }


    if ( rho < 0 || rho >= 1.0 ){
        mexErrMsgIdAndTxt("ppr_fast_grid_rho_mex:wrongArgumentsRho",
                          "ppr_fast_grid_rho_mex needs 0 <= rho < 1, not %f", rho);
    }
    if ( eps_min <= 0 || eps_min > 1.0 ){
        mexErrMsgIdAndTxt("ppr_fast_grid_rho_mex:wrongArgumentsEps",
                          "ppr_fast_grid_rho_mex needs 0 < eps <= 1, not %f", eps_min);
    }
    if ( alpha <= 0 || alpha >= 1.0 ){
        mexErrMsgIdAndTxt("ppr_fast_grid_rho_mex:wrongArgumentsAlpha",
                          "ppr_fast_grid_rho_mex needs 0 < alpha < 1, not %f", alpha);
    }
    if ( theta <= 0 || theta >= 1.0 ){
        mexErrMsgIdAndTxt("ppr_fast_grid_rho_mex:wrongArgumentsTheta",
                          "ppr_fast_grid_rho_mex needs 0 < theta < 1, not %f", theta);
    }
        
        
    sparserow G;
    G.m = mxGetM(mat);
    G.n = mxGetN(mat);
    G.ai = mxGetJc(mat);
    G.aj = mxGetIr(mat);
    G.a = mxGetPr(mat);
    G.volume = (double)(G.ai[G.m]);

    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    bestset_stats stats(&G);    

    //mxArray* soln = mxCreateDoubleMatrix(G.n,1,mxREAL); // initialized to 0
    //double* solnptr = mxGetPr(soln);

    std::vector<mwIndex> soltracker;
    
    DEBUGPRINT(
            ("ppr_grid_rho_mex: inputs handled, call compute_pagerank "
             "(%i seeds, %.3f, %g, %g) \n", seeds.size(), alpha, eps_min, theta, rho));

    //compute_pagerank(&G, seeds, alpha, eps_min, theta, stats, solnptr, soltracker);

    ppr_grid computer(&G, alpha, eps_min, theta, rho);
    
    computer.compute(seeds, stats);
    soltracker = computer.tracker;
    
    DEBUGPRINT(("ppr_grid_rho_mex: compute_pagerank DONE \n"));
    
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
    
    DEBUGPRINT(("ppr_grid_rho_mex::  done with compute_pagerank \n" ));  
    *condptr = stats.conductance;
    *cutptr = stats.cut;
    *volptr = stats.volume;     
    if (num_outputs > 0) { // sets output "bestset" to the set of best conductance
        mxArray* cassign = mxCreateDoubleMatrix(stats.cluster.size(),1,mxREAL);
        ptr_outputs[0] = cassign;
        double *ci = mxGetPr(cassign);
        for (size_t i = 0; i < stats.cluster.size(); i++) {
            ci[i] = (double)(stats.cluster[i] + 1);
        }
    }   
}
