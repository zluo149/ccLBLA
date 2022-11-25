/*
 *
 * Copyright (c) 2022
 *
 * Author: Zhiwen Luo
 * <zluo2903@gmail.com>
 *
 * HPP-ccLBLA is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * HPP-ccLBLA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

/*
 * Based on the LDA Gibbs Sampler by Xuan-Hieu Phan, and ccLDA by Michael Paul
 *	http://gibbslda.sourceforge.net/
*/

#ifndef	_MODEL_H
#define	_MODEL_H

#include "constants.h"
#include "dataset.h"

using namespace std;

// LDA model
class model {
public:
    // fixed options
    string wordmapfile;		// file that contains word map [string -> integer id]
    string trainlogfile;	// training log file
    string tassign_suffix;	// suffix for topic assignment file
    string theta_suffix;	// suffix for theta file
    string phi_suffix;		// suffix for phi file
    string psi_suffix;		// suffix for psi file
    string others_suffix;	// suffix for file containing other parameters
    string twords_suffix;	// suffix for file containing words-per-topics

    string dir;			// model directory
    string dfile;		// data file
    string model_name;		// model name
    int model_status;		// model status:
				// MODEL_STATUS_UNKNOWN: unknown status
				// MODEL_STATUS_EST: estimating from scratch
				// MODEL_STATUS_ESTC: continue to estimate the model from a previous one
				// MODEL_STATUS_INF: do inference

    dataset * ptrndata;	// pointer to training dataset object
    dataset * pnewdata; // pointer to new dataset object

    mapid2word id2word; // word map [int => string]

    // --- model parameters and variables ---
    int M; // dataset size (i.e., number of docs)
    int V; // vocabulary size
    int C; // number collections
    int K; // number of topics
    double alpha, alpha0, beta, hpu, hpv, hps, hpt, gamma0, gamma1; // LDA hyperparameters
    double countU, countS;
    int niters; // number of Gibbs sampling iterations
    int liter; // the iteration at which the model was saved
    int savestep; // saving period
    int twords; // print out top words per each topic
    int withrawstrs;

    double ** pu;
    double * pv;

    double * p; // temp variable for sampling
    int ** z; // topic assignments for words, size M x doc.size()
    int ** zx;
    int * zc;
    int ** nw; // cwt[i][j]: number of instances of word/term i assigned to topic j, size V x K
    int ** nwF;
    int *** nwc;
    int *** nwcF;
    int *** nx;
    int *** nxF;
    int ** nd; // na[i][j]: number of words in document i assigned to topic j, size M x K
    int ** ndF;
    int ** nc;
    int ** ncF;
    int * wcsum;
    int * wcsumF;
    int * ncsum;
    int * ncsumF;
    int * nwsum; // nwsum[j]: total number of words assigned to topic j, size K
    int * nwsumF;
    int ** nwcsum;
    int ** nwcsumF;
    int ** nxsum;
    int ** nxsumF;
    int * ndsum; // nasum[i]: total number of words in document i, size M
    int * ndsumF;
    double ** theta; // theta: document-topic distributions, size M x K
    double ** phi; // phi: topic-word distributions, size K x V
    double ** phiF;
    double *** sigma;
    double *** sigmaF;
    double ** psi; // psi: topic-collection distributions, size K x C
    double *probC;

    //Differential privacy
    double epsilonL;
    double ClipBound;

    // for inference only
    int inf_liter;
    int newM;
    int newV;
    int ** newz;
    int ** newzx;
    int ** newnw;
    int *** newnwc;
    int ** newnd;
    int *** newnx;
    int ** newnxsum;
    int ** newnwcsum;
    int * newnwsum;
    int * newndsum;
    double ** newtheta;
    double *** newsigma;
    double ** newphi;
    double ** newpsi;
    // --------------------------------------

    model() {
	set_default_values();
    }

    ~model();

    // set default values for variables
    void set_default_values();

    // parse command line to get options
    int parse_args(int argc, char ** argv);

    // initialize the model
    int init(int argc, char ** argv);

    // load LDA model to continue estimating or to do inference
    int load_model(string model_name);

    // save LDA model to files
    // model_name.tassign: topic assignments for words in docs
    // model_name.theta: document-topic distributions
    // model_name.phi: topic-word distributions
    // model_name.others: containing other parameters of the model (alpha, beta, M, V, K)
    int save_model(string model_name);
    int save_model_tassign(string filename);
    int save_model_tassign2(string filename);
    int save_model_xassign(string filename);
    int save_model_theta(string filename);
    int save_model_phi(string filename);
    int save_model_psi(string filename);
    int save_model_others(string filename);
    //int save_model_twords(string filename);

    // saving inference outputs
    int save_inf_model(string model_name);
    int save_inf_model_tassign(string filename);
    int save_inf_model_newtheta(string filename);
    int save_inf_model_newpsi(string filename);
    int save_inf_model_newphi(string filename);
    int save_inf_model_others(string filename);
    //int save_inf_model_twords(string filename);

    // init for estimation
    int init_est();
    int init_estc();

    // estimate LDA model using Gibbs sampling
    void estimate();
    void cuda_estimate();
    int sampling(int m, int n);
    void compute_theta();
    void compute_phi();
    void compute_psi();

    // init for inference
    int init_inf();
    // inference for new (unseen) data based on the estimated LDA model
    void inference();
    int inf_sampling(int m, int n);
    void compute_newtheta();
    void compute_newphi();
    void compute_newpsi();

    double laplaceRandom(double m, double b);
};

#endif
