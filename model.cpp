/*
 *
 * Copyright (c) 2022
 *
 * Author: Zhiwen Luo
 * <zluo2903@gmail.com>
 *
 * ccLBLA is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * ccLBLA is distributed in the hope that it will be useful, but
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cassert>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "model.h"
#include "math.h"
#include <ctime>


using namespace std;

model::~model() {
    if (p) {
	     delete p;
    }

    if (ptrndata) {
	     delete ptrndata;
    }

    if (pnewdata) {
	     delete pnewdata;
    }

    if (z) {
	     for (int m = 0; m < M; m++) {
	        if (z[m]) {
		          delete z[m];
	        }
	     }
    }

    if (zx) {
	     for (int m = 0; m < M; m++) {
	        if (zx[m]) {
		          delete zx[m];
	        }
	     }
    }

    if (zc){
      delete[] zc;
    }

    if (nw) {
	     for (int w = 0; w < V; w++) {
	        if (nw[w]) {
		          delete nw[w];
	           }
	     }
   }

   if (nx){
     for (int x = 0; x < 2; x++){
       for (int c = 0; c < C; c++)
       {
         if (nx[x][c])
           delete nx[x][c];
       }

       if(nx[x])
         delete nx[x];
     }
   }

   if (nwc) {
      for (int c = 0; c < C; c++) {
         for (int w = 0; w < V; w++) {
            if (nwc[c][w]) {
                delete nwc[c][w];
            }
         }

        if (nwc[c]) {
          delete nwc[c];
        }
    }
 }

   if (nxsum) {
     for (int c = 0; c < C; c++) {
       if(nxsum[c])
          delete nxsum[c];
    }
   }

    if (nd) {
	     for (int m = 0; m < M; m++) {
	    if (nd[m]) {
		delete nd[m];
	    }
	}
    }

    if (nwsum) {
	delete[] nwsum;
    }

    if (nwcsum) {
	for (int c = 0; c < C; c++) {
	    if (nwcsum[c]) {
		delete nwcsum[c];
	    }
	}
    }

    if (ndsum) {
	delete[] ndsum;
    }

    if (theta) {
	for (int m = 0; m < M; m++) {
	    if (theta[m]) {
		delete theta[m];
	    }
	}
    }

    if (phi) {
	for (int k = 0; k < K; k++) {
	    if (phi[k]) {
		delete phi[k];
	    }
	}
    }


    if (sigma) {
      for (int c = 0; c < C; c++) {
      if (sigma[c]) {
	     for (int k = 0; k < K; k++) {
	    if (sigma[c][k]) {
		  delete sigma[c][k];
	    }
	   }
    }}}


    if (psi) {
	for (int k = 0; k < K; k++) {
	    if (psi[k]) {
		delete psi[k];
	    }
	}
    }

    // only for inference
    if (newz) {
	for (int m = 0; m < newM; m++) {
	    if (newz[m]) {
		delete newz[m];
	    }
	}
    }

    if (newzx) {
	for (int m = 0; m < newM; m++) {
	    if (newzx[m]) {
		delete newzx[m];
	    }
	}
    }

    if (newnw) {
	for (int w = 0; w < newV; w++) {
	    if (newnw[w]) {
		delete newnw[w];
	    }
	}
    }

    if (newnx){
      for (int x = 0; x < 2; x++){
        for (int c = 0; c < C; c++)
        {
          if (newnx[x][c])
            delete newnx[x][c];
        }

        if(newnx[x])
          delete newnx[x];
      }
    }

    if (newnwc) {
	     for (int c = 0; c < C; c++) {
	        for (int w = 0; w < newV; w++) {
	           if (newnwc[c][w]) {
		             delete newnwc[c][w];
	           }
	        }

	       if (newnwc[c]) {
           delete newnwc[c];
         }
	   }
  }

    if (newnd) {
	for (int m = 0; m < newM; m++) {
	    if (newnd[m]) {
		delete newnd[m];
	    }
	}
    }

    if (newnwsum) {
	delete[] newnwsum;
    }

    if (newndsum) {
	delete[] newndsum;
    }



    if (newnwcsum) {
	for (int c = 0; c < C; c++) {
	    if (newnwcsum[c]) {
		delete newnwcsum[c];
	    }
	}
    }

    if (newnxsum) {
     for (int c = 0; c < C; c++) {
       if(newnxsum[c])
          delete newnxsum[c];
    }
   }

   if (newtheta) {
 for (int m = 0; m < newM; m++) {
     if (newtheta[m]) {
   delete newtheta[m];
     }
 }
   }

    if (newphi) {
	for (int k = 0; k < K; k++) {
	    if (newphi[k]) {
		delete newphi[k];
	    }
	}
    }


    if (newsigma) {
      for (int c = 0; c < C; c++) {
      if (newsigma[c]) {
	     for (int k = 0; k < K; k++) {
	    if (newsigma[c][k]) {
		  delete newsigma[c][k];
	    }
	   }
    }}}


    if (newpsi) {
	for (int k = 0; k < K; k++) {
	    if (newpsi[k]) {
		delete newpsi[k];
	    }
	}
    }

    if (probC) delete [] probC;

    if (pu) {
	     for (int c = 0; c < C; c++) {
	        if (pu[c]) {
		          delete pu[c];
	        }
	     }
    }

    if (pv) {
		    delete pv;
    }

}

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    trainlogfile = "trainlog.txt";
    tassign_suffix = ".tassign";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    psi_suffix = ".psi";
    others_suffix = ".others";
    twords_suffix = ".twords";

    dir = "./";
    dfile = "trndocs.dat";
    model_name = "model-final";
    model_status = MODEL_STATUS_UNKNOWN;

    ptrndata = NULL;
    pnewdata = NULL;

    M = 0;
    V = 0;
    K = 100;
    C = 0;

    alpha = 0.1;
    beta = 0.01;
    hpu = 0.1;
    hpv = 0.1;
    hps = 0.1;
    hpt = 0.1;

    countU = 0.0;
    countS = 0.0;

    gamma0 = 1.0;
    gamma1 = 1.0;

    niters = 500;
    liter = 0;
    savestep = 200;
    twords = 0;
    withrawstrs = 0;

    p = NULL;
    z = NULL;
    zx = NULL;
    nw = NULL;
    nd = NULL;
    nwsum = NULL;
    ndsum = NULL;
    theta = NULL;
    phi = NULL;
    psi = NULL;
    probC = NULL;

    pu = NULL;
    pv = NULL;

    newM = 0;
    newV = 0;
    newz = NULL;
    newzx = NULL;
    newnw = NULL;
    newnwc = NULL;
    newnd = NULL;
    newnwsum = NULL;
    newnwcsum = NULL;
    newndsum = NULL;
    newtheta = NULL;
    newphi = NULL;
    newsigma = NULL;
    newpsi = NULL;
}

int model::parse_args(int argc, char ** argv) {
  printf("parse args\n");
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {
    // call parse_args
    if (parse_args(argc, argv)) {
	return 1;
    }

    printf("init\n");

    if (model_status == MODEL_STATUS_EST) {
	// estimating the model from scratch
	if (init_est()) {
	    return 1;
	}

    } else if (model_status == MODEL_STATUS_ESTC) {
	// estimating the model from a previously estimated one
	if (init_estc()) {
	    return 1;
	}

    } else if (model_status == MODEL_STATUS_INF) {
	// do inference
	if (init_inf()) {
	    return 1;
	}
    }

    printf("init0\n");

    return 0;
}

int model::load_model(string model_name) {
    int i, j;

    printf("load model\n");

    string filename = dir + model_name + tassign_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	     printf("Cannot open file %d to load model!\n", filename.c_str());
	     return 1;
    }

    char buff[BUFF_SIZE_LONG];
    string line;

    // allocate memory for z and ptrndata
    z = new int*[M];
    zx = new int*[M];
    zc = new int[M];
    ptrndata = new dataset(M);
    ptrndata->V = V;
    printf("%d\n",M);

    for (i = 0; i < M; i++) {

	     char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	     if (!pointer) {
	        printf("Invalid word-topic assignment file, check the number of docs!\n");
	        return 1;
	     }

    	line = buff;
    	strtokenizer strtok(line, " \t\r\n");
    	int length = strtok.count_tokens();

    	vector<int> words;
    	vector<int> topics;
    	vector<int> routes;
    	vector<int> classes;
    	for (j = 0; j < length; j++) {
    	    string token = strtok.token(j);

    	    strtokenizer tok(token, ":");
    	    if (tok.count_tokens() != 4) {
    		      printf("Invalid word-topic assignment line!\n");
    		      return 1;
    	    }

    	    words.push_back(atoi(tok.token(0).c_str()));
    	    topics.push_back(atoi(tok.token(1).c_str()));
    	    routes.push_back(atoi(tok.token(2).c_str()));
    	    classes.push_back(atoi(tok.token(3).c_str()));
          zc[i] = atoi(tok.token(3).c_str());
    	}

    	// allocate and add new document to the corpus
    	document * pdoc = new document(words, zc[i]);
    	ptrndata->add_doc(pdoc, i);

    	// assign values for z
    	z[i] = new int[topics.size()];
    	zx[i] = new int[routes.size()];
    	for (j = 0; j < topics.size(); j++) {
    	    z[i][j] = topics[j];
    	    zx[i][j] = routes[j];
    	    zc[i] = classes[j];
    	};
    }



    printf("load model3\n");
    fclose(fin);

    return 0;
}

int model::save_model(string model_name) {

    if (save_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }

    if (save_model_tassign2(dir + model_name + tassign_suffix + "2")) {
	return 1;
    }

    if (save_model_xassign(dir + model_name + ".xassign")) {
	return 1;
    }

    if (save_model_others(dir + model_name + others_suffix)) {
	return 1;
    }

    if (save_model_theta(dir + model_name + theta_suffix)) {
	return 1;
    }

    if (save_model_phi(dir)) {
	return 1;
    }

    if (save_model_psi(dir + model_name + psi_suffix)) {
	     return 1;
    }

    /*
    if (twords > 0) {
	     if (save_model_twords(dir + model_name + twords_suffix)) {
	        return 1;
	     }
    }
    */

    return 0;
}

int model::save_model_tassign(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {
	for (j = 0; j < ptrndata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d:%d:%d ", ptrndata->docs[i]->words[j], z[i][j], zx[i][j], zc[i]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);

    return 0;
}

int model::save_model_tassign2(string filename) {
    int i, j;
    mapid2word::iterator it;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {
       fprintf(fout, "%d ", ptrndata->docs[i]->cat);
	     for (j = 0; j < ptrndata->docs[i]->length; j++) {
          it = id2word.find(ptrndata->docs[i]->words[j]);
          if (it != id2word.end()) {
	           fprintf(fout, "%s:%d:%d ", (it->second).c_str(), z[i][j], zx[i][j]);
          }
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);

    return 0;
}

int model::save_model_xassign(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {
	for (j = 0; j < ptrndata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d:%d ", ptrndata->docs[i]->words[j], zx[i][j], ptrndata->docs[i]->cat);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);

    return 0;
}

int model::save_model_theta(string filename) {

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }

    for (int i = 0; i < M; i++) {
	     for (int j = 0; j < K; j++) {
	        fprintf(fout, "%f ", theta[i][j]);
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);
    return 0;
}

int model::save_model_phi(string filebase) {

    string filename = filebase+"words.phi";
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file to save!\n");
	     return 1;
    }

    for (int i = 0; i < K; i++) {
	     for (int j = 0; j < V; j++) {
	        fprintf(fout, "%f ", phi[i][j]);
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);

    string sc;
    for (int c = 0; c < C; c++) {
      std::stringstream ss;
      ss << c;
      filename = filebase+"words"+ss.str().c_str()+".phi";

      fout = fopen(filename.c_str(), "w");
      if (!fout) {
	       printf("Cannot open file to save!\n");
	       return 1;
      }

      double prob;
      for (int i = 0; i < K; i++) {
	       for (int j = 0; j < V; j++) {
	          prob = sigma[c][i][j];
	          fprintf(fout, "%f ", prob);
	       }
	     fprintf(fout, "\n");
      }

      fclose(fout);
    }

    return 0;
}

int model::save_model_psi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
       return 1;
    }


    for (int i = 0; i < K; i++) {
	     for (int j = 0; j < C; j++) {
	        fprintf(fout, "%f ", psi[i][j]);
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);

    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }

    //fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);

    fclose(fout);
    return 0;
}


int model::save_inf_model(string model_name) {
    /*
    if (save_inf_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }

    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	return 1;
    }

    if (save_inf_model_newphi(dir + model_name + phi_suffix)) {
	return 1;
    }

    if (twords > 0) {
	if (save_inf_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    */

    printf("end save\n");
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	     return 1;
    }

    printf("end save2\n");
    if (save_inf_model_newpsi(dir + model_name + ".psi")) {
	     return 1;
    }

    printf("end save3\n");
    if (save_inf_model_newphi(dir + model_name + phi_suffix)) {
	     return 1;
    }

    if (save_inf_model_others(dir + model_name + others_suffix)) {
	     return 1;
    }

    printf("end\n");

    return 0;
}

int model::save_inf_model_tassign(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < pnewdata->M; i++) {
	     for (j = 0; j < pnewdata->docs[i]->length; j++) {
	        fprintf(fout, "%d:%d ", pnewdata->docs[i]->words[j], newz[i][j]);
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;
    printf("end save1\n");

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }
    printf("end save2\n");

    for (i = 0; i < newM; i++) {
	     for (j = 0; j < K; j++) {
	        fprintf(fout, "%f ", newtheta[i][j]);
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);
    return 0;
}

int model::save_inf_model_newpsi(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }

	  for (j = 0; j < K; j++) {
      for (i = 0; i < C; i++) {
	       fprintf(fout, "%f ", newpsi[j][i]);
	    }
	    fprintf(fout, "\n");
    }

    fclose(fout);
    return 0;
}

int model::save_inf_model_newphi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }

    for (int i = 0; i < K; i++) {
	     for (int j = 0; j < newV; j++) {
	        fprintf(fout, "%f ", newphi[i][j]);
	     }
	     fprintf(fout, "\n");
    }

    fclose(fout);
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	     printf("Cannot open file %s to save!\n", filename.c_str());
	     return 1;
    }

    //fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", newM);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "liter=%d\n", inf_liter);

    fclose(fout);

    return 0;
}

/*
int model::save_inf_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    if (twords > newV) {
	twords = newV;
    }
    mapid2word::iterator it;
    map<int, int>::iterator _it;

    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < newV; w++) {
	    word_prob.first = w;
	    word_prob.second = newphi[k][w];
	    words_probs.push_back(word_prob);
	}

        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);

	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    _it = pnewdata->_id2id.find(words_probs[i].first);
	    if (_it == pnewdata->_id2id.end()) {
		continue;
	    }
	    it = id2word.find(_it->second);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }

    fclose(fout);

    return 0;
}
*/


int model::init_est() {
    int m, n, w, k;
    p = new double[K];

    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile)) {
        printf("Fail to read training data!\n");
        return 1;
    }

    // + allocate memory and assign values for variables
    M = ptrndata->M;
    printf("M = %d\n", M);
    V = ptrndata->V;
    printf("V = %d\n", V);
    C = ptrndata->C;
    printf("C = %d\n", C);

    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }

    nwc = new int**[C];
    for (int c = 0; c < C; c++) {
      nwc[c] = new int*[V];
      for (w = 0; w < V; w++) {
        nwc[c][w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nwc[c][w][k] = 0;
        }
      }
    }


    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nx = new int**[2];
    for (int x = 0; x < 2; x++) {
      nx[x] = new int*[C];
      for (int c = 0; c < C; c++) {
        nx[x][c] = new int[K];
        for (k = 0; k < K; k++) {
          nx[x][c][k] = 0;
        }
      }
    }

    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }

    nwcsum = new int*[C];
    for (int c = 0; c < C; c++) {
      nwcsum[c] = new int[K];
      for (k = 0; k < K; k++) {
	     nwcsum[c][k] = 0;
      }
    }


    nxsum = new int*[C];
    for (int c = 0; c < C; c++) {
      nxsum[c] = new int[K];
      for (k = 0; k < K; k++) {
        nxsum[c][k] = 0;
      }
    }

    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	     ndsum[m] = 0;
    }

    pu = new double*[C];
    for (int c = 0; c < C; c++) {
        pu[c] = new double[K];
        for (int k = 0; k < K; k++) {
    	    pu[c][k] = (double)(k)/(K+1);
        }
    }

    for (int k = 0; k < K-1; k++){
      countU += pu[0][k];
    }


    pv = new double[V];
    for (int v = 0; v < V; v++) {
      pv[v] = min(0.1, (double)1/(V+1));
    }

    for (int v = 0; v < V-1; v++){
      countS += pv[v];
    }

    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    zx = new int*[M];
    zc = new int[M];

    for (m = 0; m < ptrndata->M; m++) {
	     int N = ptrndata->docs[m]->length;
       //printf("The document%d size is %d\n", m, ptrndata->docs[m]->length);
	     z[m] = new int[N];
	     zx[m] = new int[N];
       zc[m] = ptrndata->docs[m]->cat;

	     int c = ptrndata->docs[m]->cat;

        // initialize for z
        for (n = 0; n < N; n++) {
    	    int topic = -1;
          while (topic == -1) {
            topic = (int)(((double)random() / RAND_MAX) * K);
          }
          int route = (int)(((double)random() / RAND_MAX) * 2);

    	    z[m][n] = topic;
    	    zx[m][n] = route;

          //printf("%d",ptrndata->docs[m]->words[n]);

    	    // number of instances of word i assigned to topic j
    	    if (route == 0) {
    	     nw[ptrndata->docs[m]->words[n]][topic] += 1;
           nwsum[topic] += 1;
    	    } else {
    	     nwc[c][ptrndata->docs[m]->words[n]][topic] += 1;
           nwcsum[c][topic] += 1;
    	    }

          // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    nx[route][c][topic] += 1;
    	    nxsum[c][topic] += 1;

        }
        //printf("\n");
        // total number of words in document i
        ndsum[m] = N;
    }

    theta = new double*[M];
    for (m = 0; m < M; m++) {
      theta[m] = new double[K];
    }

    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }

    sigma = new double**[C];
    for (int c = 0; c < C; c++) {
      sigma[c] = new double*[K];
      for (k = 0; k < K; k++) {
        sigma[c][k] = new double[V];
      }
    }

    psi = new double*[K];
    for (k = 0; k < K; k++) {
        psi[k] = new double[C];
    }

    //save_model(utils::generate_model_name(100));

    return 0;
}

int model::init_estc() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
	printf("Fail to load word-topic assignmetn file of the model!\n");
	return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }

    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }

    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;

	// assign values for nw, nd, nwsum, and ndsum
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];

    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        }
        // total number of words in document i
        ndsum[m] = N;
    }


    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }

    return 0;
}


void model::estimate() {
    if (twords > 0) {
	     // print out top words per topic
	     dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    printf("Sampling %d iterations!\n", niters);

    int last_iter = liter;
    printf("iteration = %d\n", niters + last_iter);
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
	     printf("Iteration %d ...\n", liter);

    	for (int m = 0; m < M; m++) {
    	    for (int n = 0; n < ptrndata->docs[m]->length; n++) {
        		int topic = sampling(m, n);
        		z[m][n] = topic;
        	}
    	}

      if (savestep > 0) {
        if (liter % savestep == 0) {
          // saving the model
          printf("Saving the model at iteration %d ...\n", liter);
          compute_theta();
          compute_phi();
          compute_psi();
        }
      }
    }

    printf("Gibbs sampling completed!\n");
    printf("Saving the final model!\n");
    compute_theta();
    compute_phi();
    compute_psi();
    liter--;

    save_model("new_Model"); //utils::generate_model_name(-1)
    printf("saved model\n");

}

int model::sampling(int m, int n) {
    int route = zx[m][n];
    int topic = z[m][n];
    int c = ptrndata->docs[m]->cat;
    int w = ptrndata->docs[m]->words[n];

    // remove z_i from the count variables
    if (route == 0) {
      nw[w][topic] -= 1;
      nwsum[topic] -= 1;
      nx[0][c][topic] -= 1;
    }
    else {
      nwc[c][w][topic] -= 1;
      nwcsum[c][topic] -= 1;
      nx[1][c][topic] -= 1;
    }
    nd[m][topic] -= 1;
    ndsum[m] -= 1;
    nxsum[c][topic] -= 1;

    //test doc-topic math equations
    double Kalpha = K * alpha;
    double Vbeta = V * beta;

    double totalprob = 0;
    double prob;
    // sample a topic

    double constant = (hpu + ndsum[m] - nd[m][K-1])/(countU + ndsum[m] - nd[m][K-1]);

    if (route == 0) {

      for (int k = 0; k < K; k++) {

        prob = ((nw[w][k] + pv[w]) / (hps + hpt + nwsum[k]) * (hps + nwsum[k] - nw[V-1][k])/(countS + nwsum[k] - nw[V-1][k])) *
            ((nd[m][k] + pu[c][k]) / (hpu + hpv + ndsum[m]) *constant);     // (countU + ndsum[m] - nd[m][K-1]) * (hpu + ndsum[m] - nd[m][K-1]) / (hpu + hpv + ndsum[m]));


        p[k] = prob;
        totalprob += prob;
      }

    }
    else {

      for (int k = 0; k < K; k++) {

        prob = ((nwc[c][w][k] + pv[w]) / (hps + hpt + nwcsum[c][k]) * (hps + nwcsum[c][k] - nwc[c][V-1][k]) / (countS + nwcsum[c][k] - nwc[c][V-1][k])) * ///
            ((nd[m][k] + pu[c][k])  / (hpu + hpv + ndsum[m]) *constant); /// (countU + ndsum[m] - nd[m][K-1]) * (hpu + ndsum[m] - nd[m][K-1])

        p[k] = prob;
		    totalprob += prob;
      }

    }

    double u = ((double)random() / RAND_MAX) * totalprob;

    prob = 0;
    for (int k = 0; k < K; k++) {
      prob += p[k];

      if (prob > u) {
        topic = k;
        break;
      }
    }

    // sample a route
    double prob0 = 0.0;
    double prob1 = 0.0;


    prob0 = ((nw[w][topic] + pv[w]) / (hps + hpt + nwsum[topic]) * (hps + nwsum[topic] - nw[V-1][topic]) / (countS + nwsum[topic] - nw[V-1][topic])) *
                 (nx[0][c][topic] + gamma0) / (nxsum[c][topic] + gamma0 + gamma1);

    prob1 = ((nwc[c][w][topic] + pv[w]) / (hps + hpt + nwcsum[c][topic]) * (hps + nwcsum[c][topic] - nwc[c][V-1][topic]) / (countS + nwcsum[c][topic] - nwc[c][V-1][topic])) *
                 (nx[1][c][topic] + gamma1) / (nxsum[c][topic] + gamma0 + gamma1);

    u = ((double)random() / RAND_MAX) * (prob0+prob1);

    if (u > prob0) route = 1;
    else route = 0;

    // add newly estimated values to count variables
    if (route == 0) {
      nw[w][topic] += 1;
      nwsum[topic] += 1;
      nx[0][c][topic] += 1;
    }
    else {
      nwc[c][w][topic] += 1;
      nwcsum[c][topic] += 1;
      nx[1][c][topic] += 1;
    }
    //again

    nd[m][topic] += 1;
    ndsum[m] += 1;
    nxsum[c][topic] += 1;

    zx[m][n] = route;

    return topic;
}

void model::compute_theta() {
  for (int m = 0; m < M; m++) {
    int c = zc[m];
    double constant = (hpu + ndsum[m] - nd[m][K-1]) / (countU + ndsum[m] - nd[m][K-1]);
    for (int k = 0; k < K; k++) {

      theta[m][k] = (nd[m][k] + pu[c][k]) / (hpu + hpv + ndsum[m]) *constant;
    }
  }
}

void model::compute_phi() {
  for (int k = 0; k < K; k++) {
     double constant1 = (hps + nwsum[k] - nw[V-1][k]) / (countS + nwsum[k] - nw[V-1][k]);
  	 for (int w = 0; w < V; w++) {

        phi[k][w] = (nw[w][k] + pv[w]) / (hps + hpt + nwsum[k]) * constant1;

  	    for (int c = 0; c < C; c++) {

          sigma[c][k][w] = (nwc[c][w][k] + pv[w]) / (countS + nwcsum[c][k] - nwc[c][V-1][k]) * (hps + nwcsum[c][k] - nwc[c][V-1][k]) / (hps + hpt + nwcsum[c][k]);
        }
  	 }
  }
}

void model::compute_psi() {
    for (int k = 0; k < K; k++) {
    	for (int c = 0; c < C; c++) {
	       psi[k][c] = (nx[1][c][k] + gamma1) / (nxsum[c][k] + gamma0 + gamma1);
	    }
    }
}

int model::init_inf() {
    // estimating the model from a previously estimated one
    int m, n, w, k;
    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
	     printf("Fail to load word-topic assignmetn file of the model!\n");
	     return 1;
    }
    printf("init load model\n");

    M = ptrndata->M;
    printf("M = %d\n", M);
    V = ptrndata->V;
    printf("V = %d\n", V);
    C = 15;
    printf("C = %d\n", C);

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }

    nwc = new int**[C];
    for (int c = 0; c < C; c++) {
      nwc[c] = new int*[V];
      for (w = 0; w < V; w++) {
        nwc[c][w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nwc[c][w][k] = 0;
        }
      }
    }


    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K+1];
        for (k = 0; k < K+1; k++) {
    	    nd[m][k] = 0;
        }
    }

    nx = new int**[2];
    for (int x = 0; x < 2; x++) {
      nx[x] = new int*[C];
      for (int c = 0; c < C; c++) {
        nx[x][c] = new int[K];
        for (k = 0; k < K; k++) {
          nx[x][c][k] = 0;
        }
      }
    }

    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	     nwsum[k] = 0;
    }

    nwcsum = new int*[C];
    for (int c = 0; c < C; c++) {
      nwcsum[c] = new int[K];
      for (k = 0; k < K; k++) {
	       nwcsum[c][k] = 0;
      }
    }

    nxsum = new int*[C];
    for (int c = 0; c < C; c++) {
      nxsum[c] = new int[K];
      for (k = 0; k < K; k++) {
        nxsum[c][k] = 0;
      }
    }

    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	     ndsum[m] = 0;
    }

    pu = new double*[C];
    for (int c = 0; c < C; c++) {
        pu[c] = new double[K];
        for (int k = 0; k < K; k++) {
    	    pu[c][k] = (double)(k+1)/(K+1);
        }
    }

    pv = new double[V];
    for (int v = 0; v < V; v++) {
      pv[v] = (double)v/(V+1);
    }

    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }

    sigma = new double**[C];
    for (int c = 0; c < C; c++) {
      sigma[c] = new double*[K];
      for (k = 0; k < K; k++) {
        sigma[c][k] = new double[V];
      }
    }

    psi = new double*[K];
    for (k = 0; k < K; k++) {
        psi[k] = new double[C];
    }

    for (int k = 0; k < K-1; k++){
      countU += pu[0][k];
    }

    for (int v = 0; v < V-1; v++){
      countS += pv[v];
    }

    for (m = 0; m < ptrndata->M; m++) {

        int c = zc[m];
	      int N = ptrndata->docs[m]->length;

	      // assign values for nw, nd, nwsum, and ndsum
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    int route= zx[m][n];

    	    // number of instances of word i assigned to topic j
    	    if (route == 0) {
    	     nw[w][topic] += 1;
           nwsum[topic] += 1;
    	    } else {
    	     nwc[c][w][topic] += 1;
           nwcsum[c][topic] += 1;
    	    }

          // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    nx[route][c][topic] += 1;
    	    nxsum[c][topic] += 1;

        }
        // total number of words in document i
        ndsum[m] = N;
    }


    // read new data for inference
    pnewdata = new dataset;
    if (withrawstrs) {
      printf("withrawstrs\n");
	    if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	    }
    }else {
	     if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	       }
   }

    newM = pnewdata->M;
    printf("newM = %d\n", newM);
    newV = pnewdata->V;
    printf("newV = %d\n", newV);

    newnw = new int*[newV];
    for (w = 0; w < newV; w++) {
        newnw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnw[w][k] = 0;
        }
    }

    newnwc = new int**[C];
    for (int c = 0; c < C; c++) {
      newnwc[c] = new int*[newV];
      for (w = 0; w < newV; w++) {
        newnwc[c][w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnwc[c][w][k] = 0;
        }
      }
    }

    newnd = new int*[newM];
    for (m = 0; m < newM; m++) {
        newnd[m] = new int[K+1];
        for (k = 0; k < K+1; k++) {
    	    newnd[m][k] = 0;
        }
    }

	  newnx = new int**[2];
	  for (int x = 0; x < 2; x++) {
      newnx[x] = new int*[C];
      for (int c = 0; c < C; c++) {
        newnx[x][c] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnx[x][c][k] = 0;
        }
      }
    }

    newnwsum = new int[K];
    for (k = 0; k < K; k++) {
	     newnwsum[k] = 0;
    }

    newnwcsum = new int*[C];
    for (int c = 0; c < C; c++) {
      newnwcsum[c] = new int[K];
      for (k = 0; k < K; k++) {
	     newnwcsum[c][k] = 0;
      }
    }

    newndsum = new int[newM];
    for (m = 0; m < newM; m++) {
	     newndsum[m] = 0;
    }


    newnxsum = new int*[C];
    for (int c = 0; c < C; c++) {
      newnxsum[c] = new int[K];
      for (k = 0; k < K; k++) {
	       newnxsum[c][k] = 0;
      }
    }

    probC = new double[C];
    for (int c = 0; c < C; c++)
      probC[c] = 0.0;

    srandom(time(0)); // initialize for random number generation
    newz = new int*[newM];
    newzx = new int*[newM];

    double total = 0.0;
    for (m = 0; m < pnewdata->M; m++) {
      int c = pnewdata->_docs[m]->cat;
	    int N = pnewdata->docs[m]->length;
      probC[c] += 1;
      total += 1;

	    newz[m] = new int[N];
	    newzx[m] = new int[N];

	    // assign values for nw, nd, nwsum, and ndsum
      for (n = 0; n < N; n++) {
  	    int w = pnewdata->docs[m]->words[n];
  	    int _w = pnewdata->_docs[m]->words[n];
  	    int topic = (int)(((double)random() / RAND_MAX) * K);
  	    int route = (int)(((double)random() / RAND_MAX) * 2);
  	    newz[m][n] = topic;
  	    newzx[m][n] = route;

        if (route == 0) {
    	     newnw[_w][topic] += 1;
           newnwsum[topic] += 1;
    	  }else{
    	     newnwc[c][_w][topic] += 1;
           newnwcsum[c][topic] += 1;
    	  }

  	    newnd[m][topic] += 1;
  	    newnx[route][c][topic] += 1;
  	    newnxsum[c][topic] += 1;
      }
      newndsum[m] = N;
    }

    for (int c = 0; c<C; c++){
      probC[c] /= total;
      printf("probC[%d] = %f\n", c, probC[c]);
    }

    double temp = 0.0;
    for (int c = 0; c<C; c++){
      temp += probC[c];
      printf("all probC = %f\n",  temp);
    }

    printf("init inf6\n");

    newtheta = new double*[newM];
    for (m = 0; m < newM; m++) {
        newtheta[m] = new double[K];
    }

    newphi = new double*[K];
    for (k = 0; k < K; k++) {
        newphi[k] = new double[newV];
    }

    newsigma = new double**[C];
    for (int c = 0; c < C; c++) {
      newsigma[c] = new double*[K];
      for (k = 0; k < K; k++) {
        newsigma[c][k] = new double[newV];
      }
    }

    newpsi = new double*[K];
    for (k = 0; k < K; k++) {
        newpsi[k] = new double[C];
    }

    return 0;
}

void model::inference() {
    if (twords > 0) {
	     // print out top words per topic
	     dataset::read_wordmap(dir + wordmapfile, &id2word);
    }
    printf("Sampling %d iterations for inference!\n", niters);

    for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
	     printf("Iteration %d ...\n", inf_liter);
    	// for all newz_i
    	for (int m = 0; m < newM; m++) {
    	    for (int n = 0; n < pnewdata->docs[m]->length; n++) {
    		// (newz_i = newz[m][n])
    		// sample from p(z_i|z_-i, w)
    		int topic = inf_sampling(m, n);
    		newz[m][n] = topic;
    	    }
    	}
    }

    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    compute_newtheta();
    compute_newphi();
    compute_newpsi();
    //printf("end\n");
    compute_phi();
    compute_psi();
    //printf("end\n");
    inf_liter--;
    save_inf_model(dfile);
    printf("end\n");
}

int model::inf_sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = newz[m][n];
    int route = newzx[m][n];
    int c = pnewdata->_docs[m]->cat;
    //printf("C = %d\n", c);
    int w = pnewdata->docs[m]->words[n];
    int _w = pnewdata->_docs[m]->words[n];
    if (route == 0) {
      newnw[_w][topic] -= 1;
      newnwsum[topic] -= 1;
      newnx[0][c][topic] -= 1;
    }
    else {
      newnwc[c][_w][topic] -= 1;
      newnwcsum[c][topic] -= 1;
      newnx[1][c][topic] -= 1;
    }
    newnd[m][topic] -= 1;
    newndsum[m] -= 1;
    newnxsum[c][topic] -= 1;

    double Kalpha = K * alpha;
    double Vbeta = V * beta;

    double gamma0 = 1.0;
    double gamma1 = 1.0;

    double totalprob = 0.0;
    double prob = 0.0;

    // sample a topic
    double constant = (hpu + newndsum[m] - newnd[m][K-1]) / (countU + newndsum[m] - newnd[m][K-1]);

    if (route == 0) {
      for (int k = 0; k < K; k++) {

        prob = ((nw[w][k] + newnw[_w][k] + pv[w]) / (hps + hpt + nwsum[k] + newnwsum[k]) * (hps + nwsum[k] + newnwsum[k] - nw[V-1][k] - newnw[newV-1][k])/(countS + nwsum[k] + newnwsum[k] - nw[V-1][k] - newnw[newV-1][k])) *
            ((newnd[m][k] + pu[c][k]) / (hpu + hpv + newndsum[m]) *constant);

        p[k] = prob;
        totalprob += prob;
      }
    }
    else {
      for (int k = 0; k < K; k++) {

        prob = ((nwc[c][w][k] + newnwc[c][_w][k] + pv[w]) / (hps + hpt + nwcsum[c][k] + newnwcsum[c][k]) * (hps + nwcsum[c][k] + newnwcsum[c][k] - nwc[c][V-1][k] - newnwc[c][newV-1][k])/(countS + nwcsum[c][k] + newnwcsum[c][k] - nwc[c][V-1][k] - newnwc[c][newV-1][k])) *
            ((newnd[m][k] + pu[c][k]) / (hpu + hpv + newndsum[m]) *constant);

        p[k] = prob;
		    totalprob += prob;
      }
    }

    double u = ((double)random() / RAND_MAX) * totalprob;

    //printf("u %f\n", u);
    prob = 0;
    for (int k = 0; k < K; k++) {
      prob += p[k];

      if (prob > u) {
        topic = k;
        break;
      }
    }

    // sample a route
   double prob0 = 0.0;
   double prob1 = 0.0;

   prob0 = ((nw[w][topic] + newnw[_w][topic] + pv[w]) / (hps + hpt + nwsum[topic] + newnwsum[topic]) * (hps + nwsum[topic] + newnwsum[topic] - nw[V-1][topic] - newnw[newV-1][topic]) / (countS + nwsum[topic] + newnwsum[topic] - nw[V-1][topic] - newnw[newV-1][topic])) *
                (nx[0][c][topic] + newnx[0][c][topic] +  gamma0) / (nxsum[c][topic] + newnxsum[c][topic] + gamma0 + gamma1);

   prob1 = ((nwc[c][w][topic] + newnwc[c][_w][topic] + pv[w]) / (hps + hpt + nwcsum[c][topic] + newnwcsum[c][topic]) * (hps + nwcsum[c][topic] + newnwcsum[c][topic] - nwc[c][V-1][topic] - newnwc[c][newV-1][topic]) / (countS + nwcsum[c][topic] + newnwcsum[c][topic] - nwc[c][V-1][topic] - newnwc[c][newV-1][topic])) *
                (nx[1][c][topic] + newnx[1][c][topic] + gamma1) / (nxsum[c][topic] + newnxsum[c][topic] + gamma0 + gamma1);


    u = ((double)random() / RAND_MAX) * (prob0+prob1);

    if (u > prob0) route = 1;
    else route = 0;

    if (route == 0) {
      newnw[_w][topic] += 1;
      newnwsum[topic] += 1;
      newnx[0][c][topic] += 1;
    }
    else {
      newnwc[c][_w][topic] += 1;
      newnwcsum[c][topic] += 1;
      newnx[1][c][topic] += 1;
    }
    newnd[m][topic] += 1;
    newndsum[m] += 1;
    newnxsum[c][topic] += 1;

    newzx[m][n] = route;

    return topic;
}

void model::compute_newtheta() {

  for (int m = 0; m < newM; m++) {
     int c = pnewdata->_docs[m]->cat;
     double constant = (hpu + newndsum[m] - newnd[m][K-1]) / (countU + newndsum[m] - newnd[m][K-1]);
     for (int k = 0; k < K; k++) {

       newtheta[m][k] = (newnd[m][k] + pu[c][k]) / (hpu + hpv + newndsum[m]) * constant;
     }

  }
}

void model::compute_newphi() {
    map<int, int>::iterator it;
    for (int k = 0; k < K; k++) {
	     for (int w = 0; w < newV; w++) {
	        it = pnewdata->_id2id.find(w);
	        if (it != pnewdata->_id2id.end()) {

              newphi[k][w] = (nw[it->second][k] + newnw[w][k] + pv[w]) / (countS + nwsum[k] + newnwsum[k] - nw[V-1][k] - newnw[newV-1][k]) * (hps + nwsum[k] + newnwsum[k] - nw[V-1][k] - newnw[newV-1][k]) / (hps + hpt + nwsum[k] + newnwsum[k]);
	        }

          for (int c = 0; c < C; c++) {

             newsigma[c][k][w] = (nwc[c][(it->second)][k] + newnwc[c][w][k] + pv[w]) / (countS + nwcsum[c][k] + newnwcsum[c][k] - nwc[c][V-1][k] - newnwc[c][newV-1][k]) * (hps + nwcsum[c][k] + newnwcsum[c][k] - nwc[c][V-1][k] - newnwc[c][newV-1][k]) / (hps + hpt + nwcsum[c][k] + newnwcsum[c][k]);

          }
	    }
    }
}

void model::compute_newpsi() {

    for (int k = 0; k < K; k++) {
	     for (int c = 0; c < C; c++) {
	        newpsi[k][c] = (nx[1][c][k] + newnx[1][c][k] + gamma1) / (nxsum[c][k] + newnxsum[c][k] + gamma0 + gamma1);
	     }
    }
}
