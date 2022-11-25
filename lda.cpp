/*
 *
 * Copyright (c) 2009 Semantic Frontiers Group
 * University of Illinois at Urbana-Champaign
 * http://apfel.ai.uiuc.edu/
 *
 * Author: Michael Paul <mpaul39@gmail.com>
 *
 * ccLDA is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * ccLDA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/*
 * Based on the LDA Gibbs Sampler by Xuan-Hieu Phan
 *	http://gibbslda.sourceforge.net/
*/

#include "model.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void show_help();

int main(int argc, char ** argv) {
    model lda;

    if (lda.init(argc, argv)) {
	     show_help();
	     return 1;
    }

    if (lda.model_status == MODEL_STATUS_EST || lda.model_status == MODEL_STATUS_ESTC) {
	     // parameter estimation
	     lda.estimate();
       //lda.cuda_estimate();
    }

    if (lda.model_status == MODEL_STATUS_INF) {
	     // do inference
	     lda.inference();
    }

    return 0;
}

void show_help() {
    printf("Command line usage:\n");
    printf("\tlda -est -alpha <double> -beta <double> -ntopics <int> -niters <int> -savestep <int> -twords <int> -dfile <string>\n");
    printf("\tlda -estc -dir <string> -model <string> -niters <int> -savestep <int> -twords <int>\n");
    printf("\tlda -inf -dir <string> -model <string> -niters <int> -twords <int> -dfile <string>\n");
    // printf("\tlda -inf -dir <string> -model <string> -niters <int> -twords <int> -dfile <string> -withrawdata\n");
}
