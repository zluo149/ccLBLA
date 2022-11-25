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

#ifndef _UTILS_H
#define _UTILS_H

#include <string>

using namespace std;

class model;

class utils {
public:
    // parse command line arguments
    static int parse_args(int argc, char ** argv, model * pmodel);
    
    // read and parse model parameters from <model_name>.others
    static int read_and_parse(string filename, model * model); 
    
    // generate the model name for the current iteration
    // iter = -1 => final model
    static string generate_model_name(int iter);  
    
    // sort    
    static void sort(vector<double> & probs, vector<int> & words);
    static void quicksort(vector<pair<int, double> > & vect, int left, int right);
};

#endif

