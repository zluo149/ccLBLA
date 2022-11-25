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


#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "strtokenizer.h"
#include "dataset.h"

using namespace std;

int dataset::write_wordmap(string wordmapfile, mapword2id * pword2id) {
    FILE * fout = fopen(wordmapfile.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to write!\n", wordmapfile.c_str());
	return 1;
    }

    mapword2id::iterator it;
    fprintf(fout, "%d\n", pword2id->size());
    for (it = pword2id->begin(); it != pword2id->end(); it++) {
	fprintf(fout, "%s %d\n", (it->first).c_str(), it->second);
    }

    fclose(fout);

    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapword2id * pword2id) {
    pword2id->clear();

    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	return 1;
    }

    char buff[BUFF_SIZE_SHORT];
    string line;

    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);

    for (int i = 0; i < nwords; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;

	strtokenizer strtok(line, " \t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}

	pword2id->insert(pair<string, int>(strtok.token(0), atoi(strtok.token(1).c_str())));
    }

    fclose(fin);

    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapid2word * pid2word) {
    pid2word->clear();

    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	return 1;
    }

    char buff[BUFF_SIZE_SHORT];
    string line;

    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);

    for (int i = 0; i < nwords; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;

	strtokenizer strtok(line, " \t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}

	pid2word->insert(pair<int, string>(atoi(strtok.token(1).c_str()), strtok.token(0)));
    }

    fclose(fin);

    return 0;
}

int dataset::read_trndata(string dfile, string wordmapfile) {
    mapword2id word2id;

    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", dfile.c_str());
	return 1;
    }

    mapword2id::iterator it;
    char buff[BUFF_SIZE_LONG];
    string line;

    // get the number of documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }

    // allocate memory for corpus
    if (docs) {
	deallocate();
    } else {
	docs = new document*[M];
    }

    // set number of words to zero
    V = 0;
    C = 0;

    for (int i = 0; i < M; i++) {
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();

	if (length <= 0) {
	    printf("Invalid (empty) document!\n");
	    deallocate();
	    M = V = 0;
	    return 1;
	}

	// get doc collection

	string c = strtok.token(0);
	int cat = atoi(c.c_str());

	if (cat > C) C = cat;

	// allocate new document
	document * pdoc = new document(length - 1, cat);

	for (int j = 1; j < length; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., new word
		pdoc->words[j - 1] = word2id.size();
		word2id.insert(pair<string, int>(strtok.token(j), word2id.size()));
	    } else {
		pdoc->words[j - 1] = it->second;
	    }
	}

	// add new doc to the corpus
	add_doc(pdoc, i);
    }

    fclose(fin);

    // write word map to file
    if (write_wordmap(wordmapfile, &word2id)) {
	return 1;
    }

    // update number of words
    V = word2id.size();

    C = C+1;

    return 0;
}

int dataset::read_newdata(string dfile, string wordmapfile) {
    mapword2id word2id;
    map<int, int> id2_id;

    read_wordmap(wordmapfile, &word2id);
    if (word2id.size() <= 0) {
	     printf("No word map available!\n");
	     return 1;
    }

    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	     printf("Cannot open file %s to read!\n", dfile.c_str());
	     return 1;
    }

    mapword2id::iterator it;
    map<int, int>::iterator _it;
    char buff[BUFF_SIZE_LONG];
    string line;

    // get number of new documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	     printf("No document available!\n");
	     return 1;
    }

    // allocate memory for corpus
    if (docs) {
	     deallocate();
    } else {
	     docs = new document*[M];
    }

    _docs = new document*[M];

    // set number of words to zero
    V = 0;
    C = 0;

    for (int i = 0; i < M; i++) {
	     fgets(buff, BUFF_SIZE_LONG - 1, fin);
	     line = buff;
	     strtokenizer strtok(line, " \t\r\n");
	     int length = strtok.count_tokens();

       // get doc collection

       string c = strtok.token(0);
       int cat = atoi(c.c_str());

       if (cat > C) C = cat;

	vector<int> doc;
	vector<int> _doc;
	for (int j = 1; j < length; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., word unseen in training data
		// do anything? (future decision)
	    } else {
		int _id;
		_it = id2_id.find(it->second);
		if (_it == id2_id.end()) {
		    _id = id2_id.size();
		    id2_id.insert(pair<int, int>(it->second, _id));
		    _id2id.insert(pair<int, int>(_id, it->second));
		} else {
		    _id = _it->second;
		}

		doc.push_back(it->second);
		_doc.push_back(_id);
	    }
	}

	// allocate memory for new doc
	document * pdoc = new document(doc, cat);
	document * _pdoc = new document(_doc, cat);

	// add new doc
	add_doc(pdoc, i);
	_add_doc(_pdoc, i);
    }

    fclose(fin);

    // update number of new words
    V = id2_id.size();

    C = C+1;

    return 0;
}

int dataset::read_newdata_withrawstrs(string dfile, string wordmapfile) {
    mapword2id word2id;
    map<int, int> id2_id;

    read_wordmap(wordmapfile, &word2id);
    if (word2id.size() <= 0) {
	printf("No word map available!\n");
	return 1;
    }

    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", dfile.c_str());
	return 1;
    }

    mapword2id::iterator it;
    map<int, int>::iterator _it;
    char buff[BUFF_SIZE_LONG];
    string line;

    // get number of new documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }

    // allocate memory for corpus
    if (docs) {
	deallocate();
    } else {
	docs = new document*[M];
    }
    _docs = new document*[M];

    // set number of words to zero
    V = 0;

    for (int i = 0; i < M; i++) {
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();

	// get collection

	string c = strtok.token(0);
	char c0[c.length()];
	for (int ci = 0; ci < c.length(); ci++)
	  c0[ci] = c[ci];
	int cat = atoi(c0);

	vector<int> doc;
	vector<int> _doc;
	for (int j = 1; j < length - 1; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., word unseen in training data
		// do anything? (future decision)
	    } else {
		int _id;
		_it = id2_id.find(it->second);
		if (_it == id2_id.end()) {
		    _id = id2_id.size();
		    id2_id.insert(pair<int, int>(it->second, _id));
		    _id2id.insert(pair<int, int>(_id, it->second));
		} else {
		    _id = _it->second;
		}

		doc.push_back(it->second);
		_doc.push_back(_id);
	    }
	}

	// allocate memory for new doc
	document * pdoc = new document(doc, line, cat);
	document * _pdoc = new document(_doc, line, cat);

	// add new doc
	add_doc(pdoc, i);
	_add_doc(_pdoc, i);
    }

    fclose(fin);

    // update number of new words
    V = id2_id.size();

    return 0;
}
