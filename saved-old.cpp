
// Just for the sample -- populate the intput data set
all_alloc_t build_input() {
   all_alloc_t vvi;

   for(int i = 0; i < NPH; i++) {
      alloc_t vi;
      for(int j = 0; j < M; j++) {
         vi.push_back(j+1);
      }
      vvi.insert(vi);
   }
   return vvi;
}

int main2() {
    all_alloc_t vvi = {};
    queue<int> vi;
    construct_alloc(vvi,0);
    //cout << vi << endl;
    cout << endl << vvi << endl;
    cout << "allocation-size " << vvi.size() << endl;
    return 0;
}

void epsilon_move_rule34(all_alloc_t &children, const alloc_t startpoint, all_alloc_t &valid_actions, bool cond) {
    int i, j;
    valid_actions.clear(); /* Alternative raise an assertion */

    for (i = 0; i < NPH; i++) {
        for (j = 0; j < NPH; j++) {
            if (i != j) {
                if (startpoint[i] > 1) {
                    /* Transfer a core from i to j */
                    alloc_t new_point = startpoint;
                    new_point[i]--;
                    new_point[j]++;
                    children.insert(new_point);

                    /* Enq/Store the list of valid actions */
                    if (cond) {
                        alloc_t action (NPH,0);
                        action[i] = -1;
                        action[j] = 1;
                        valid_actions.insert(action);
                    }
                }
            }
        }
    }
}

/* rule 1 expansion */
void expand_rule1(all_alloc_t &fbidden, set<alloc_t> &fset, set<alloc_t> &aset, int actv_id) {
    bool no_child = true;
    for (set<alloc_t>::iterator it = fset.begin(); it != fset.end(); ) {
        alloc_t tmp = *it;
        // cout << tmp << endl;
        fbidden.insert(tmp); // Insert all the elements of frontier set into fbidden set.

        /* Expand all the elements of the frontier */
        for (int idx = 0; idx < NPH; idx++) {
            alloc_t tmp2  = tmp;   
            if (tmp2[idx] > 1) {
                tmp2[idx]--;
                aset.insert(tmp2);
                no_child &= false;
            }
        }

        it = fset.erase(it);
    }
    //cout << "Forbidden Set = " << fbidden << endl;
    if (!no_child) {
        expand_rule1(fbidden,aset,fset,++actv_id);
    }
    // cout << "Exited : " <<actv_id << endl;
}

/* rule 2 expansion */
void expand_rule2(all_alloc_t &fbidden, set<alloc_t> &fset, set<alloc_t> &aset, int actv_id) {
    bool no_child = true;
    for (set<alloc_t>::iterator it = fset.begin(); it != fset.end(); ) {
        alloc_t tmp = *it;
        // cout << tmp << endl;
        fbidden.insert(tmp); // Insert all the elements of frontier set into fbidden set.

        /* Expand all the elements of the frontier */
        for (int idx = 0; idx < NPH; idx++) {
            alloc_t tmp2  = tmp;   
            if (tmp2[idx] < M) {
                tmp2[idx]++;
                aset.insert(tmp2);
                no_child &= false;
            }
        }

        it = fset.erase(it);
    }
    //cout << "Forbidden Set = " << fbidden << endl;
    if (!no_child) {
        expand_rule2(fbidden,aset,fset,++actv_id);
    }
    // cout << "Exited : " <<actv_id << endl;
}


void apply_action(alloc_t &point, const alloc_t action) {
    int i = 0;
    for (auto it = point.begin(); it != point.end(); it++,i++) {
        if ((*it + action[i] < NPH) && (*it + action[i] > 0))
            *it += action[i];
    }
}


// bool no_child = true;
    // for (set<alloc_t>::iterator it = fset.begin(); it != fset.end(); ) {
    //     alloc_t tmp = *it;
    //     // cout << tmp << endl;
    //     fbidden.insert(tmp); // Insert all the elements of frontier set into fbidden set.

    //     /* Expand all the elements of the frontier (Rule 1) */
    //     for (int idx = 0; idx < NPH; idx++) {
    //         alloc_t tmp2  = tmp;   
    //         if (tmp2[idx] > 1) {
    //             tmp2[idx]--;
    //             aset.insert(tmp2);
    //             no_child &= false;
    //         }
    //     }

    //     /* Expand all the elements of the frontier (Rule 3) */
    //     epsilon_move2_rule3(aset,init);

    //     it = fset.erase(it);
    // }
    // //cout << "Forbidden Set = " << fbidden << endl;
    // if (!no_child) {
    //     apply_action_rule13(fbidden,aset,fset,++actv_id);
    // }