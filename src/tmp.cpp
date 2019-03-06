  // ET  
    // switch (bench) {
    //     case BENCH_LACE_CILKSORT : {
    //         y = 1/(0.004254273510905619*x - 0.0008304322052644017);
    //         break;
    //     }
    //     case BENCH_LACE_DFS : {
    //         y = 1/(0.0052912772458832795*x + 0.002118298041182201);
    //         break;
    //     }
    //     case BENCH_LACE_PI : {
    //         y = 1/(0.003688517486111831*x + 0.0011143885733635485);
    //         break;
    //     }
    //     case BENCH_LACE_QUEENS : {
    //         y = 1/(0.003959440926636103*x + 0.001000559321893281);
    //         break;
    //     }
    //     case BENCH_LACE_FIB : {
    //         y = 1/(0.005266330309541947*x + 0.0021923143697295602);
    //         break;
    //     }
    //     case BENCH_PRSC_BLACKSCHOLES : {
    //         y = 1/(0.009683975849832148*log(x) - 0.002893833328927599);
    //         break;
    //     }

    //     case BENCH_PRSC_STREAMCLUSTER:{
    //         y = 1/(0.0024794976904501873*log(x) - 0.0007496878381248601);
    //         break;
    //     }

    //     case BENCH_PRSC_CANNEAL : {
    //         y = 1/(2.153404950919569e-05*log(x) + 0.0001677868522003634);
    //         break;
    //     }

    //     // case BENCH_PRSC_SWAPTIONS : {

    //     // }

    //     case BENCH_PRSC_FLUIDANIMATE : {
    //         y = 1/(0.0021638291897166048*log(x) - 0.001202365402456029);
    //         break;
    //     }

    //     case BENCH_PRSC_BODYTRACK : {
    //         y = 1/(0.003533669349872262*log(x) - 0.0027172901605883966);
    //         break;
    //     }

    //     case BENCH_PRSC_DEDUP : {
    //         y = 1/(0.0002555000090127519*log(x) - 7.198705194815577e-05);
    //         break;
    //     }
    //     default : {
    //         cerr << "(ET) Unrecognized Benchmark : " << bench << endl;
    //         exit(EXIT_FAILURE);
    //     }
    // }

// void construct_alloc(all_alloc2_t &vvi, \
//                      const vector<int> &bench, \
//                      double deadline, \
//                      int ph,\
//                      alloc2_t &opt_point,\
//                      double &opt_pkp_power,\
//                      double &opt_exec_time) {
//     alloc2_t vi2;
//     static long unsigned int cnt = 0;
//     double et, pkp;
//     long unsigned int r   = M;
//     long unsigned int dec;
//     long unsigned int j   = 0;

//     if (ph == NPH) {
//         // cout << "deadline (calloc) : " << deadline << endl;
//         //cout << "ph : " << ph << ", invoc : " << cnt << endl;
//         dec = cnt++;

        
//         /* Resolve cnt into M-radix number */
//         for (j = 0; j < NPH; j++) {
//             //cout << dec%r + 1 << ",";
//             phase_t phinfo;
//             phinfo.alloc = dec%r + 1;
//             phinfo.bench_id = bench[j];

//             // if (cnt == 800083) {
//             //     cout << "phinfo:"<<phinfo<<endl;
//             // }
//             // cout << "Bench Created : " << phinfo.bench_id << endl;

//             vi2.push_back(phinfo);
//             dec = dec/r;
//         }
//         // cout << "cnt="<<cnt<<","<<vi2<<endl;
//         /* Oracle Evaluation */
//         et = compute_execution_time(vi2);
//         // cout << "construct alloc et " << et << "\n";
//         if (et <= deadline && et >= 0) {
//             pkp = compute_pkpower(vi2);
//             if (pkp <= opt_pkp_power) {
//                 opt_point     = vi2;
//                 opt_pkp_power = pkp;
//                 opt_exec_time = et;
//             }
//             // cout << "construct alloc pkp" << pkp << "\n";
//         }

//         /* Add an extreme point if no feasible point is found */
//         // if (opt_point.size() < NPH) {
//         //     cout << "NO FEASIBLE POINT FOUND" << endl;
//         //     for(unsigned int i = 0; i < NPH; i++) {
//         //         phase_t phinfo;
//         //         phinfo.alloc = ULIM;
//         //         phinfo.bench_id = bench[i];
//         //         opt_point.push_back(phinfo);
//         //     }
//         // }
//         if (et >= 0) {
//             vvi.insert(vi2);
//         }
//         //cout << "}\n";
//         return;
//     }
//     for (int m = 0; m < ULIM; m++) {
//         construct_alloc(vvi,bench,deadline,ph+1,opt_point,opt_pkp_power,opt_exec_time);
//     }
// }