#include <cmath>
#include <vector>
#include <random>
#include "ptss_StatisticalDistribution.hpp"
extern "C" {
#include <ecotools/cpu_uarch.h>
}

#define SAMPLE_SIZE 5000

int main() {
    vector<PTSSETDistType*> ets;
    std::ofstream ofs;
    ofs.open("DistTest.csv",ios::out | ios::trunc);

    for (int i = 0; i < 8; i++) {
        // cout << "Generating Dist "<< i << endl;
        phase_t ph(0,i+1);
        PTSSETDistType *gen = new PTSSETDistType(ph);
        if (gen != NULL) {
            ets.push_back(gen);
        }
        else {
            // cout << "Cannot create new distribution" << endl;
            PRINTERROR("Cannot create new distribution");
        }
        
        if (i < 7)
            ofs << "Alloc-"<<i<<",";
        else
            ofs << "Alloc-"<<i<<std::endl;
    }

    
    for (int j = 0; j < SAMPLE_SIZE; j++) {
        
        for (int i = 0; i < 8; i++) {
            if (i < 7)
                ofs << ets[i]->sample() << ",";
            else
                ofs << ets[i]->sample() << endl;
        }
    }

    for (int i = 0; i < 8; i++) {
        delete ets[i];
    }
    ofs.close();
    return 0;
}