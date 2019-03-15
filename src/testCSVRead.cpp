#include <iostream>
#include <vector>
#include <stdexcept>
#include <rapidcsv.h>

int main() {
  rapidcsv::Document doc("/home/amaity/Dropbox/workspace/islped_bonmip/workloads-exp1/wkld_4.csv",rapidcsv::LabelParams(-1,-1));
  try {
    std::vector<double> deadlines = doc.GetRow<double>(100);
    for (int i = 0; i < deadlines.size(); i++) {
        std::cout << i+1 <<","<< deadlines[i] << std::endl;
    }
  } catch(std::out_of_range) {
    std::cout << "Reached End of CSV File" << std::endl;
  }
}
