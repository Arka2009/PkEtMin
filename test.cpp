#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  int A[] = {0, 2, 3, 1, 4,67,9,8,67,2,3,5,5};
  //vector<int> A = {0, 2, 3, 1, 4,67,9,8,67,2,3,5,5};
  const int N = sizeof(A) / sizeof(int);

  cout << "Index of max element: "
       << distance(A, max_element(A, A + N))
       << endl;

  return 0;
}
