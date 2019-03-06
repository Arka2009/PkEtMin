#include <iomanip>
#include <iostream>
#include <vector>
#include <vector>
#include <nlopt.hpp>
#include <gurobi_c++.h>

using namespace std;

int main(int argc, char *argv[]) {
  try {
    GRBEnv env = GRBEnv();

    GRBModel model = GRBModel(env);

    double a = 1e-3;
    double b = 1e-3;
    double a2 = 3e-3;
    double b2 = 2e-3;
    double D  = 450.0;

    GRBVar x = model.addVar(1.0, 16.0, 0.0, GRB_INTEGER, "x");
    GRBVar y = model.addVar(1.0, 16.0, 0.0, GRB_INTEGER, "y");
    GRBVar t = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t");

    // Set objective: maximize x + y + 2 z

    model.setObjective(1*t, GRB_MINIMIZE);

    // Add constraint: x + 2 y + 3 z <= 4

    model.addConstr(x - t <= 0, "c0");

    // Add constraint: x + y >= 1

    model.addConstr(y - t <= 0, "c1");

    // Add constraint

    model.addQConstr( (1/(a*x)+b)) + (1/(a2*y+b2)) - D <= 0, "c2");
    // model.addQConstr( x*x + y*x + t*t - D <= 0, "c2");

    // Optimize model

    model.optimize();

    cout << x.get(GRB_StringAttr_VarName) << " "
         << x.get(GRB_DoubleAttr_X) << endl;
    cout << y.get(GRB_StringAttr_VarName) << " "
         << y.get(GRB_DoubleAttr_X) << endl;
    cout << t.get(GRB_StringAttr_VarName) << " "
         << t.get(GRB_DoubleAttr_X) << endl;

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}

