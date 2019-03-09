#include "ptss_bonmin.hpp"
#include "ptss_config.hpp"

bool MyTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,\
                       Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style) {
  n = NPH+1; //number of variable
  m = 3;     //number of constraints
  nnz_jac_g = 7;//number of non zeroes in Jacobian
  nnz_h_lag = 2;//number of non zeroes in Hessian of Lagrangean
  index_style = TNLP::FORTRAN_STYLE;
  return true;
}