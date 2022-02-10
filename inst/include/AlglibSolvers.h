// AlglibSolvers.h: biosensors.usc glue
//
// Copyright (C) 2019 - 2021  Juan C. Vidal and Marcos Matabuena
//
// This file is part of biosensors.usc.
//
// biosensors.usc is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// biosensors.usc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with biosensors.usc.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _ALGLIB_SOLVERS_h
#define _ALGLIB_SOLVERS_h

#include <stdlib.h>
#include <string>
#include <RcppArmadillo.h>
#include "stdafx.h"
#include "optimization.h"
#include "solvers.h"


namespace bio {

  enum QP_Solver {Quick, BLEIC, DENSE_AUL};

  std::string mat2string(arma::mat A) {
    std::string str("");
    if (A.is_empty()) {
      str += "[]";
    } else {
      str += "[";
      for (arma::uword i=0; i < A.n_rows; i++) {
        str += "[";
        for (arma::uword j=0; j < A.n_cols; j++) {
          str += std::to_string(A(i,j));
          if (j < A.n_cols-1)
            str += ",";
        }
        str += "]";
        if (i < A.n_rows-1)
          str += ",";
      }
      str += "]";
    }
    return str;
  }

  std::string vec2string(arma::vec x) {
    std::string str("");
    if (x.is_empty()) {
      str += "[]";
    } else {
      str += "[";
      for (arma::uword i=0; i < x.n_elem; i++) {
        str += std::to_string(x(i));
        if (i < x.n_elem-1)
          str += ",";
      }
      str += "]";
    }
    return str;
  }



  std::string ct2string(arma::uword n) {
    std::string str("[");
    for (arma::uword i=0; i < n; i++) {
        str += "-1";
        if (i < n-1)
          str += ",";
    }
    str += "]";
    return str;
  }

  arma::vec real_1d_array2vec(alglib::real_1d_array x) {
    arma::vec v(x.length());
    // double* content = x.getcontent();
    for (int i=0; i < x.length(); i++) {
      v(i) = x.operator()(i);
    }
    return v;
  }

  arma::mat real_2d_array2mat(alglib::real_2d_array x) {
    arma::mat a(x.rows(),x.cols());
    for (int i=0; i < x.rows(); i++) {
      for (int j=0; j < x.cols(); j++) {
        a(i,j) = x.operator()(i,j);
      }
    }
    return a;
  }

  arma::vec quadprog(arma::mat A, arma::vec b, arma::mat C, arma::vec d,
                     arma::vec lb = arma::vec(), arma::vec ub = arma::vec(), arma::vec x0 = arma::vec()) {
    // arma::vec quadprog(QP_Solver type, arma::mat A, arma::vec b, arma::mat C, arma::vec d,
    //                    arma::vec lb = arma::vec(), arma::vec ub = arma::vec(), arma::vec x0 = arma::vec()) {
    alglib::minqpstate state;
    alglib::minqpreport rep;
    alglib::real_1d_array x;
    // Create solver, set quadratic/linear terms
    alglib::minqpcreate(A.n_cols, state);

    // convert armadillo to alglib structures
    alglib::real_2d_array as = mat2string(A).c_str();
    alglib::minqpsetquadraticterm(state, as);

    alglib::real_1d_array bs = vec2string(b).c_str();
    alglib::minqpsetlinearterm(state, bs);

    alglib::real_2d_array cs = mat2string(arma::join_horiz(C,d)).c_str();
    alglib::integer_1d_array cts = ct2string(C.n_rows).c_str();
    alglib::minqpsetlc(state, cs, cts);

    if (!lb.is_empty() && !ub.is_empty()) {
      alglib::real_1d_array lbs = vec2string(lb).c_str();
      alglib::real_1d_array ubs = vec2string(ub).c_str();
      alglib::minqpsetbc(state, lbs, ubs);
    }

    if (!x0.is_empty()) {
      alglib::real_1d_array x0s = vec2string(x0).c_str();
      alglib::minqpsetstartingpoint(state, x0s);
    }

    // Set scale of the parameters.
    // It is strongly recommended that you set scale of your variables.
    // Knowing their scales is essential for evaluation of stopping criteria
    // and for preconditioning of the algorithm steps.
    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
    //
    // NOTE: for convex problems you may try using minqpsetscaleautodiag()
    //       which automatically determines variable scales.
    arma::vec scale(A.n_cols);
    scale.fill(10);
    // alglib::real_1d_array s = "[1,1]";
    alglib::real_1d_array s = vec2string(scale).c_str();
    alglib::minqpsetscale(state, s);
    // alglib::minqpsetscaleautodiag(state);

    // Solve problem with DENSE-AUL solver.
    // alglib::minqpsetalgodenseaul(state, 1.0e-6, 1.0e+3, 5);
    // alglib::minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
    alglib::minqpsetalgodenseipm(state, 1.0e-5);

    //
    // switch(type) {
    // case Quick:
    //   // Solve problem with QuickQP solver.
    //   minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
    //   break;
    // case BLEIC:
    //   // Solve problem with BLEIC-based QP solver.
    //    alglib::minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
    //   break;
    // case DENSE_AUL:
    //   // Solve problem with DENSE-AUL solver.
    //   alglib::minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5);
    //   break;
    // default:
    //   // Solve problem with DENSE-AUL solver.
    //   alglib::minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5);
    //   break;
    // }

    alglib::minqpoptimize(state);

    alglib::minqpresults(state, x, rep);
    // printf("POST OPT: %d\n", int(rep.terminationtype));

    return real_1d_array2vec(x);
  }

  arma::mat linear_solver(arma::mat A, arma::mat B) {
    alglib::real_2d_array As = mat2string(A).c_str();
    alglib::real_2d_array Bs = mat2string(B).c_str();
    // alglib::real_2d_array As = mat2string(A.t() * A).c_str();
    // alglib::real_2d_array Bs = mat2string(A.t() * B).c_str();
    alglib::real_2d_array Xs;
    alglib::densesolverreport rep;
    alglib::ae_int_t info;


    // alglib::real_2d_array LUAs = mat2string(A).c_str();
    // alglib::integer_1d_array pivots;
    // alglib::rmatrixlu(LUAs, A.n_rows, A.n_cols, pivots);
    // alglib::rmatrixmixedsolvem (As, LUAs, pivots, A.n_rows, Bs, B.n_cols, info, rep, Xs);
    // // rmatrixlusolvem (As, LUAs, pivots, A.n_rows, Bs, B.n_cols, info, rep, Xs);
    // printf("LU Solution \n");
    // for (arma::uword i=0; i < Xs.rows(); i++) {
    //   for (arma::uword j=0; j < Xs.cols(); j++) {
    //     printf("%f ", Xs.operator()(i,j));
    //   }
    //   printf("\n");
    // }
    alglib::rmatrixsolvem (As, A.n_rows, Bs, B.n_cols, true, info, rep, Xs);
    // printf("\nSolution \n");
    // for (arma::uword i=0; i < Xs.rows(); i++) {
    //   for (arma::uword j=0; j < Xs.cols(); j++) {
    //     printf("%f ", Xs.operator()(i,j));
    //   }
    //   printf("\n");
    // }

    return real_2d_array2mat(Xs);
  }

}


#endif
