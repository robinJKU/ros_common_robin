#ifndef BS_COEFFS_H
#define BS_COEFFS_H

#include <iostream>
#include <stdio.h>
#include <string.h> /* strlen */
#include <cmath>
#include <vector>
#include <exception>
#include <algorithm> // sort()
#include "Eigen/Dense" // MatrixXd
#include "Eigen/Core" // VectorXf

#define NUMEPS 1E-8 //!< approx. zero
#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names)) //!< number of element in field_names array of strings

using namespace std;


/**
 * B-spline coefficient class.
 */
class Bs_coeffs {
    std::vector<Eigen::MatrixXd> coeff_mat_vec; //!<  vector of coefficient matrices (for each non-zero interval)
    unsigned int degree;    //!<  maximum degree of curve
    Eigen::RowVectorXd knots;  //!<  knots vector
    Eigen::RowVectorXd get_coeff_vec(unsigned int j, unsigned int d, unsigned int interval);
    Eigen::RowVectorXd get_unique_knots();
    void generate_coeff_mat_vec();
    Eigen::RowVectorXd poly_mult(Eigen::RowVectorXd& poly1, Eigen::RowVectorXd& poly2);
public:
    Bs_coeffs();
    Bs_coeffs(unsigned int degree, Eigen::RowVectorXd& knots);
    Bs_coeffs(unsigned int degree, Eigen::RowVectorXd& knots, std::vector<Eigen::MatrixXd>& coeff_mat_vec);
    unsigned int get_interval(double parameter);
    unsigned int get_unique_interval(double parameter);
    Eigen::MatrixXd* get_coeff_mat(double parameter);
    unsigned int get_degree();
    std::vector<Eigen::MatrixXd> get_coeff_mat_vec();
    Eigen::RowVectorXd get_knots();
};


#endif
