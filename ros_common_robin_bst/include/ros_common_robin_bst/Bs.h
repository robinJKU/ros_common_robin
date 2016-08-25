/*
 * B-spline toolbox
 * Alexander Reiter, Dominik Kaserer
 * Institute of Robotics, Johannes Kepler University Linz
 * November 2014 - May 2016
 */
// $Id$
/**
 * @file bst.cpp
 * Implementation file of the B-spline toolbox.
 * TODO: 
 * - evaluate use of sparse matrices, #include <Eigen/Sparse>, Eigen::SparseMatrix<double>
 * - only store non-zero parts of coefficient matrices (for less zero-multiplications and lower memory usage)
 * - more and more appropriate exceptions
 *
 * @brief Can use "brief" tag to explicitly generate comments for file documentation.
 *
 * @author Alexander Reiter
 * @author Dominik Kaserer
 * @version 16.05
 */
// $Log$

#ifndef BS_H
#define BS_H

#include <iostream>
#include <stdio.h>
#include <string.h> /* strlen */
#include <cmath>
#include <vector>
#include <exception>
#include <algorithm> // sort()
#include "Eigen/Dense" // MatrixXd
#include "Eigen/Core" // VectorXf

#include "Bs_coeffs.h"

#define NUMEPS 1E-8 //!< approx. zero
#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names)) //!< number of element in field_names array of strings

using namespace std;

/**
 * B-spline class.
 */
class Bs {
    Bs_coeffs coeffs;
    Eigen::RowVectorXd ctrl_pts;    //!<  control points vector
    double par_start;   //!<  start of parameter interval
    double par_end; //!<  end of parameter interval
    Eigen::RowVectorXd linspace(double a, double b, unsigned int n);
    Eigen::RowVectorXd get_generic_knots(unsigned int degree, unsigned int size_ctrl_pts);
    double normalize(double par);
    void normalize_row_vec(Eigen::RowVectorXd& vec);
    Eigen::VectorXd get_pow_vec(double par, int highest_pow);
    Eigen::VectorXd get_pow_vec(double par);
    Eigen::VectorXd get_exp_vec(int highest_exp);
public:
    Bs();
    Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts);
    Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts, double par_start, double par_end);
    Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts, Eigen::RowVectorXd& knots);
    Bs(unsigned int degree, Eigen::MatrixXd& approx_data, unsigned int n_ctrl_pts);
    Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts, double par_start, double par_end, Bs_coeffs coeffs);
    void set_par_start(double par_start);
    void set_par_end(double par_end);
    double get_value(double parameter);
    double get_value(double parameter, unsigned int derivative);
    Eigen::RowVectorXd get_values(Eigen::RowVectorXd& parameters);
    Eigen::RowVectorXd get_values(Eigen::RowVectorXd& parameters, Eigen::VectorXd& derivatives);
    Eigen::RowVectorXd to_unique(Eigen::RowVectorXd& vec);
    double get_par_start();
    double get_par_end();
    Eigen::RowVectorXd get_ctrl_pts();
    Bs_coeffs* get_coeffs();
};

#endif
