#include "Bs.h"


/**
 * Generates linearly spaced vectors.
 * Equivalent to MATLAB's linspace function.
 *
 * @param[in] a starting value
 * @param[in] b end value
 * @param[in] n number of values
 * @return resulting row vector
 */
Eigen::RowVectorXd Bs::linspace(double a, double b, unsigned int n) {
    Eigen::RowVectorXd vec = Eigen::RowVectorXd::Zero(n);
    for(size_t i = 0; i < n; i++) {
        vec(i) = a + (b - a) * i/(n - 1);
    }
    return vec;
}

/**
 * Computes generic knot distribution such that
 * knots = [(degree zeros), linspace(0,1,m-2*(degree + 1)), (degree ones)].
 *
 * @param[in] degree maximum degree of curve
 * @param[in] size_ctrl_pts number of control points
 * @return vector of knots
 */
Eigen::RowVectorXd Bs::get_generic_knots(unsigned int degree, unsigned int size_ctrl_pts) {
    unsigned int m = size_ctrl_pts + degree + 1; // number of knots
    Eigen::RowVectorXd knots = Eigen::RowVectorXd::Ones(m);
    knots << Eigen::RowVectorXd::Zero(degree), linspace(0., 1., m - 2*degree);
//     knots << linspace(0., 1., m);
    return knots;
}

/**
 * (Empty) default constructor.
 */
Bs::Bs() {
    return;
}

/**
 * Constructor.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] ctrl_pts row vector of control points
 */
Bs::Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts) {
    this->ctrl_pts = ctrl_pts;
    this->par_start = 0.;
    this->par_end = 1.;
    Eigen::RowVectorXd knots = get_generic_knots(degree, (unsigned int) ctrl_pts.size());
    this->coeffs = Bs_coeffs(degree, knots);
}

/**
 * Constructor.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] ctrl_pts row vector of control points
 * @param[in] par_start lower bound of parameters
 * @param[in] par_end upper bound of parameters
 */
Bs::Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts, double par_start, double par_end) {
    this->ctrl_pts = ctrl_pts;
    this->par_start = par_start;
    this->par_end = par_end;
    Eigen::RowVectorXd knots = get_generic_knots(degree, (unsigned int) ctrl_pts.size());
    this->coeffs = Bs_coeffs(degree, knots);
}

/**
 * Constructor.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] ctrl_pts row vector of control points
 * @param[in] knots row vector of knots
 */
Bs::Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts, Eigen::RowVectorXd& knots) {
    if(ctrl_pts.size() != knots.size() - degree - 1) {
        //throw x_wrong_range;
    }
    this->ctrl_pts = ctrl_pts;
    this->par_start = knots(0);
    this->par_end = knots.tail(1)(0); // first element of tail block(= last element)
    normalize_row_vec(knots);
    this->coeffs = Bs_coeffs(degree, knots);
}

/**
 * Approximation constructor.
 *
 * Approximation data matrix column description:
 * | tk | q^(der) | der | 
 * where
 * tk: knot
 * q^(der): der-th derivative of q at t_k
 * der: order of derivative
 *
 * Approximation is done by computing the least-squares solution
 * (using QR decomposition) of the linear system:
 * | R_0_n(t0) ... R_e_n(t0) | | d0 |   |   q0   |
 * | R_0_n(t1) ... R_e_n(t1) | | d1 | = |   q1   |
 * | ...                     | | .. |   |   ..   |
 * | R_0_n(te) ... R_e_n(te) | | de |   |   qe   |
 * | R^(1)_x_y(tz) ...       |          | qx^(1) |    <-- first order derivative
 * | R^(2)_a_b(tc) ...       |          | qa^(2) |    <-- second order derivative
 * resulting in the control points of the corresponding spline.
 * The solution defaults to interpolation in case of a full-rank matrix.
 * See http://eigen.tuxfamily.org/dox-devel/group__LeastSquares.html for more information.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] approx_data data used for appoximation (see description above)
 * @param[in] n_ctrl_pts   number of control points
 */
Bs::Bs(unsigned int degree, Eigen::MatrixXd& approx_data, unsigned int n_ctrl_pts) {
    // normalize parameters
    this->par_start = approx_data.col(0).minCoeff();
    this->par_end = approx_data.col(0).maxCoeff();
    
    // generate knots vector
    Eigen::RowVectorXd knots = get_generic_knots(degree, n_ctrl_pts);
    Eigen::RowVectorXd val_vec = Eigen::VectorXd();
    Eigen::RowVectorXd par_vec = Eigen::VectorXd();
    for(int i = 0;  i < approx_data.rows(); i++) {
        if((int) approx_data(i,2) == 0) { 
            val_vec.conservativeResize(1,val_vec.cols() + 1);
            par_vec.conservativeResize(1,par_vec.cols() + 1);
            val_vec(par_vec.cols() - 1) = approx_data(i,1);
            par_vec(par_vec.cols() - 1) = approx_data(i,0);
        }
    }
    /*
    double val_span = val_vec.maxCoeff() -  val_vec.minCoeff();
    Eigen::RowVectorXd a_par_vec = Eigen::VectorXd::Ones(approx_data.rows());
    par_vec(0) = 0.;
    for(int i = 1;  i < approx_data.rows(); i++) {
        a_par_vec(i) = a_par_vec(i-1) + 4./val_span;
    }
    Eigen::RowVectorXd knots = Eigen::VectorXd::Ones(n_ctrl_pts + degree + 2);
    double span = (approx_data.rows() + 1)/((double) n_ctrl_pts - degree + 1);
    knots << Eigen::VectorXd::Zero(n_ctrl_pts + 1);
    for(unsigned int j = 1; j <= n_ctrl_pts - degree; j++) {
        int i = (int) j*(int) span;
        double alpha = (double) j*span - i;
        knots(degree + j) = (1.-alpha)*a_par_vec(i - 1) + alpha*a_par_vec(i);
    }
    */
    //cout << knots << endl << endl;
    // generate coefficients
    this->coeffs = Bs_coeffs(degree, knots);
    
    // assemble approximation matrix
    Eigen::MatrixXd mat(approx_data.rows(), n_ctrl_pts);
    for(int i = 0; i < approx_data.rows(); i++) {
		//cout << approx_data(i,0) << " (" << this->par_start << ", " << this->par_end << "), ";
        double parameter = normalize(approx_data(i,0));
        Eigen::VectorXd exp_vec = Eigen::VectorXd::Ones(degree + 1); // vector of exponents for derivatives
    
        for(unsigned int j = degree; j > degree - (int) approx_data(i,2); j--) {
            Eigen::VectorXd this_exp_vec = Bs::get_exp_vec(j);
            exp_vec = exp_vec.cwiseProduct(this_exp_vec);
        }
        Eigen::VectorXd pow_vec = get_pow_vec(parameter, degree - (int) approx_data(i,2));
		//cout << parameter << ", " << approx_data(i,2) << "\t" << pow_vec.transpose() << endl << endl;
        mat.row(i) = 1./pow((double) (this->par_end - this->par_start), (int) approx_data(i,2))*((*(coeffs.get_coeff_mat(parameter)))*exp_vec.cwiseProduct(pow_vec)).transpose();
    }
    Eigen::VectorXd rhs = approx_data.col(1) - mat.col(0)*val_vec(0) - mat.rightCols<1>()*val_vec.rightCols<1>();
    mat = mat.block(0,1,mat.rows(),mat.cols()-2);
    // QR-based pseudo-inverse solution of mat*ctrl_pts = vals
    this->ctrl_pts = Eigen::VectorXd::Zero(n_ctrl_pts);
    this->ctrl_pts(0) = val_vec(0);
    this->ctrl_pts.block(0,1,n_ctrl_pts-2,1) = mat.fullPivHouseholderQr().solve(rhs).transpose();
    this->ctrl_pts.tail(1) = val_vec.rightCols<1>();
}

/**
 * Warm start constructor.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] ctrl_pts row vector of control points
 * @param[in] par_start lower bound of parameters
 * @param[in] par_end upper bound of parameters
 * @param[in] coeffs coefficient data
 */
Bs::Bs(unsigned int degree, Eigen::RowVectorXd& ctrl_pts, double par_start, double par_end, Bs_coeffs coeffs) {
    this->ctrl_pts = ctrl_pts;
    this->par_start = par_start;
    this->par_end = par_end;
    this->coeffs = coeffs;
}

/**
 * Reduces the entries of vec to unique values.
 * Expects a vector with entries sorted in ascending order.
 *
 * @param[in] vec vector to be reduced
 * @return reduced vector
 */
Eigen::RowVectorXd Bs::to_unique(Eigen::RowVectorXd& vec) {
    Eigen::RowVectorXd red(vec.size());
    red(0) = vec(0);
    int j = 1;
    for(int i = 1; i < vec.size(); i++) {
        if(fabs(vec(i) - vec(i-1)) > NUMEPS) {
            red(j) = vec(i);
            j++;
        }
    }
    red.conservativeResize(1, j);
    return red;
}

/**
 * Normalizes par in the interval [par_start, par_end]
 *
 * @param[in] par parameter to be normalized
 * @return vector of parameter powers
 */
double Bs::normalize(double par) {
    // robust calculation
    return min(max(0., (par - this->par_start)/(this->par_end - this->par_start) ), 1.);
}

/**
 * Normalized all elements of vec in the interval [par_start, par_end]
 *
 * @param[in] vec parameter vector to be normalized
 */
void Bs::normalize_row_vec(Eigen::RowVectorXd& vec) {
    for(int i = 0; i < vec.size(); i++) {
        vec(i) = normalize(vec(i));
    }
}

/**
 * Returns a power column vector of the given parameter such that
 * power vector = [1; par^1; ...; par^(degree-1); par^degree]
 *
 * @param[in] par parameter
 * @return vector of parameter
 */
Eigen::VectorXd Bs::get_pow_vec(double par) {
    return get_pow_vec(par, this->coeffs.get_degree());
}

/**
 * Returns a power column vector of length (this->coeffs.get_degree() + 1)
 * of the given parameter such that
 * power vector = [0, 0, ...; par^(highest_pow-1); par^highest_pow]
 *
 * @param[in] par parameter
 * @param[in] highest_pow highest power in parameter vector
 * @return vector of parameter
 */
Eigen::VectorXd Bs::get_pow_vec(double par, int highest_pow) {
    Eigen::VectorXd pow_vec = Eigen::VectorXd::Zero(this->coeffs.get_degree() + 1);
    int power = highest_pow;
    int i = 0;
    while( power >= 0) {
        pow_vec(i) = pow((double) par, (int) power);
        power--;
        i++;
    }
    return pow_vec;
}

/**
 * Returns a column vector of length (this->coeffs.get_degree() + 1)
 * of exponents such that
 * power exp_vec = [highest_exp; highest_exp-1; ... 1; 0; 0]
 *
 * @param[in] highest_exp highest exponent
 * @return exponent vector
 */
Eigen::VectorXd Bs::get_exp_vec(int highest_exp) {
    Eigen::VectorXd exp_vec = Eigen::VectorXd::Zero(this->coeffs.get_degree() + 1);
    int exp = highest_exp;
    int i = 0;
    while(exp >= 0) {
        exp_vec(i) = exp;
        exp--;
        i++;
    }
    return exp_vec;
}

/**
 * Returns the curve's value at the given parameter.
 *
 * @param[in] parameter
 * @return value at current parameter
 */
double Bs::get_value(double parameter) {
    return get_value(parameter, 0);
}

/**
 * Returns the curve's value at the given parameters.
 *
 * @param[in] parameters
 * @return value at current parameter
 */
Eigen::RowVectorXd Bs::get_values(Eigen::RowVectorXd& parameters) {
    printf("muuuh\n");
    Eigen::RowVectorXd values = Eigen::RowVectorXd::Zero(parameters.size());
    for(int i = 0; i < parameters.size(); i++) {
        values(i) = get_value(parameters(i), 0);
    }
    return values;
}

/**
 * Returns the curve's derivative at the given parameter.
 *
 * @param[in] parameter
 * @param[in] derivative order of derivative
 * @return value at current parameter
 */
double Bs::get_value(double parameter, unsigned int derivative) {
    if(this->coeffs.get_degree() < derivative) {
        return 0;
    }
    parameter = this->normalize(parameter);
    Eigen::VectorXd exp_vec = Eigen::VectorXd::Ones(this->coeffs.get_degree() + 1); // vector of exponents for derivatives
    
    for(unsigned int i = this->coeffs.get_degree(); i > this->coeffs.get_degree() - derivative; i--) {
        Eigen::VectorXd this_exp_vec = Bs::get_exp_vec(i);
        exp_vec = exp_vec.cwiseProduct(this_exp_vec);
    }
    Eigen::VectorXd pow_vec = get_pow_vec(parameter, this->coeffs.get_degree() - derivative);
	return 1./pow((double) (this->par_end - this->par_start), (int) derivative)*ctrl_pts*(*(coeffs.get_coeff_mat(parameter)))*exp_vec.cwiseProduct(pow_vec);
}


/**
 * Returns the curve's derivative at the given parameters.
 *
 * @param[in] parameters
 * @param[in] derivatives vector of orders of derivative
 * @return value at current parameter
 */
Eigen::RowVectorXd Bs::get_values(Eigen::RowVectorXd& parameters, Eigen::VectorXd& derivatives) {
    Eigen::RowVectorXd values = Eigen::RowVectorXd::Zero(parameters.size());
    for(int i = 0; i < parameters.size(); i++) {
        values(i) = get_value(parameters(i), (int) derivatives(i));
    }
    return values;
//     printf("moeh new_for\n");
//     Eigen::MatrixXd mat(this->coeffs.get_degree() + 1, parameters.size());
//     mat.setZero();
//     for(int i = 0; i < parameters.size(); i++) {
//         unsigned int derivative = (int) derivatives(i);
//         if(this->coeffs.get_degree() < derivative) { // result will be 0, skip calculation
//             mat.col(i) = Eigen::VectorXd::Zero(this->coeffs.get_degree() + 1);
//             continue;
//         }
//         printf("moeh_ got der %d\n", i);
//         double parameter = this->normalize(parameters(i));
//         printf("moeh_ got par %d\n", i);
//         Eigen::VectorXd exp_vec = Eigen::VectorXd::Ones(this->coeffs.get_degree() + 1); // vector of exponents for derivatives
//         printf("moeh_ got exp_vec %d\n", i);
//         for(unsigned int j = this->coeffs.get_degree(); j > this->coeffs.get_degree() - derivative; j--) {
//             Eigen::VectorXd this_exp_vec = Bs::get_exp_vec(j);
//             exp_vec = exp_vec.cwiseProduct(this_exp_vec);
//         }
//         Eigen::VectorXd pow_vec = 1./pow((double) (this->par_end - this->par_start), (int) derivative)*get_pow_vec(parameter, this->coeffs.get_degree() - derivative);
//         mat.col(i) = (*(coeffs.get_coeff_mat(parameter)))*exp_vec*pow_vec;
//     }
//     return ctrl_pts*mat;
}

/**
 * Sets start parameter.
 *
 * @param[in] par_start    start parameter
 */
void Bs::set_par_start(double par_start) {
    this->par_start = par_start;
}

/**
 * Sets end parameter.
 *
 * @param[in] par_end    start parameter
 */
void Bs::set_par_end(double par_end) {
    this->par_end = par_end;
}

/**
 * Returns lower bound of parameters.
 *
 * @return lower bound of parameters
 */
double Bs::get_par_start() {
    return this->par_start;
}

/**
 * Returns upper bound of parameters.
 *
 * @return upper bound of parameters
 */
double Bs::get_par_end() {
    return this->par_end;
}

/**
 * Returns vector of control poitns.
 *
 * @return vector of control poitns
 */
Eigen::RowVectorXd Bs::get_ctrl_pts() {
    return this->ctrl_pts;
}

/**
 * Returns coefficient structure.
 *
 * @return coefficient structure
 */
Bs_coeffs* Bs::get_coeffs() {
    return &(this->coeffs);
}

void init_Bs() {
    //cout << "Test successful" << endl;
    printf("bla\n");
}
