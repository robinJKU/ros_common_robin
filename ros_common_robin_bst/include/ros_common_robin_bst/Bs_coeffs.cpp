#include "Bs_coeffs.h"


/**
 * Returns the coefficient vector of the resulting polynomial of
 * the multiplication of poly1 and poly2.
 * Coefficient order starts with lower orders, e.g.
 * [1, 2, 3] means 1 + 2x + 3x^2
 *
 * @param[in] poly1 coefficient vector a polynomial
 * @param[in] poly2 coefficient vector another polynomial
 * @return resulting coefficient vector
 */
Eigen::RowVectorXd Bs_coeffs::poly_mult(Eigen::RowVectorXd& poly1, Eigen::RowVectorXd& poly2) {
    Eigen::RowVectorXd poly = Eigen::RowVectorXd::Zero((int) poly1.cols() + (int) poly2.cols() - 1);  // init convolution result vector with zeros
    for (int i = 0; i < poly1.cols(); i++) {
        for (int j = 0; j < poly2.cols(); j++) {
            poly(i+j) += poly1(i) * poly2(j);
        }
    }
    return poly;
}

/**
 * Computes spline basis function N_j^d by means of de Boor-Cox-Mansfield algorithm
 * and returns coefficient vector
 *
 * @param[in] j index of basis function (zero-based)
 * @param[in] d (local) degree of basis function
 * @param[in] interval interval number (zero-based)
 * @return resulting vector of polynomial coefficients
 */
Eigen::RowVectorXd Bs_coeffs::get_coeff_vec(unsigned int j, unsigned int d, unsigned int interval) {
    Eigen::RowVectorXd coeff_vec = Eigen::RowVectorXd::Zero(d+1);
    if(d == 0) { // simple case
        if((knots(j) <= knots(interval) && knots(interval) < knots(j+1)) ) {
            coeff_vec(0) = 1.;
        } else {
            coeff_vec(0) = 0.;
        }
    } else { // general case, d >= 1
        double left_den = knots(j+d) - knots(j);
        if(fabs(left_den) > NUMEPS) { // 0/0 := 0
            Eigen::RowVectorXd left = Eigen::RowVectorXd::Zero(2);
            left(0) = 1.;
            left(1) = -knots(j);
            Eigen::RowVectorXd right = get_coeff_vec(j, d-1, interval);
            left = poly_mult(left, right);
            left.conservativeResize(1, d+1);
            coeff_vec = coeff_vec + left/left_den;
        }
        double right_den = knots(j+d+1) - knots(j+1);
        if(fabs(right_den) > NUMEPS) { // 0/0 := 0
            Eigen::RowVectorXd right = Eigen::RowVectorXd::Zero(2);
            right(0) = -1.;
            right(1) = knots(j+d+1);
            Eigen::RowVectorXd left = get_coeff_vec(j+1, d-1, interval);
            right = poly_mult(right, left);
            right.conservativeResize(1, d+1);
            coeff_vec = coeff_vec + right/right_den;
        }
    }
    return coeff_vec;
}

/**
 * Returns the unique values of the sorted knots vector of *this.
 *
 * @return vector of unique knots
 */
Eigen::RowVectorXd Bs_coeffs::get_unique_knots() {
    Eigen::RowVectorXd unique_knots = Eigen::RowVectorXd::Zero(1);
    unique_knots << this->knots(0);
    for(int i = 1; i < knots.size(); i++) {
        if(fabs(this->knots(i) - unique_knots(unique_knots.size()-1)) > NUMEPS) {
            unique_knots.conservativeResize(unique_knots.size()+1);
            unique_knots(unique_knots.size()-1) = knots(i);
        }
    }
    return unique_knots;
}

/**
 * Generates coefficient matrices for all non-zero knot intervals.
*/
void Bs_coeffs::generate_coeff_mat_vec() {
    Eigen::RowVectorXd unique_knots = this->get_unique_knots();
    this->coeff_mat_vec = std::vector<Eigen::MatrixXd>(unique_knots.size()-1);
    for(unsigned int i = 0; i < this->coeff_mat_vec.size(); i++) { // for all unique intervals
        this->coeff_mat_vec[i] = Eigen::MatrixXd::Zero(this->knots.size() - degree - 1, degree+1);
    }
    for(unsigned int i = 0; i < this->knots.size()-1; i++) { // for all intervals
        for(unsigned int j = 0; j < this->knots.size() - degree - 1; j++) {
            this->coeff_mat_vec[get_unique_interval(knots(i))].row(j) << get_coeff_vec(j, degree, i);
        }
        // skip knots with multiplicity > 1
        int current_unique_interval = get_unique_interval(knots(i));
        while(i < knots.size() && current_unique_interval == get_unique_interval(knots(i))) {
            i++;
        }
        i--;
    }
}

/**
 * (Empty) default Constructor.
 *
 */
Bs_coeffs::Bs_coeffs() {
    return;
}

/**
 * Constructor.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] knots knots vector
 */
Bs_coeffs::Bs_coeffs(unsigned int degree, Eigen::RowVectorXd& knots) {
    this->degree = degree;
    this->knots = knots;
    generate_coeff_mat_vec();
}

/**
 * Warm start constructor.
 *
 * @param[in] degree maximum degree of curve
 * @param[in] knots knots vector
 * @param[in] coeff_mat_vec vector of coefficient matrices
 */
Bs_coeffs::Bs_coeffs(unsigned int degree, Eigen::RowVectorXd& knots, std::vector<Eigen::MatrixXd>& coeff_mat_vec) {
    this->degree = degree;
    this->knots = knots;
    this->coeff_mat_vec = coeff_mat_vec;
}

/**
 * Returns interval number of vector of unique knots for current parameter.
 *
 * @param[in] parameter parameter
 * @return interval number
 */
unsigned int Bs_coeffs::get_unique_interval(double parameter) {
    if(parameter < 0 || parameter > 1) {
        //throw x_wrong_range;
    }
    Eigen::RowVectorXd unique_knots = get_unique_knots();
    int interval = 0;
    while(interval < unique_knots.size() - 1) {
        if(unique_knots(interval) > parameter) return interval - 1;
        interval++;
    }
    return (unsigned int) (unique_knots.size() - 2);
}

/**
 * Returns interval number of vector of knots for current parameter.
 *
 * @param[in] parameter parameter
 * @return interval number
 */
unsigned int Bs_coeffs::get_interval(double parameter) {
    if(parameter < 0 || parameter > 1) {
        //throw x_wrong_range;
    }
    int interval = 0;
    while(interval < knots.size() - 1) {
        if(knots(interval) > parameter) return interval - 1;
        interval++;
    }
    return (unsigned int) (knots.size() - 1);
}

/**
 * Returns a pointer to the coefficient matrix for the given parameter.
 *
 * @param[in] parameter parameter
 * @return pointer to coefficient matrix
 */
Eigen::MatrixXd* Bs_coeffs::get_coeff_mat(double parameter) {
    unsigned int interval;
    try {
        interval = get_unique_interval(parameter);
    } catch (exception& e) {
        cout << e.what() << endl;
        return NULL;
    }
    return &(coeff_mat_vec[interval]);
}

/**
 * Returns degree of spline.
 *
 * @return degree of spline
 */
unsigned int Bs_coeffs::get_degree() {
    return this->degree;
}

/**
 * Returns vector of coefficient matrices.
 *
 * @return vector of coefficient matrices
 */
std::vector<Eigen::MatrixXd> Bs_coeffs::get_coeff_mat_vec() {
    return this->coeff_mat_vec;
}

/**
 * Returns vector of knots.
 *
 * @return vector of knots
 */
Eigen::RowVectorXd Bs_coeffs::get_knots() {
    return this->knots;
}

