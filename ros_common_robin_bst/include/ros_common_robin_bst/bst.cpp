#include "Eigen/Dense" // MatrixXd
#include "Bs.h"
#include "mex.h" // MATLAB mex

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>, Eigen::RowMajor> MapType; //!< maps row major C-array to Eigen::Matrix
typedef Eigen::Map<Eigen::Matrix<double,1,Eigen::Dynamic> > RowMapType; //!< maps C-array to Eigen::Matrix row vector
typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > ColMapType; //!< maps C-array to Eigen::Matrix column vector
 
void mexFunction(
        int          nlhs,
        mxArray      *plhs[],
        int          nrhs,
        const mxArray *prhs[]
        );

/**
 * Gateway function for MATLAB.
 *
 * @param[in] nlhs Number of expected output mxArrays
 * @param[out] plhs Array of pointers to the expected output mxArrays
 * @param[in] nrhs Number of input mxArrays
 * @param[in] prhs Array of pointers to the input mxArrays.
 */
void mexFunction(
        int          nlhs,
        mxArray      *plhs[],
        int          nrhs,
        const mxArray *prhs[]
        )
{
    
    /* Check for proper number of arguments */
    
    if (nrhs < 2 || nrhs > 4) { // incorrect number of inputs
        mexErrMsgIdAndTxt("bst:nargin",
                "This function requires between two and four inputs.");
        return;
    }
    if(nlhs > 1) { // incorrect number of outputs
        mexErrMsgIdAndTxt("bst:wrong_number_outputs",
                "This function requires at most one output argument.");
    }
    // correct number of arguments
    
    // create spline from given data
    Bs spline;
    if(mxIsStruct(prhs[0]) && (nrhs == 2 || nrhs == 3)) { // warm start
        if(mxGetField(prhs[0], 0, "degree") == NULL || 
                mxGetField(prhs[0], 0, "ctrl_pts") == NULL ||
                mxGetField(prhs[0], 0, "par_start") == NULL || 
                mxGetField(prhs[0], 0, "par_end") == NULL ||
                mxGetField(prhs[0], 0, "knots") == NULL ||
                mxGetField(prhs[0], 0, "coeffs") == NULL) {
            mexErrMsgIdAndTxt("bst:missing_field",
                    "The first parameter needs to have the fields 'degree', 'ctrl_pts', 'par_start', 'par_end', 'knots' and 'coeffs'.");
            return;
        }
        if(!mxIsNumeric(mxGetField(prhs[0], 0, "degree")) || mxGetNumberOfElements(mxGetField(prhs[0], 0, "degree")) > 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                "The first parameter's field 'degree' needs to be the (scalar) degree of the spline.");
            return;
        }
        if(!mxIsNumeric(mxGetField(prhs[0], 0, "ctrl_pts")) || mxGetNumberOfElements(mxGetField(prhs[0], 0, "ctrl_pts")) < 2) {
            mexErrMsgIdAndTxt("bst:nargin",
                "The first parameter's field 'ctrl_pts' needs to be the vector of control points of length > 1.");
            return;
        }
        if(!mxIsNumeric(mxGetField(prhs[0], 0, "par_start")) || mxGetNumberOfElements(mxGetField(prhs[0], 0, "par_start")) > 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                "The first parameter's field 'par_start' needs to be the (scalar) lower bound of the spline parameters.");
            return;
        }
        if(!mxIsNumeric(mxGetField(prhs[0], 0, "par_end")) || mxGetNumberOfElements(mxGetField(prhs[0], 0, "par_end")) > 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                "The first parameter's field 'par_end' needs to be the (scalar) upper bound of the spline parameters.");
            return;
        }
        if(!mxIsNumeric(mxGetField(prhs[0], 0, "knots")) || mxGetNumberOfElements(mxGetField(prhs[0], 0, "knots")) < 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                "The first parameter's field 'knots' needs to be the vector of knots of length > 1.");
            return;
        }
        if(!mxIsNumeric(mxGetField(prhs[0], 0, "coeffs")) || mxGetNumberOfDimensions(mxGetField(prhs[0], 0, "coeffs")) > 3 || mxGetNumberOfDimensions(mxGetField(prhs[0], 0, "coeffs")) < 2 || mxGetNumberOfElements(mxGetField(prhs[0], 0, "coeffs")) < 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                "The first parameter's field 'coeffs' needs to a non-empty 3D array of doubles.");
            return;
        }
        std::vector<Eigen::MatrixXd> coeff_mat_vec;
        const int* coeffs_dim = mxGetDimensions(mxGetField(prhs[0], 0, "coeffs"));
        double* coeffs_data = mxGetPr(mxGetField(prhs[0], 0, "coeffs"));
        if(mxGetNumberOfDimensions(mxGetField(prhs[0], 0, "coeffs")) == 2) {
			coeff_mat_vec.push_back(Eigen::Map<Eigen::MatrixXd>(coeffs_data, coeffs_dim[0], coeffs_dim[1]));
		} else {
			for(int i = 0; i < coeffs_dim[2]; i++) {
				coeff_mat_vec.push_back(Eigen::Map<Eigen::MatrixXd>(coeffs_data+i*coeffs_dim[0]*coeffs_dim[1],coeffs_dim[0],coeffs_dim[1]));
			}
		}
        spline = Bs::Bs((unsigned int) mxGetScalar(mxGetField(prhs[0], 0, "degree")),
                (Eigen::RowVectorXd) RowMapType(mxGetPr(mxGetField(prhs[0], 0, "ctrl_pts")),mxGetNumberOfElements(mxGetField(prhs[0], 0, "ctrl_pts"))),
                *mxGetPr(mxGetField(prhs[0], 0, "par_start")),
                *mxGetPr(mxGetField(prhs[0], 0, "par_end")),
                Bs_coeffs::Bs_coeffs((unsigned int) mxGetScalar(mxGetField(prhs[0], 0, "degree")),
                                     (Eigen::RowVectorXd) RowMapType(mxGetPr(mxGetField(prhs[0], 0, "knots")),mxGetNumberOfElements(mxGetField(prhs[0], 0, "knots"))),
                                     coeff_mat_vec)); // vector of coefficient matrices
        if(!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) < 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                    "The second argument needs to be a finite vector of parameters with size > 0.");
            return;
        }
        Eigen::VectorXd der;
        if(nrhs == 2) {
            der = Eigen::VectorXd::Zero(mxGetNumberOfElements(prhs[1]));
        } else { // nrhs == 3, vector of derivative orders is available
            der = (Eigen::VectorXd) ColMapType(mxGetPr(prhs[2]), mxGetNumberOfElements(prhs[2]));
            if(!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2]) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The third argument needs to be a finite vector of derivatives with size > 0.");
                return;
            }
            if(mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2])) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The length of the vector of parameters (%d) is not equal the length of the vector of orders of derivatives (%d)", mxGetNumberOfElements(prhs[1]), mxGetNumberOfElements(prhs[2]));
                return;
            }
        }
        Eigen::RowVectorXd res = spline.get_values((Eigen::RowVectorXd) RowMapType(mxGetPr(prhs[1]), mxGetNumberOfElements(prhs[1])), der);
        
        // return result
        plhs[0] = mxCreateDoubleMatrix((int) res.size(), 1, mxREAL);
        for(int i = 0; i < res.size(); i++) {
            mxGetPr(plhs[0])[i] = res(i);
        }
        return;
    }
    if(!mxIsNumeric(prhs[0]) || mxGetNumberOfElements(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("bst:nargin",
                "The first argument needs to be the (scalar) degree of the spline.");
        return;
    }
    if(nrhs == 2) { // create spline of degree on interval [0,1]
        if(mxIsNumeric(prhs[1])) {
            if(mxGetNumberOfElements(prhs[1]) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The second argument needs to be the finite vector of control points with size > 0.");
                return;
            }
            spline = Bs::Bs((unsigned int) mxGetScalar(prhs[0]), // degree
                    (Eigen::RowVectorXd) RowMapType(mxGetPr(prhs[1]),mxGetNumberOfElements(prhs[1]))); // ctrl_pts
        } else {
            mexErrMsgIdAndTxt("bst:wrong_type",
                    "The second parameter needs to be the vector of control points.");
            return;
        }
    } else if(nrhs == 3) {
        if(mxIsNumeric(prhs[1])) { // second argument is a numeric vector not an approximation structure
            if(mxGetNumberOfElements(prhs[1]) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The second argument needs to be the finite vector of control points with size > 0.");
                return;
            }
            unsigned int degree = (unsigned int) mxGetScalar(prhs[0]);
            if(mxGetNumberOfElements(prhs[2]) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The third argument needs to be the finite vector of knots with size > 0.");
                return;
            }
            Eigen::RowVectorXd knots = (Eigen::RowVectorXd) RowMapType(mxGetPr(prhs[2]),mxGetNumberOfElements(prhs[2]));
            /*for(int i = 1; i <= degree; i++) {
                if(fabs(knots(0) - knots(i)) > NUMEPS) {
                    mexErrMsgIdAndTxt("bst:knots_problem",
                        "The first (degree + 1) knots must be equal.");
                    return;
                }
            }
            for(int i = knots.size() - degree + 1; i < knots.size(); i++) {
                if(fabs(knots(knots.size() - degree) - knots(i)) > NUMEPS) {
                    mexErrMsgIdAndTxt("bst:knots_problem",
                        "The last (degree + 1) knots must be equal.");
                    return;
                }
            }*/
            for(int i = 1; i < knots.size(); i++) {
                if(knots(i) - knots(i-1) < -NUMEPS) {
                    mexErrMsgIdAndTxt("bst:knots_problem",
                        "The knots vector must be strictly monotonic.");
                    return;
                }
            }
            if((int) mxGetNumberOfElements(prhs[2]) - (int) mxGetScalar(prhs[0]) - 1 != (int) mxGetNumberOfElements(prhs[1])) {
                mexErrMsgIdAndTxt("bst:wrong_parameter_relation",
                        "Wrong size. The following equation must hold:\n    number of knots (%d) - degree (%d) - 1 = number of control points (%d).", mxGetNumberOfElements(prhs[2]), (unsigned int) mxGetScalar(prhs[0]), mxGetNumberOfElements(prhs[1]));
                return;
            }
            spline = Bs::Bs(degree, // degree
                    (Eigen::RowVectorXd) RowMapType(mxGetPr(prhs[1]),mxGetNumberOfElements(prhs[1])), // ctrl_pts
                    knots); // knots
        } else { // second argument is an approximation structure
            if(mxGetField(prhs[1], 0, "par") == NULL || mxGetField(prhs[1], 0, "val") == NULL || mxGetField(prhs[1], 0, "der") == NULL) {
                mexErrMsgIdAndTxt("bst:approx_field",
                        "The second parameter needs to be a structure with the fields 'par', 'val' and 'der'");
                return;
            }
            if(!mxIsNumeric(mxGetField(prhs[1], 0, "par")) || mxGetNumberOfElements(mxGetField(prhs[1], 0, "par")) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The field 'par' of the second argument needs to be the finite vector of parameters with size > 0.");
                return;
            }
            if(!mxIsNumeric(mxGetField(prhs[1], 0, "val")) || mxGetNumberOfElements(mxGetField(prhs[1], 0, "val")) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The field 'val' of the second argument needs to be the finite vector of values with size > 0.");
                return;
            }
            if(!mxIsNumeric(mxGetField(prhs[1], 0, "der")) || mxGetNumberOfElements(mxGetField(prhs[1], 0, "der")) < 1) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The field 'der' of the second argument needs to be the finite vector of derivative orders with size > 0.");
                return;
            }
            if(mxGetNumberOfElements(mxGetField(prhs[1], 0, "der")) != mxGetNumberOfElements(mxGetField(prhs[1], 0, "val")) ||
                    mxGetNumberOfElements(mxGetField(prhs[1], 0, "val")) != mxGetNumberOfElements(mxGetField(prhs[1], 0, "par"))) {
                mexErrMsgIdAndTxt("bst:nargin",
                        "The lengths of the fields of the second argument are not equal.");
                return;
            }
            Eigen::MatrixXd approx_data(mxGetNumberOfElements(mxGetField(prhs[1], 0, "val")),3);
            approx_data.col(0) = (Eigen::VectorXd) ColMapType(mxGetPr(mxGetField(prhs[1], 0, "par")), mxGetNumberOfElements(mxGetField(prhs[1], 0, "par")));
            approx_data.col(1) = (Eigen::VectorXd) ColMapType(mxGetPr(mxGetField(prhs[1], 0, "val")), mxGetNumberOfElements(mxGetField(prhs[1], 0, "val")));
            approx_data.col(2) = (Eigen::VectorXd) ColMapType(mxGetPr(mxGetField(prhs[1], 0, "der")), mxGetNumberOfElements(mxGetField(prhs[1], 0, "der")));
            
            spline = Bs::Bs((unsigned int) mxGetScalar(prhs[0]), // degree
                    approx_data,
                    (unsigned int) mxGetScalar(prhs[2])); // n_ctrl_pts
        }
    } else if(nrhs == 4) {
        if(!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) < 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                    "The second argument needs to be the finite vector of control points with size > 0.");
            return;
        }
        if(!mxIsNumeric(prhs[0]) || mxGetNumberOfElements(prhs[2]) > 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                    "The third argument needs to be the (scalar) lower bound of the spline parameters.");
            return;
        }
        if(!mxIsNumeric(prhs[0]) || mxGetNumberOfElements(prhs[3]) > 1) {
            mexErrMsgIdAndTxt("bst:nargin",
                    "The fourth argument needs to be the (scalar) upper bound of the spline parameters.");
            return;
        }
        if(mxGetScalar(prhs[2]) >= mxGetScalar(prhs[3])) {
            mexErrMsgIdAndTxt("bst:wrong_parameter_relation",
                    "The upper bound of the spline parameters (%f) needs to be greater than the lower bound (%f).", mxGetScalar(prhs[3]), mxGetScalar(prhs[2]));
            return;
        }
        spline = Bs::Bs((unsigned int) mxGetScalar(prhs[0]), // degree
                (Eigen::RowVectorXd) RowMapType(mxGetPr(prhs[1]),mxGetNumberOfElements(prhs[1])), // ctrl_pts
                mxGetScalar(prhs[2]), mxGetScalar(prhs[3])); // par_start and par_end
    }
    
    mwSize dims[2] = {1, 1};
    const char *field_names[] = {"degree", "ctrl_pts", "par_start", "par_end", "knots", "coeffs"};
    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);
    
    mxSetField(plhs[0], 0, "par_start", mxCreateDoubleScalar(spline.get_par_start()));
    mxSetField(plhs[0], 0, "par_end", mxCreateDoubleScalar(spline.get_par_end()));
    Eigen::RowVectorXd ctrl_pts = spline.get_ctrl_pts();
    mxArray* ctrl_pts_out = mxCreateDoubleMatrix((int) ctrl_pts.size(), 1, mxREAL);
    for(int i = 0; i < ctrl_pts.size(); i++) {
        mxGetPr(ctrl_pts_out)[i] = ctrl_pts(i);
    }
    mxSetField(plhs[0], 0, "ctrl_pts", ctrl_pts_out);
    
    Bs_coeffs* coeffs = spline.get_coeffs();
    mxSetField(plhs[0], 0, "degree", mxCreateDoubleScalar(coeffs->get_degree()));
    Eigen::RowVectorXd knots = coeffs->get_knots();
    mxArray* knots_out = mxCreateDoubleMatrix((int) knots.size(), 1, mxREAL);
    for(int i = 0; i < knots.size(); i++) {
        mxGetPr(knots_out)[i] = knots(i);
    }
    mxSetField(plhs[0], 0, "knots", knots_out);
    
    std::vector<Eigen::MatrixXd> coeff_mat_vec = coeffs->get_coeff_mat_vec();
    int coeff_out_size[3];
    coeff_out_size[0] = (int) coeff_mat_vec[0].rows();
    coeff_out_size[1] = (int) coeff_mat_vec[0].cols();
    coeff_out_size[2] = (int) coeff_mat_vec.size();
    mxArray* coeffs_out = mxCreateNumericArray(3,coeff_out_size, mxDOUBLE_CLASS,mxREAL);
    double* coeffs_ptr = (double*) mxGetPr(coeffs_out);
    for(int i = 0; i < coeff_out_size[2]; i++) {
        for(int k = 0; k < coeff_out_size[1]; k++) {
            for(int j = 0; j < coeff_out_size[0]; j++) {
                coeffs_ptr[j + k*coeff_out_size[0] + i*coeff_out_size[0]*coeff_out_size[1]] = coeff_mat_vec[i](j,k);
            }
        }
    }
    /*mxArray* coeffs_out = mxCreateCellMatrix(1, (int) coeff_mat_vec.size());
    for(int i = 0; i < coeff_mat_vec.size(); i++) {
        mxArray* mat = mxCreateDoubleMatrix((int) coeff_mat_vec[i].rows(), (int) coeff_mat_vec[i].cols(), mxREAL);
        for(int j = 0; j < coeff_mat_vec[i].rows()*coeff_mat_vec[i].cols(); j++) {
            mxGetPr(mat)[j] = coeff_mat_vec[i](j);
        }
        mxSetCell(coeffs_out, i, mat);
    }*/
    mxSetField(plhs[0], 0, "coeffs", coeffs_out);
    
    return;
}