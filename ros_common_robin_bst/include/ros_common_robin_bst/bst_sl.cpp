#include <iostream>
#include <vector>
#include "Eigen/Dense" // MatrixXd
#include "Bs.h"

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>, Eigen::RowMajor> MapType; //!< maps row major C-array to Eigen::Matrix
typedef Eigen::Map<Eigen::Matrix<double,1,Eigen::Dynamic> > RowMapType; //!< maps C-array to Eigen::Matrix row vector
typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > ColMapType; //!< maps C-array to Eigen::Matrix column vector

#ifdef __cplusplus
extern "C" { // use the C fcn-call standard for all functions
#endif       // defined within this scope
    
#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  bst_sl
    
    /*
     * Need to include simstruc.h for the definition of the SimStruct and
     * its associated macro definitions.
     */
#include "simstruc.h"
    
    /*====================*
     * S-function methods *
     *====================*/
    
    /* Function: mdlInitializeSizes ===============================================
     * Abstract:
     *    The sizes information is used by Simulink to determine the S-function
     *    block's characteristics (number of inputs, outputs, states, etc.).
     */
    static void mdlInitializeSizes(SimStruct *S)
    {
        /* See sfuntmpl_doc.c for more details on the macros below */
        
        if (!ssSetNumOutputPorts(S, 1)) return;
        
        ssSetNumPWork(S, 2); // reserve elements in the pointers vector
        
        ssSetNumSFcnParams(S, 2);
        if (ssGetSFcnParamsCount(S) != 2) {
            return;
        }
        
        //ssSetOutputPortWidth(S, 0, 4);
        const mxArray* der_param = ssGetSFcnParam(S, 1);
        if(mxIsNumeric(der_param) && mxGetNumberOfDimensions(der_param) == 2) {
            ssSetOutputPortWidth(S, 0, mxGetNumberOfElements(der_param));
        } else {
            return;
        }
            
        ssSetNumContStates(S, 0);
        ssSetNumDiscStates(S, 0);
        
        if (!ssSetNumInputPorts(S, 1)) return;
        ssSetInputPortWidth(S, 0, 1);
        ssSetInputPortDirectFeedThrough(S, 0, 1);
        
        ssSetNumSampleTimes(S, 1);
        ssSetNumRWork(S, 0);
        ssSetNumIWork(S, 0);
        ssSetNumModes(S, 0); 
        ssSetNumNonsampledZCs(S, 0);
        
        ssSetOptions(S, 0);
    }
    
    
    
    
    /* Function: mdlInitializeSampleTimes =========================================
     * Abstract:
     *    This function is used to specify the sample time(s) for your
     *    S-function. You must register the same number of sample times as
     *    specified in ssSetNumSampleTimes.
     */
    static void mdlInitializeSampleTimes(SimStruct *S)
    {
        ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
        ssSetOffsetTime(S, 0, 0.0);
        
    }
    
#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START)
    /* Function: mdlStart =======================================================
     * Abstract:
     *    This function is called once at start of model execution. If you
     *    have states that should be initialized once, this is the place
     *    to do it.
     */
    static void mdlStart(SimStruct *S)
    {
        const mxArray* spline_struct = ssGetSFcnParam(S, 0);
        
        if(!mxIsStruct(spline_struct)) { // check if parameter is a struct
            return;
        }
        if(mxGetField(spline_struct, 0, "degree") == NULL || // check if all required fields exist
                mxGetField(spline_struct, 0, "ctrl_pts") == NULL ||
                mxGetField(spline_struct, 0, "par_start") == NULL ||
                mxGetField(spline_struct, 0, "par_end") == NULL ||
                mxGetField(spline_struct, 0, "knots") == NULL ||
                mxGetField(spline_struct, 0, "coeffs") == NULL) {
            return;
        }
        // check if fields are valid
        if(!mxIsNumeric(mxGetField(spline_struct, 0, "degree")) || mxGetNumberOfElements(mxGetField(spline_struct, 0, "degree")) > 1) {
            mexErrMsgTxt("spline.degree must be numeric and of length 1.\n");
            return;
        }
        if(!mxIsNumeric(mxGetField(spline_struct, 0, "ctrl_pts")) || mxGetNumberOfElements(mxGetField(spline_struct, 0, "ctrl_pts")) < 2) {
            mexErrMsgTxt("spline.ctrl_pts must be numeric and of length >2.\n");
            return;
        }
        if(!mxIsNumeric(mxGetField(spline_struct, 0, "par_start")) || mxGetNumberOfElements(mxGetField(spline_struct, 0, "par_start")) > 1) {
            mexErrMsgTxt("spline.par_start must be numeric and of length 1.\n");
            return;
        }
        if(!mxIsNumeric(mxGetField(spline_struct, 0, "par_end")) || mxGetNumberOfElements(mxGetField(spline_struct, 0, "par_end")) > 1) {
            mexErrMsgTxt("spline.par_end must be numeric and of length 1.\n");
            return;
        }
        if(!mxIsNumeric(mxGetField(spline_struct, 0, "knots")) || mxGetNumberOfElements(mxGetField(spline_struct, 0, "knots")) < 1) {
            mexErrMsgTxt("spline.par_end must be numeric and non-empty.\n");
            return;
        }
        if(!mxIsNumeric(mxGetField(spline_struct, 0, "coeffs")) || mxGetNumberOfDimensions(mxGetField(spline_struct, 0, "coeffs")) != 3 || mxGetNumberOfElements(mxGetField(spline_struct, 0, "coeffs")) < 1) {
            mexErrMsgTxt("spline.par_end must be numeric, of dimension 3 and non-empty.\n");
            return;
        }
        // parameter is valid
        std::vector<Eigen::MatrixXd> coeff_mat_vec;
        const int* coeffs_dim = mxGetDimensions(mxGetField(spline_struct, 0, "coeffs"));
        double* coeffs_data = mxGetPr(mxGetField(spline_struct, 0, "coeffs"));
        for(int i = 0; i < coeffs_dim[2]; i++) {
            coeff_mat_vec.push_back(Eigen::Map<Eigen::MatrixXd>(&coeffs_data[i*coeffs_dim[0]*coeffs_dim[1]],coeffs_dim[0],coeffs_dim[1]));
        }
        ssGetPWork(S)[0] = (void *) new Bs((unsigned int) mxGetScalar(mxGetField(spline_struct, 0, "degree")),
                (Eigen::RowVectorXd) RowMapType(mxGetPr(mxGetField(spline_struct, 0, "ctrl_pts")),mxGetNumberOfElements(mxGetField(spline_struct, 0, "ctrl_pts"))),
                *mxGetPr(mxGetField(spline_struct, 0, "par_start")),
                *mxGetPr(mxGetField(spline_struct, 0, "par_end")),
                Bs_coeffs::Bs_coeffs((unsigned int) mxGetScalar(mxGetField(spline_struct, 0, "degree")),
                (Eigen::RowVectorXd) RowMapType(mxGetPr(mxGetField(spline_struct, 0, "knots")),mxGetNumberOfElements(mxGetField(spline_struct, 0, "knots"))),
                coeff_mat_vec)); // vector of coefficient matrices
        const mxArray* der_param = ssGetSFcnParam(S, 1);
        double* der_mask = mxGetPr(der_param);
        ssGetPWork(S)[1] = (void*) new Eigen::RowVectorXd(mxGetNumberOfElements(der_param));
        Eigen::RowVectorXd* der = (Eigen::RowVectorXd*) ssGetPWork(S)[1];
        
        for(size_t i = 0; i < mxGetNumberOfElements(der_param); i++) {
            (*der)[i] = der_mask[i];
        }  
    }
#endif /*  MDL_START */
    
    /* Function: mdlOutputs =======================================================
     * Abstract:
     *    In this function, you compute the outputs of your S-function
     *    block. Generally outputs are placed in the output vector, ssGetY(S).
     */
    static void mdlOutputs(SimStruct *S, int_T tid) {   
        InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
        
        real_T  *y = ssGetOutputPortRealSignal(S,0);
        int_T y_size = ssGetOutputPortWidth(S,0);
        
        Bs* spline = (Bs*) ssGetPWork(S)[0];
        // check parameter
        if(spline->get_par_start() > *uPtrs[0] || spline->get_par_end() < *uPtrs[0]) {
            mexErrMsgTxt("Input is out of spline range\n");
            return;
        }
        
        Eigen::RowVectorXd params = Eigen::RowVectorXd::Ones(y_size)*(*uPtrs[0]);
        Eigen::VectorXd* der = (Eigen::VectorXd*) ssGetPWork(S)[1];
        Eigen::RowVectorXd res = spline->get_values(params, *der);
        for(int_T i = 0; i < y_size; i++) {
            y[i] = res(i);
        }
    }                                                
    
    /* Function: mdlTerminate =====================================================
     * Abstract:
     *    In this function, you should perform any actions that are necessary
     *    at the termination of a simulation.  For example, if memory was
     *    allocated in mdlStart, this is the place to free it.
     */
    static void mdlTerminate(SimStruct *S)
    {
        Bs* spline = (Bs*) ssGetPWork(S)[0];
        delete spline;
        Eigen::RowVectorXd* der = (Eigen::RowVectorXd*) ssGetPWork(S)[1];
        delete der;
    }                                             
    /*======================================================*
     * See sfuntmpl_doc.c for the optional S-function methods *
     *======================================================*/
    
    /*=============================*
     * Required S-function trailer *
     *=============================*/
    
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
    
#ifdef __cplusplus
} // end of extern "C" scope
#endif