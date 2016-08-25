% to be built with VS2010 as VS2015 is too strict
delete bst.mex*
mex bst.cpp Bs.cpp Bs_coeffs.cpp
delete bst_sl.mex*
mex bst_sl.cpp Bs.cpp Bs_coeffs.cpp