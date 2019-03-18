/* Octave include*/
#include <octave/oct.h>
#include <iostream>
//#include <octave/ov.h>

/* mSim include */
#include "Convoluteurf.h"

DEFUN_DLD (calcbtc_oct, args, nargout , "BTC = calcbtc_oct(URF,LF); \n \n"){

	
	octave_value_list retval;
	retval.resize (1);

	int nargin = args.length ();
	if (nargin != 2){
		std::cout << "Wrong number of inputs" << std::endl;
		print_usage ();
		retval(0)=0;
		return retval;
	}


	Matrix URF;
	Matrix LF;

	//assign input arguments to matrices
	URF = args(0).matrix_value();
	LF = args(1).matrix_value();
	int Ns = args(1).rows();
	int Nt = args(1).columns();
	Matrix BTC(Ns,Nt);

	clcbtc(BTC.fortran_vec(), URF.fortran_vec(), LF.fortran_vec(), Ns, Nt);
	retval(0)=BTC;

	return octave_value_list(retval);

}