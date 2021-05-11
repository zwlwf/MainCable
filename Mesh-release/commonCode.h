#ifndef elemMMatrix_fun
virtual amat::Matrix<double> elemMMatrix(MassMode mmode= LUMPM )  
{
	std::cout<<"elemMMatrix Not implemented yet for elemType "<<typeName()<<std::endl;
}
#endif

#ifndef get6DofDisp_fun

virtual void get6DofDisp( double xi, double* elemU, //input
						double d6[], double pos[], //ouput
					 double xi_y, double xi_z )  // optional
{
	std::cout<<"get6DofDisp() function not implemented yet for elemType "<<typeName()<<std::endl;
}
#endif

#ifndef applyNodeForce_fun
virtual void applyNodeForce(double xi, double f_in[], double *f_out) // f_in is with size 6x1
{
	std::cout<<"applyNodeForce() function not implemented yet for elemType "<<typeName()<<std::endl;
}
#endif

#ifndef calcVolume_fun
virtual void calcVolume()
{
	std::cout<<"calcVolume() function not implemented yet for elemType "<<typeName()<<std::endl;
}
#endif

#ifndef elementGlobalForce_fun
virtual amat::Matrix<double> elementGlobalForce( amat::Matrix<double> u )
{
	std::cout<<"elementGlobalForce Not implemented for yet elemType "<<typeName()<<std::endl;
}
#endif

#ifndef elementForce_fun
virtual amat::Matrix<double> elementForce( amat::Matrix<double> u )
{
	std::cout<<"elementForce Not implemented yet for elemType "<<typeName()<<std::endl;
}
#endif
