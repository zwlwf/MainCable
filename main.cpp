#include "MainCable.h"

int main(int argc, char** argv)
{
	MainCable mc;
	if ( argc<2 )
		mc.initialize("config");
	else
		mc.initialize(argv[1]);

	/*
	// check the H0 vs y(mid)
	double H0;
	char outName[512];
	int ind=1;
	for(H0=4.0e8; H0>=4.0e7; H0-=4.0e7)
	{
		mc.solveYZ(H0);
		amat::Matrix<double> & yy = mc.getSolutionY();
		std::cout<<"H0=" <<H0<<" "<<mc.middleSpanY()<<std::endl;
		sprintf(outName,"out_%d.vtk",ind++);
		mc.printStatus( outName );
	}
	*/

	mc.guessH0( );
	mc.solveH0();
	mc.printStatus("out.vtk");
	mc.printAPDL("node_modified.mac");

	/*
	MainCable mc_back;
	mc_back.initialize("config_back");
	mc_back.guessH0( mc.getH0() );
	mc_back.solveYZ();
	mc_back.printAPDL("node_modified2.mac");
	*/
	return 0;
}
