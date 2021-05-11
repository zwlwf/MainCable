#include <string>
#include <cstdio>
#include <sstream>
std::string num2str(int i);
std::string num2str( double i );
std::string num2str( double i, char const * fmt );
std::string removeComments(const std::string &);
std::string configStream( const std::string filename);

#include <ctime>
#include <iostream>

class mytime
{
	private:
	clock_t tStart;
	clock_t tEnd;

	public:
	mytime() { tStart = std::clock(); }

	void tic()
	{   tStart= std::clock(); }

	void toc()
	{
		tEnd = std::clock();
		clock_t tspan =tEnd-tStart ;
		if ( tspan <=10*CLOCKS_PER_SEC)
			std::cout<<1000.0*tspan/CLOCKS_PER_SEC<<" ms escaped!\n";
		else if ( tspan <=3600*CLOCKS_PER_SEC)
			std::cout<<tspan/CLOCKS_PER_SEC<<" s escaped!\n";
		else 
			std::cout<<tspan/CLOCKS_PER_SEC/60<<" mins escaped!\n";
	}
}; 
