#include "fitsreader.h"
#include "fitswriter.h"
#include "image.h"

#include "units/imagecoordinates.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <boost/optional.hpp>

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout <<
			"Syntax: apbeam [options] <input> <output>\n\n"
			"Will create an output file with a simple Westerbork beam for the given input beam.\n";
		return 0;
	}
	
	boost::optional<double> frequency;
	
	int argi = 1;
	while(argi < argc && argv[argi][0] == '-')
	{
		std::string p(&argv[argi][1]);
		if(p == "frequency")
		{
			++argi;
			frequency = atoi(argv[argi])*1e6;
		}
		else throw std::runtime_error("Bad parameter");
		++argi;
	}
	
	const char* inpFilename = argv[argi];
	const char* outFilename = argv[argi+1];
	
	FitsReader reader(inpFilename);
	size_t width = reader.ImageWidth(), height = reader.ImageHeight();
	Image image(width, height);

	// pb = cos**6(beta*freq(MHz)*angle(degrees))
	// where beta = 0.0629 for f < 500 MHz, and 0.065 for f > 500 MHz
	if(!frequency)
		frequency = reader.Frequency();
	double freqMHz = frequency.get()*1e-6;
	double beta;
	if(freqMHz < 500)
		beta = 0.0629;
	else
		beta = 0.065;
	std::cout << "Making beam with freq=" << freqMHz << " MHz\n";
	
	// In the equation, angle should be in degrees, but we calculate it in rad so absorp the conversion in beta:
	beta *= M_PI/180.0;

	double* pbPtr = image.data();
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, reader.PixelSizeX(), reader.PixelSizeY(), width, height, l, m);
			ImageCoordinates::LMToRaDec(l, m, reader.PhaseCentreRA(), reader.PhaseCentreDec(), ra, dec);
			double angle = ImageCoordinates::AngularDistance(ra, dec, reader.PhaseCentreRA(), reader.PhaseCentreDec());
			
			double cosTerm = cos(beta*freqMHz*angle);
			*pbPtr = cosTerm*cosTerm*cosTerm*cosTerm*cosTerm*cosTerm;
			if(x == width/2)
				std::cout << *pbPtr << ' ';
			++pbPtr;
		}
	}
	FitsWriter writer(reader);
	writer.Write(outFilename, image.data());
}
