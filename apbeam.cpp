#include "fitsreader.h"
#include "fitswriter.h"
#include "image.h"

#include "units/angle.h"
#include "units/imagecoordinates.h"
#include "units/radeccoord.h"
#include "units/ncpprojection.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <boost/optional.hpp>

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout <<
			"\tSyntax: apbeam [options] <input> <outbeam> <outweight>\n"
			"This tool creates an output file with a simple Westerbork beam for the given input beam.\n"
			"options:\n"
			"\t-frequency <value in MHz>\n";
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
	const char* outBeamFilename = argv[argi+1];
	const char* outWeightFilename = argv[argi+2];
	
	FitsReader reader(inpFilename);
	size_t width = reader.ImageWidth(), height = reader.ImageHeight();
	Image beam(width, height), weight(width, height);

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
	if(reader.ProjectionType() == FitsReader::NCPProjection)
		std::cout << "Image is in deprecated NCP projection.\n";
	std::cout <<
		"Making beam with freq=" << freqMHz << " MHz\n"
		"Pixelscale: " << Angle::ToNiceString(reader.PixelSizeX()) << " x " << Angle::ToNiceString(reader.PixelSizeY()) << '\n' <<
		"Phase centre: " << RaDecCoord::RaDecToString(reader.PhaseCentreRA(), reader.PhaseCentreDec()) << '\n';
	
	// In the equation, angle should be in degrees, but we calculate it in rad so absorp the conversion in beta:
	//beta *= 180.0/M_PI; // it's also cos in degrees so not necessary

	double* pbPtr = beam.data();
	double* wPtr = weight.data();
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, reader.PixelSizeX(), reader.PixelSizeY(), width, height, l, m);
			if(reader.ProjectionType() == FitsReader::SINProjection)
			{
				ImageCoordinates::LMToRaDec(l, m, reader.PhaseCentreRA(), reader.PhaseCentreDec(), ra, dec);
			}
			else {
				NCPProjection::LMToRaDec(l, m, reader.PhaseCentreRA(), reader.PhaseCentreDec(), ra, dec);
			}
			double angle = ImageCoordinates::AngularDistance(ra, dec, reader.PhaseCentreRA(), reader.PhaseCentreDec());
			if(x==0 && y==0)
				std::cout << "Max angle: " << Angle::ToNiceString(angle) << '\n';
			
			double cosTerm = cos(beta*freqMHz*angle);
			*pbPtr = cosTerm*cosTerm*cosTerm*cosTerm*cosTerm*cosTerm;
			*wPtr = (*pbPtr) * (*pbPtr);
			++pbPtr;
			++wPtr;
		}
	}
	FitsWriter writer(reader);
	writer.Write(outBeamFilename, beam.data());
	writer.Write(outWeightFilename, weight.data());
}
