#include "fitsreader.h"
#include "fitswriter.h"

#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>
#include <string>

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout << "Syntax: applybeam [-not-squared / -is-weight] <inpfits> <beamfits> <outfits>\n";
		return 0;
	}
	
	bool squared = true, isWeight = false;
	int argi = 1;
	while(argi < argc && argv[argi][0] == '-')
	{
		std::string p(&argv[argi][1]);
		if(p == "not-squared")
		{
			squared = false;
		}
		else if(p == "is-weight")
		{
			isWeight = true;
		}
		else throw std::runtime_error("Bad parameter");
		++argi;
	}
	
	const char *inpFits = argv[argi];
	const char *beamFits = argv[argi+1];
	const char *outFits = argv[argi+2];
	
	FitsReader inpReader(inpFits);
	size_t
		width = inpReader.ImageWidth(),
		height = inpReader.ImageHeight();
	FitsReader beamReader(beamFits);
	if(beamReader.ImageWidth() != width || beamReader.ImageHeight() != height)
		throw std::runtime_error("Beam and image do not have same size!");
	
	std::vector<double> inpImage(width*height), beamImage(width*height);
	
	inpReader.Read<double>(&inpImage[0]);
	beamReader.Read<double>(&beamImage[0]);
	
	std::vector<double>::iterator beamIter = beamImage.begin();
	for(std::vector<double>::iterator i=inpImage.begin(); i!=inpImage.end(); ++i)
	{
		if(std::fabs(*beamIter) < 1e-2)
			*i = std::numeric_limits<double>::quiet_NaN();
		else if(isWeight)
			*i /= sqrt(*beamIter);
		else if(squared)
			*i /= *beamIter * *beamIter;
		else
			*i /= *beamIter;
		++beamIter;
	}
	
	FitsWriter writer(inpReader);
	writer.Write<double>(outFits, &inpImage[0]);
}
