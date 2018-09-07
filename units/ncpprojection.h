#ifndef NCP_PROJECTION_H
#define NCP_PROJECTION_H

#include <algorithm>
#include <cmath>
#include <vector>

/**
 * This class performs LM coordinate transforms for the (deprecated) NCP projection.
 * As defined in AIPS memo 27 (citing Brouw 1971):
 *   L = cos dec sin deltaRa
 *   M = (cos dec0 - cos dec cos deltaRa) / sin dec0
 * (with deltaRa = ra - ra0)
 */
class NCPProjection
{
public:
	template<typename T>
	static void RaDecToLM(T ra, T dec, T phaseCentreRa, T phaseCentreDec, T &destL, T &destM)
	{
		const T
			deltaAlpha = ra - phaseCentreRa,
			sinDeltaAlpha = sin(deltaAlpha),
			cosDeltaAlpha = cos(deltaAlpha),
			cosDec = cos(dec),
			sinDec0 = sin(phaseCentreDec),
			cosDec0 = cos(phaseCentreDec);
		
		destL = cosDec * sinDeltaAlpha;
		destM = (cosDec0 - cosDec*cosDeltaAlpha) / sinDec0;
	}
	
	template<typename T>
	static void LMToRaDec(T l, T m, T phaseCentreRa, T phaseCentreDec, T &destRa, T &destDec)
	{
		const T
			cosDec0 = cos(phaseCentreDec),
			sinDec0 = sin(phaseCentreDec),
			deltaAlpha = atan2(l, cosDec0 - m*sinDec0),
			cosDeltaAlpha = cos(deltaAlpha);
			
		destRa = phaseCentreRa + deltaAlpha;
		if(phaseCentreDec < 0.0)
			destDec = -acos((cosDec0 - m*sinDec0) / cosDeltaAlpha);
		else
			destDec = acos((cosDec0 - m*sinDec0) / cosDeltaAlpha);
	}

	NCPProjection() = delete;
};

#endif

