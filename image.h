#ifndef IMAGE_H
#define IMAGE_H

#include <cmath>
#include <cstring>

#include "uvector.h"

class Image
{
public:
	typedef double* iterator;
	typedef const double* const_iterator;
	
	Image() : _data(), _width(0), _height(0) { }
	Image(size_t width, size_t height);
	Image(size_t width, size_t height, double initialValue);
	
	~Image();
	
	Image(const Image&) = default;
	Image& operator=(const Image&) = default;
	Image& operator=(double value);
	
	Image(Image&& source) = default;
	Image& operator=(Image&& source) = default;
	
	double* data() { return _data.data(); }
	const double* data() const { return _data.data(); }
	
	size_t Width() const { return _width; }
	size_t Height() const { return _height; }
	size_t size() const { return _width * _height; }
	bool empty() const { return _width == 0 || _height == 0; }
	
	iterator begin() { return _data.begin(); }
	const_iterator begin() const { return _data.begin(); }
	
	iterator end() { return _data.end(); }
	const_iterator end() const { return _data.end(); }
	
	const double& operator[](size_t index) const { return _data[index]; }
	double& operator[](size_t index) { return _data[index]; }
	
	
	Image& operator*=(double factor);
	Image& operator*=(const Image& other);
	Image& operator/=(double factor)
	{ return (*this) *= 1.0/factor; }
	
	void reset();
	
	/** Cut-off the borders of an image.
	 * @param outWidth Should be &lt;= inWidth.
	 * @param outHeight Should be &lt;= inHeight.
	 */
	static void Trim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	template<typename T>
	static void TrimBox(T* output, size_t x1, size_t y1, size_t boxWidth, size_t boxHeight, const T* input, size_t inWidth, size_t inHeight);
	
	/** Extend an image with zeros, complement of Trim.
	 * @param outWidth Should be &gt;= inWidth.
	 * @param outHeight Should be &gt;= inHeight.
	 */
	static void Untrim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	static double Median(const double* data, size_t size)
	{
		ao::uvector<double> copy;
		return median_with_copy(data, size, copy);
	}
	
	static double MAD(const double* data, size_t size);
	
	double Sum() const;
	double Average() const;
	
	double Min() const;
	double Max() const;
	
	double StdDevFromMAD() const { return StdDevFromMAD(_data.data(), _data.size()); }
	static double StdDevFromMAD(const double* data, size_t size)
	{
		// norminv(0.75) x MAD
		return 1.48260221850560 * MAD(data, size);
	}
	
	static double RMS(const double* data, size_t size)
	{
		double sum = 0.0;
		for(size_t i=0; i!=size; ++i)
			sum += data[i]*data[i];
		return sqrt(sum/size);
	}
	
	void Negate()
	{
		for(double& d : *this)
			d = -d;
	}
private:
	ao::uvector<double> _data;
	size_t _width, _height;
	
	static double median_with_copy(const double* data, size_t size, ao::uvector<double>& copy);
};

#endif
