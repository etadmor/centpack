////////////////////////////////////////////////////////////////////////////////
//
//  doublearray2d.h -- headers of class doublearray2d
//                     build from the array classes provided by C. Anderson
//                     (UCLA)
//
// data structures used by Central Hyperbolic Solver (C) Jorge Balbas and
// Eitan Tadmor, Nov. 2004
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray2d__
#define __doublearray2d__

#include<iostream>

using namespace std;

class doublearray2d
{
	friend ostream& operator << (ostream& outstream, const doublearray2d& d);	
	
	private:
	
	double* dataPtr;
	long index1Begin;
	long index1End;
	long index1Size;
	long index2Begin;
	long index2End;
	long index2Size;
	int internalAlloc;
	
	public:
	
	doublearray2d();
	doublearray2d(long size1, long size2);
	doublearray2d(const doublearray2d& d);
	~doublearray2d();
	void initialize(long size1, long size2);
	void initialize(const doublearray2d& d);
	double&  operator()(long i1, long i2);
	const double& operator()(long i1, long i2) const;
	double* getDataPointer();
	
	void setIndex1Begin(long i);
	long getIndex1Begin() const;
	long getIndex1End() const;
	long getIndex1Size() const;
	
	void setIndex2Begin(long i);
	long getIndex2Begin() const;
	long getIndex2End() const;
	long getIndex2Size() const;
	
	void operator=(const doublearray2d& d);
	
	void setToValue(double d);
	void addValue(double d);
};
#endif
