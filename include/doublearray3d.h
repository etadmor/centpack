////////////////////////////////////////////////////////////////////////////////
//
//  doublearray3d.h -- headers of class doublearray3d
//                     build from the array classes provided by C. Anderson
//                     (UCLA)
//
// data structures used by Central Hyperbolic Solver (C) Jorge Balbas and
// Eitan Tadmor, Nov. 2004
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray3d__
#define __doublearray3d__

#include<iostream>

using namespace std;

class doublearray3d
{
	friend ostream& operator << (ostream& outstream, const doublearray3d& d);	
	
	private:
	
	double* dataPtr;
	long    index1Begin;
	long    index1End;
	long    index1Size;
	long    index2Begin;
	long    index2End;
	long    index2Size;
	long    index3Begin;
	long    index3End;
	long    index3Size;
	int     internalAlloc;
	
	public:
	
	doublearray3d();
	doublearray3d(long size1, long size2, long size3);
	doublearray3d(const doublearray3d& d);
	~doublearray3d();
	void initialize(long size1, long size2, long size3);
	void initialize(const doublearray3d& d);
	double&  operator()(long i1, long i2, long i3);
	const double& operator()(long i1, long i2, long i3) const;
	double* getDataPointer();
	
	void setIndex1Begin(long i);
	long getIndex1Begin() const;
	long getIndex1End() const;
	long getIndex1Size() const;
	
	void setIndex2Begin(long i);
	long getIndex2Begin() const;
	long getIndex2End() const;
	long getIndex2Size() const;
	
	void setIndex3Begin(long i);
	long getIndex3Begin() const;
	long getIndex3End() const;
	long getIndex3Size() const;
	
	void operator=(const doublearray3d& d);
	
	void setToValue(double d);
	void addValue(double d);
};

#endif
