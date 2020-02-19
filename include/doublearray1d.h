////////////////////////////////////////////////////////////////////////////////
//
//  doublearray1d.h -- headers of class doublearray1d
//                     build from the array classes provided by C. Anderson
//                     (UCLA)
//
// data structures used by Central Hyperbolic Solver (C) Jorge Balbas and
// Eitan Tadmor, Nov. 2004
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray1d__
#define __doublearray1d__

#include<iostream>

using namespace std;

class doublearray1d
{
	friend ostream& operator << (ostream& outstream, const doublearray1d& d);	
	
	private:
	
	double* dataPtr;
	long    index1Begin;
	long    index1End;
	long    index1Size;
	int     internalAlloc;
	
	public:
	
	doublearray1d();
	doublearray1d(long size);
	doublearray1d(const doublearray1d& d);
	~doublearray1d();
	void initialize(long m);
	void initialize(const doublearray1d& d);
	double&  operator()(long i1);
	const double&  operator()(long i1) const;
	double* getDataPointer();
	
	void setIndex1Begin(long i);
	long getIndex1Begin() const;
	long getIndex1End() const;
	long getIndex1Size() const;
    
	void resize(long newSize);
	void operator=(const doublearray1d& d);
	
	void setToValue(double d);
	void addValue(double d);
};
#endif
