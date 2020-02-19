////////////////////////////////////////////////////////////////////////////////
//
//  doublearray1d.cc -- implementation of class doublearray1d
//                      build from the array classes provided by C. Anderson
//                      (UCLA)
//
// data structures used by CentPack (C) 2006 Jorge Balbas and
// Eitan Tadmor
//
////////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<iomanip>
#include <cstdio>
#include"doublearray1d.h"

using namespace std;

doublearray1d::doublearray1d()
{
	dataPtr = 0;
	internalAlloc = 0;
	index1Size = 0;
	index1Begin = 0;
	index1End = 0;
}

doublearray1d::doublearray1d(long size)
{
	dataPtr = 0;
	internalAlloc = 0; 
	initialize(size);
}

doublearray1d::doublearray1d(const doublearray1d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	dataPtr = new double[index1Size];
	internalAlloc = 1; 
	
	long i;
	for(i = 0; i < index1Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray1d::~doublearray1d()
{
	if(internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray1d::initialize(long m)
{
	if(internalAlloc == 1) 
    {
		if(index1Size != m)
		{
			delete [] dataPtr;
			dataPtr = new double[m];
		}
    }
	
    else 
    {
		if(dataPtr == 0)
		{
			dataPtr = new double[m];
			internalAlloc = 1;
		}
    }
    
    index1Size = m;
    index1Begin = 0;
    index1End = index1Begin + (index1Size - 1);
}

void doublearray1d::initialize(const doublearray1d& d)
{
    if(internalAlloc == 1) 
    {
		if(index1Size != d.index1Size)
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size];
		}
    }
    
	else 
    {
		if(dataPtr == 0)
		{
			dataPtr = new double[d.index1Size];
			internalAlloc = 1;
		}
    }
    
    index1Size = d.index1Size;
    index1Begin = d.index1Begin;
    index1End = d.index1End;
    
    long i;
    
    for(i = 0; i < index1Size; i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}

double& doublearray1d::operator()(long i1)
{
	return *(dataPtr +  (i1 - index1Begin));
}

const double& doublearray1d::operator()(long i1) const
{
	return *(dataPtr +  (i1 - index1Begin));
}

double* doublearray1d::getDataPointer()
{
	return dataPtr;
}

void doublearray1d::setIndex1Begin(long i)
{
	index1Begin = i;
	index1End   = index1Begin + (index1Size - 1);
}

long doublearray1d::getIndex1Begin() const
{
	return index1Begin;
}

long doublearray1d::getIndex1End() const
{
	return index1End;
}

long doublearray1d::getIndex1Size() const
{
	return index1Size;
}

void doublearray1d::resize(long newSize)
{
	long i;
	double*  newDataPtr = new double[newSize];
	double*  tmpDataPtr;

	if(newSize > index1Size) 
	{
		for(i = 0; i < index1Size; i++)
			newDataPtr[i] = dataPtr[i];
	}
	
	else
	{
		for(i = 0; i < newSize; i++)
			newDataPtr[i] = dataPtr[i];
	}

	index1Size = newSize;
	tmpDataPtr = dataPtr;
	dataPtr    = newDataPtr;

	if(internalAlloc == 1) delete [] tmpDataPtr;
	internalAlloc = 1;

	index1End = index1Begin + (index1Size - 1);
}

void doublearray1d::operator=(const doublearray1d& d)
{
    initialize(d.index1Size);
    
    long i;
    for(i = 0; i < d.index1Size; i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}

void doublearray1d::setToValue(double val)
{
	long i;
	for(i = 0; i < index1Size; i++)
	{
		dataPtr[i] =  val;
	}
}


void doublearray1d::addValue(double val)
{
    long i;
    
    for(i = 0; i < index1Size; i++)
    {
	    dataPtr[i] += val;
    }
}

ostream& operator << (ostream& out_stream, const doublearray1d& d)
{
	double outvalue;

	cout.setf(ios::scientific);
	cout.setf(ios::floatfield);
	cout.precision(16);
	
	long i;
	for(i = d.index1Begin; i <= d.index1End; i++)
	{
		outvalue=d.dataPtr[i];
		
		if(outvalue < 0 )
			out_stream << setprecision(16)<<outvalue <<" ";
		else
			out_stream << " " << setprecision(16)<<outvalue <<" ";
	}
	
	return out_stream;
}
