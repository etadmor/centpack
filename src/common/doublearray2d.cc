////////////////////////////////////////////////////////////////////////////////
//
//  doublearray2d.cc -- implementation of class doublearray2d
//                      build from the array classes provided by C. Anderson
//                      (UCLA)
//
// data structures used by Central Hyperbolic Solver (C) Jorge Balbas and
// Eitan Tadmor, Nov. 2004
//
////////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdio>
#include"doublearray2d.h"

using namespace std;

doublearray2d::doublearray2d()
{
	dataPtr = 0;
	internalAlloc = 0;
	index1Size = 0;
	index1Begin = 0;
	index1End = 0;
	index2Size = 0;
	index2Begin = 0;
	index2End = 0;
}

doublearray2d::doublearray2d(long size1, long size2)
{
	dataPtr = 0;
	internalAlloc = 0; 
	initialize(size1, size2);
}

doublearray2d::doublearray2d(const doublearray2d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	index2Size = d.index2Size;
	index2Begin = d.index2Begin;
	index2End = d.index2End;
	
	dataPtr = new double[index1Size*index2Size];
	internalAlloc = 1; 
	
	long i;
	for(i = 0; i < index1Size*index2Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray2d::~doublearray2d()
{
	if(internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray2d::initialize(long size1, long size2)
{
	if(internalAlloc == 1) 
	{
		if((index1Size != size1) || (index2Size != size2))
		{
			delete [] dataPtr;
			dataPtr = new double[size1*size2];
		}
	}
	
	else
	{
		if(dataPtr == 0)
		{
			dataPtr = new double[size1*size2];
			internalAlloc  = 1;
		}
	}
	
	index1Size = size1;
	index1Begin = 0;
	index1End = index1Begin + (index1Size - 1);
	index2Size = size2;
	index2Begin = 0;
	index2End = index2Begin + (index2Size - 1);
	
}

void doublearray2d::initialize(const doublearray2d& d)
{
    if(internalAlloc == 1) 
    {
		if((index1Size != d.index1Size) ||(index2Size != d.index2Size))
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size*d.index2Size];
		}
    }
    
	else 
    {
		if(dataPtr == 0)
		{
			dataPtr = new double[d.index1Size*d.index2Size];
			internalAlloc = 1;
		}
    }
    
    index1Size = d.index1Size;
    index1Begin = d.index1Begin;
    index1End = d.index1End;
    index2Size = d.index2Size;
    index2Begin = d.index2Begin;
    index2End = d.index2End;
    
    long i;
    
    for(i = 0; i < index1Size*index2Size; i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}

double& doublearray2d::operator()(long i1, long i2)
{
	return *(dataPtr + (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
}

const double& doublearray2d::operator()(long i1, long i2) const
{
	return *(dataPtr + (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
}

double* doublearray2d::getDataPointer()
{
	return dataPtr;
}

void doublearray2d::setIndex1Begin(long i)
{
	index1Begin = i;
	index1End = index1Begin + (index1Size - 1);
}

long doublearray2d::getIndex1Begin() const
{
	return index1Begin;
}

long doublearray2d::getIndex1End() const
{
	return index1End;
}

long doublearray2d::getIndex1Size() const
{
	return index1Size;
}

void doublearray2d::setIndex2Begin(long i)
{
	index2Begin = i;
	index2End = index2Begin + (index2Size - 1);
}

long doublearray2d::getIndex2Begin() const
{
	return index2Begin;
}

long doublearray2d::getIndex2End() const
{
	return index2End;
}

long doublearray2d::getIndex2Size() const
{
	return index2Size;
}

void doublearray2d::operator=(const doublearray2d& d)
{
    if(index1Size*index2Size == 0)
	    initialize(d.index1Size,d.index2Size);
	
    long i;
    
    for(i = 0; i < d.index1Size*d.index2Size; i++)
	    dataPtr[i] = d.dataPtr[i];
}

void doublearray2d::setToValue(double val)
{
	long i;
	for(i = 0; i < index1Size*index2Size; i++)
	{
		dataPtr[i] =  val;
	}
}


void doublearray2d::addValue(double val)
{
	long i;
	
	for(i = 0; i < index1Size*index2Size; i++)
	{
		dataPtr[i] += val;
	}
}

ostream& operator<< (ostream& out_stream, const doublearray2d& d)
{
	double outvalue;

	cout.setf(ios::scientific, ios::floatfield);
	cout.precision(16);
	
	long i,j,k;
	for(i = d.index1Begin; i <= d.index1End; i++)
	{
		for(j = d.index2Begin; j <= d.index2End; j++)
		{
			k=(j - d.index2Begin) + (i - d.index1Begin)*d.index2Size;
			outvalue=d.dataPtr[k] ;
		
			if(outvalue < 0 )
				out_stream << setprecision(16)<<outvalue <<" ";
			else
				out_stream <<" "<< setprecision(16)<<outvalue <<" ";
		}
	
		out_stream<<endl;
	}
	
	return out_stream;
}
