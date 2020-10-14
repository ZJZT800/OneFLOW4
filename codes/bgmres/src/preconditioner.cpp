
#include "preconditioner.h"
#include "solution.h"
#include "util.h"
#include "UCom.h"
#include <cmath>

/** ************************************************************************
 * Base constructor  for the Preconditioner class. 
 *
 *
 * @param size The number of grid points used in the approximation.
 * ************************************************************************ */
Preconditioner::Preconditioner(int number)
{
	setN(number);
}

/** ************************************************************************
 *	Copy constructor  for the Preconditioner class. 
 *
 *	@param oldCopy The Preconditioner class member to make a copy of.
 * ************************************************************************ */
Preconditioner::Preconditioner(const Preconditioner& oldCopy)
{
	setN(oldCopy.getN());
}

/** ************************************************************************
 *	Destructor for the Preconditioner class. 
 *  ************************************************************************ */
Preconditioner::~Preconditioner()
{
	;
	//ArrayUtils<double>::deltwotensor(vector);
}

/** ************************************************************************
 * The method to solve the system of equations associated with the
 * preconditioner.
 * 
 * Returns the value Solution class that is the solution to the
 * preconditioned system.  For the moment we just assume that the
 * preconditioner is the identity matrix. (This needs to be changed!)
 *
 * @param vector The Solution or right hand side of the system.
 * @return A Solution class member that is the solution to the
 *         preconditioned system.
 * ************************************************************************ */
Solution Preconditioner::solve2(const Solution& current)
{
	Solution multiplied(current);
	int cId;
	int innerlupe;
	for (cId = 0; cId < getN(); cId++)
	{
		for (innerlupe = 0; innerlupe < Rank.COLNUMBER; innerlupe++)
		{
			multiplied(cId, innerlupe) = current.getEntry(cId, innerlupe);
		}
	}
	return(multiplied);
}