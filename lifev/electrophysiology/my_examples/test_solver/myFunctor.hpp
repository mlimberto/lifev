#ifndef _MYFUNCTOR_HPP_
#define _MYFUNCTOR_HPP_

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

template <typename Return_Type > 
class myFunctor
{
public:
	typedef Return_Type return_Type; 

	return_Type operator() (const VectorSmall<3> spaceCoord)
	{
		Real x = spaceCoord[0];
		Real y = spaceCoord[1];
		Real z = spaceCoord[2];
	
		return myFunction(0,x,y,z,0);
	}

	myFunctor( return_Type (*newFunction)(const Real& , const Real &x, const Real &y, const Real &z, const ID&))
			: myFunction(newFunction) {}

	~myFunctor() {}

private:

	return_Type (*myFunction) (const Real& , const Real &x, const Real &y, const Real &z, const ID&);

};


} 

#endif