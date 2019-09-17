#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

// HEADERS
#include "fdaPDE.h"
#include "mesh_objects.h"

// CLASSES
class IntegratorTriangleP2
{
	public:
		static const UInt               ORDER  = 1;
		//Number of nodes required
		static const UInt               NNODES = 3;
		//Point locations of the nodes
		static const std::vector<Point> NODES;
		static const std::vector<Real>  WEIGHTS;
};

class IntegratorTriangleP4
{
	public:
		static const UInt               ORDER  = 2;
		//Number of nodes required
		static const UInt 		NNODES = 6;
		//Point locations of the nodes
		static const std::vector<Point> NODES;
		static const std::vector<Real> 	WEIGHTS;
	};

class IntegratorTetrahedronP2
{
	public:
		static const UInt 		ORDER  = 1;
		//Number of nodes required
		static const UInt 		NNODES = 4;
		//Point locations of the nodes
		static const std::vector<Point> NODES;
		static const std::vector<Real> 	WEIGHTS;
};

class IntegratorTetrahedronP1
{
	public:
		static const UInt 		ORDER  = 1;
		//Number of nodes required
		static const UInt 		NNODES = 1;
		//Point locations of the nodes
		static const std::vector<Point> NODES;
		static const std::vector<Real> 	WEIGHTS;
};

//#include "integration_imp.hpp"
#endif
