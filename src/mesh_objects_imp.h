//#include "mesh_objects.hpp"
#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

template <UInt NNODES>
void Element<NNODES, 2, 2>::computeProperties()
{
        /*
	* 		     _ *3
	*                    /|\
	* 	        d2  /   \
	* 		   /	 \
	* 		  /	  \
	* 		 *_________\*
	* 		1    d1   /  2
 	*/
	Element<NNODES, 2, 2> & t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1]);   //Point d2 = t[2] - t[0]; reimplementare sottrazione

	// M_J = [d1, d2]
	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)

	// det_J = det(M_J)
	detJ_ = M_J_(0,0) * M_J_(1,1) - M_J_(1,0) * M_J_(0,1);

	// idet - inverse det_J
	Real idet = 1. / detJ_;

	// M_invJ = M_J^(-1)
	M_invJ_(0,0) =  idet * M_J_(1,1);	// idet * [  (y3-y1)	 -(x3-x1) ]
	M_invJ_(1,0) = -idet * M_J_(1,0);	//        [ -(y2-y1)    (x2-x1) ]
	M_invJ_(0,1) = -idet * M_J_(0,1);	//  =  M_invJ
	M_invJ_(1,1) =  idet * M_J_(0,0);	//	è la trasposta di quella della Sangalli (Ael)

	// Metric matrix
	metric_ = M_invJ_ * M_invJ_.transpose();
}

template <UInt NNODES>
Eigen::Matrix<Real, 3, 1> Element<NNODES, 2, 2>::getBaryCoordinates(const Point & point) const
{
	Element<NNODES, 2, 2>     t = *this;
	Eigen::Matrix<Real, 3, 1> lambda;
	Eigen::Matrix<Real, 4, 1> bary_coef;
//Real                      eps = 2.2204e-016,
//Real                      tolerance = 10000 * eps;

	//cout << "primovert  " << t[0](0) << endl;
	//cout << "primovert  " << t[0](1) << endl;

	bary_coef[0] = t[0][0]-t[2][0];  //x1-x3
	bary_coef[1] = t[1][0]-t[2][0];  //x2-x3
	bary_coef[2] = t[0][1]-t[2][1];  //y1-y3
	bary_coef[3] = t[1][1]-t[2][1];  //y2-y3
	//cout << baryccoef << endl;

	Real detT = bary_coef[0]*bary_coef[3] - bary_coef[1]*bary_coef[2];
	bary_coef = bary_coef / detT;
	//cout << "detT  " << detT << endl;

	//Compute barycentric coordinates for the point
	Real x_diff_third = point[0] - t[2][0];
	Real y_diff_third = point[1] - t[2][1];

	lambda[0] = ( bary_coef[3]*x_diff_third - bary_coef[1]*y_diff_third),
	lambda[1] = (-bary_coef[2]*x_diff_third + bary_coef[0]*y_diff_third),
	lambda[2] = 1 - lambda[0] - lambda[1];

	return lambda;
}

template <UInt NNODES>
bool Element<NNODES, 2, 2>::isPointInside(const Point & point) const
{
	Real eps       = 2.2204e-016;
	Real tolerance = 10 * eps;

	Eigen::Matrix<Real, 3, 1> lambda = getBaryCoordinates(point);

	return ((-tolerance <= lambda[0]) &&
	        (-tolerance <= lambda[1]) &&
		(-tolerance <= lambda[2]) );
}

// TO BE FIXED: if one dir -1, try with others
template <UInt NNODES>
int Element<NNODES, 2, 2>::getPointDirection(const Point & point) const
{
	Real eps       = 2.2204e-016;
	Real tolerance = 10 * eps;

	Eigen::Matrix<Real, 3, 1> lambda = getBaryCoordinates(point);

	//Find the minimum coordinate (if negative stronger straight to the point searched)
	int min_index;
	if (lambda[0] <= lambda[1] && lambda[0] <= lambda[2])
		min_index = 0;
	else if (lambda[1] <= lambda[0] && lambda[1] <= lambda[2])
		min_index = 1;
	else
		min_index = 2;

	if(lambda[min_index] < -tolerance)
		return min_index;
	else
		return -1;
}

template <UInt NNODES>
void Element<NNODES, 2, 2>::print(std::ostream & out) const
{
	out << "Triangle - id " << id_ << "- nodes: ";
	for (UInt i = 0; i < NNODES; ++i)
		out << points_[i].getId() << "  ";
	out << std::endl;
}

//----------------------------------------------------------------------------//
//Implementation myDim = 2, nDim = 3
template <UInt NNODES>
void Element<NNODES, 2, 3>::computeProperties()
{
        /*
	* 		     _ *3
	*                    /|\
	* 	        d2  /   \
	* 		   /	 \
	* 		  /	  \
	* 		 *_________\*
	* 		1    d1   /  2
 	*/
	Element<NNODES, 2, 3> &t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]);

	// M_J = [d1, d2]
	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(2,0) = d1[2];			// (z2-z1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)
	M_J_(2,1) = d2[2];			// (z3-z1)

	G_J_=M_J_.transpose()*M_J_;

	// det_J = det(G_J)
	detJ_ = G_J_(0,0) * G_J_(1,1) - G_J_(1,0) * G_J_(0,1);

	// det_J = det(M_J)
	Real idet = 1. / detJ_;

	// NOT NEEED
	//M_invJ_(0,0) =  idet * M_J_(1,1);	// (y3-y1)	-(x3-x1)
	//M_invJ_(1,0) = -idet * M_J_(1,0);	// -(y2-y1) (x2-x1)
	//M_invJ_(0,1) = -idet * M_J_(0,1);	//
	//M_invJ_(1,1) =  idet * M_J_(0,0);	//	è la trasposta di quella della Sangalli (Ael)

	// Metric matrix
	metric_(0,0) =  idet * G_J_(1,1);
	metric_(1,0) = -idet * G_J_(1,0);
	metric_(0,1) = -idet * G_J_(0,1);
	metric_(1,1) =  idet * G_J_(0,0);
}

template <UInt NNODES>
Eigen::Matrix<Real, 3, 1> Element<NNODES, 2, 3>::getBaryCoordinates(const Point & point) const
{
	Element<NNODES, 2, 3>     t = *this;
	Eigen::Matrix<Real, 3, 1> lambda;
	Real                      detJ_point;
	Eigen::Matrix<Real, 3, 2> M_J_point;
	Eigen::Matrix<Real, 2, 2> G_J_point;

	for (int k = 0; k < 2; ++k)
	{
		Point d1(t[1][0]-point[0], t[1][1]-point[1], t[1][2]-point[2]);
		Point d2(t[2][0]-point[0], t[2][1]-point[1], t[2][2]-point[2]);

		M_J_point(0,0) = d1[0];			// (x2-x1)
		M_J_point(1,0) = d1[1];			// (y2-y1)
		M_J_point(2,0) = d1[2];			// (z2-z1)
		M_J_point(0,1) = d2[0];			// (x3-x1)
		M_J_point(1,1) = d2[1];			// (y3-y1)
		M_J_point(2,1) = d2[2];			// (z3-z1)

		G_J_point = M_J_point.transpose()*M_J_point;

		detJ_point = G_J_point(0,0)*G_J_point(1,1) - G_J_point(1,0)*G_J_point(0,1);

		lambda[k] = std::sqrt(detJ_point)/t.getArea();
	}

	// Partition of uniy property
	lambda[2] = 1-lambda[0]-lambda[1];

	return lambda;
}

// We solve 3 scalar equation in 2 unknowns(u,v)
// u*(P1-P0)+v*(P2-P0)=P-P0
// if the system ins solveable, P is in the plane (P1,P2,P0), if in addition
// u,v>=0 and u+v<=1 then P is inside the triangle
template <UInt NNODES>
bool Element<NNODES, 2, 3>::isPointInside(const Point & point) const
{
	//Real eps = 2.2204e-016;
	//tolerance = 10 * eps;

	// First: check consistency trough Rouchè-Capelli theorem

	Element<NNODES, 2, 3> t = *this;

	Eigen::Matrix<Real, 3, 2> A;
	Eigen::Matrix<Real, 3, 1> b;
	Eigen::Matrix<Real, 2, 1> sol;
	Eigen::Matrix<Real, 3, 1> err;

	A(0,0) = t[1][0]-t[0][0];
	A(0,1) = t[2][0]-t[0][0];
	A(1,0) = t[1][1]-t[0][1];
	A(1,1) = t[2][1]-t[0][1];
	A(2,0) = t[1][2]-t[0][2];
	A(2,1) = t[2][2]-t[0][2];

	b(0)   = point[0]-t[0][0];
	b(1)   = point[1]-t[0][1];
	b(2)   = point[2]-t[0][2];

	sol = A.colPivHouseholderQr().solve(b);
	err = A*sol-b;

	Real tolerance = (A(0,0)*A(0,0) + A(1,0)*A(1,0) + A(2,0)*A(2,0) +
	 		  A(0,1)*A(0,1) + A(1,1)*A(1,1) + A(2,1)*A(2,1))/4;

	if((err(0)*err(0) + err(1)*err(1) + err(2)*err(2)) < tolerance)
	{
		return ((sol(0)+sol(1)<=1) && (sol(0)>=0) && (sol(1)>=0));
	}
	else
	{
		return 0;
	}
}

/*
template <UInt NNODES>
int Triangle<NNODES,2,3>::getPointDirection(const Point& point) const
{

	//da implementare
	std::cerr<<"ancora da implementare";
}
*/

template <UInt NNODES>
void Element<NNODES, 2, 3>::print(std::ostream & out) const
{
	out << "Triangle - id " << id_ << "- nodes: ";
	for (UInt i = 0; i < NNODES; ++i)
		out << points_[i].getId() << "  ";
	out << std::endl;
}

//----------------------------------------------------------------------------//
//Implementation myDim = 3, nDim = 3

template <UInt NNODES>
void Element<NNODES,3,3>::computeProperties()
{
	Element<NNODES, 3, 3> &t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]);
	Point d3(t[3][0]-t[0][0], t[3][1]-t[0][1], t[3][2]-t[0][2]);

	// M_J = [d1, d2, d3]
	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(2,0) = d1[2];			// (z2-z1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)
	M_J_(2,1) = d2[2];			// (z3-z1)
	M_J_(0,2) = d3[0];			// (x4-x1)
	M_J_(1,2) = d3[1];			// (y4-y1)
	M_J_(2,2) = d3[2];			// (z4-z1)

	// detMJ_ = det(M_J_)
	Real detMJ_ = 	M_J_(0,0) * (M_J_(1,1) * M_J_(2,2) - M_J_(1,2) * M_J_(2,1)) -
			M_J_(0,1) * (M_J_(1,0) * M_J_(2,2) - M_J_(1,2) * M_J_(2,0)) +
			M_J_(0,2) * (M_J_(1,0) * M_J_(2,1) - M_J_(1,1) * M_J_(2,0));

	// idetMJ = inverse det(M_J_)
	Real idetMJ = 1. / detMJ_;

	// M_invJ = M_J^(-1)
	M_invJ_(0,0) =  idetMJ * (M_J_(1, 1) * M_J_(2, 2) - M_J_(1, 2) * M_J_(2, 1));
	M_invJ_(0,1) =  idetMJ * (M_J_(0, 2) * M_J_(2, 1) - M_J_(0, 1) * M_J_(2, 2));
	M_invJ_(0,2) =  idetMJ * (M_J_(0, 0) * M_J_(2, 2) - M_J_(0, 2) * M_J_(2, 0));
	M_invJ_(1,0) =  idetMJ * (M_J_(1, 2) * M_J_(2, 0) - M_J_(1, 0) * M_J_(2, 2));
	M_invJ_(1,1) =  idetMJ * (M_J_(0, 0) * M_J_(2, 2) - M_J_(0, 2) * M_J_(2, 0));
	M_invJ_(1,2) =  idetMJ * (M_J_(1, 0) * M_J_(0, 2) - M_J_(0, 0) * M_J_(1, 2));
	M_invJ_(2,0) =  idetMJ * (M_J_(1, 0) * M_J_(2, 1) - M_J_(2, 0) * M_J_(1, 1));
	M_invJ_(2,1) =  idetMJ * (M_J_(2, 0) * M_J_(0, 1) - M_J_(0, 0) * M_J_(2, 1));
	M_invJ_(2,2) =  idetMJ * (M_J_(0, 0) * M_J_(1, 1) - M_J_(1, 0) * M_J_(0, 1));

	G_J_ = M_J_.transpose()*M_J_;

	// det_ = det(G_J_)
	detJ_ = G_J_(0,0) * (G_J_(1,1) * G_J_(2,2) - G_J_(1,2) * G_J_(2,1)) -
		G_J_(0,1) * (G_J_(1,0) * G_J_(2,2) - G_J_(1,2) * G_J_(2,0)) +
		G_J_(0,2) * (G_J_(1,0) * G_J_(2,1) - G_J_(1,1) * G_J_(2,0));

	// idet = inverse det(G_J_)
	Real idet = 1. / detJ_;

	// Metric mtix
	metric_(0,0) =  idet * (G_J_(1, 1) * G_J_(2, 2) - G_J_(1, 2) * G_J_(2, 1));
	metric_(0,1) =  idet * (G_J_(0, 2) * G_J_(2, 1) - G_J_(0, 1) * G_J_(2, 2));
	metric_(0,2) =  idet * (G_J_(0, 0) * G_J_(2, 2) - G_J_(0, 2) * G_J_(2, 0));
	metric_(1,0) =  idet * (G_J_(1, 2) * G_J_(2, 0) - G_J_(1, 0) * G_J_(2, 2));
	metric_(1,1) =  idet * (G_J_(0, 0) * G_J_(2, 2) - G_J_(0, 2) * G_J_(2, 0));
	metric_(1,2) =  idet * (G_J_(1, 0) * G_J_(0, 2) - G_J_(0, 0) * G_J_(1, 2));
	metric_(2,0) =  idet * (G_J_(1, 0) * G_J_(2, 1) - G_J_(2, 0) * G_J_(1, 1));
	metric_(2,1) =  idet * (G_J_(2, 0) * G_J_(0, 1) - G_J_(0, 0) * G_J_(2, 1));
	metric_(2,2) =  idet * (G_J_(0, 0) * G_J_(1, 1) - G_J_(1, 0) * G_J_(0, 1));

	// Volume computation
	Eigen::Matrix<Real,4,4> m;
	m(0,0) = 1;
	m(0,1) = 1;
	m(0,2) = 1;
	m(0,3) = 1;
	m(1,0) = t[0][0];
	m(1,1) = t[1][0];
	m(1,2) = t[2][0];
	m(1,4) = t[3][0];
	m(2,0) = t[0][1];
	m(2,1) = t[1][1];
	m(2,2) = t[2][1];
	m(2,3) = t[3][1];
	m(3,0) = t[0][2];
	m(3,1) = t[1][2];
	m(3,2) = t[2][2];
	m(3,3) = t[3][2];

	Volume_= 1./6*
	std::abs(m(0,3) * m(1,2) * m(2,1) * m(3,0) - m(0,2) * m(1,3) * m(2,1) * m(3,0) -
         	 m(0,3) * m(1,1) * m(2,2) * m(3,0) + m(0,1) * m(1,3) * m(2,2) * m(3,0) +
         	 m(0,2) * m(1,1) * m(2,3) * m(3,0) - m(0,1) * m(1,2) * m(2,3) * m(3,0) -
         	 m(0,3) * m(1,2) * m(2,0) * m(3,1) + m(0,2) * m(1,3) * m(2,0) * m(3,1) +
         	 m(0,3) * m(1,0) * m(2,2) * m(3,1) - m(0,0) * m(1,3) * m(2,2) * m(3,1) -
         	 m(0,2) * m(1,0) * m(2,3) * m(3,1) + m(0,0) * m(1,2) * m(2,3) * m(3,1) +
         	 m(0,3) * m(1,1) * m(2,0) * m(3,2) - m(0,1) * m(1,3) * m(2,0) * m(3,2) -
         	 m(0,3) * m(1,0) * m(2,1) * m(3,2) + m(0,0) * m(1,3) * m(2,1) * m(3,2) +
         	 m(0,1) * m(1,0) * m(2,3) * m(3,2) - m(0,0) * m(1,1) * m(2,3) * m(3,2) -
         	 m(0,2) * m(1,1) * m(2,0) * m(3,3) + m(0,1) * m(1,2) * m(2,0) * m(3,3) +
         	 m(0,2) * m(1,0) * m(2,1) * m(3,3) - m(0,0) * m(1,2) * m(2,1) * m(3,3) -
         	 m(0,1) * m(1,0) * m(2,2) * m(3,3) + m(0,0) * m(1,1) * m(2,2) * m(3,3));
}

template <UInt NNODES>
Eigen::Matrix<Real, 4, 1> Element<NNODES, 3, 3>::getBaryCoordinates(const Point & point) const
{
	Element<NNODES, 3, 3>     t = *this;
	Eigen::Matrix<Real, 4, 1> lambda;
	Eigen::Matrix<Real, 3, 3> M_J_point;
	Eigen::Matrix<Real, 3, 1> rhs;
	Eigen::Matrix<Real, 3, 1> sol;

	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1], t[1][2]-t[0][2]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1], t[2][2]-t[0][2]);
	Point d3(t[3][0]-t[0][0], t[3][1]-t[0][1], t[3][2]-t[0][2]);

	M_J_point(0,0) = d1[0];			// (x2-x1)
	M_J_point(1,0) = d1[1];			// (y2-y1)
	M_J_point(2,0) = d1[2];			// (z2-z1)
	M_J_point(0,1) = d2[0];			// (x3-x1)
	M_J_point(1,1) = d2[1];			// (y3-y1)
	M_J_point(2,1) = d2[2];			// (z3-z1)
	M_J_point(0,2) = d3[0];			// (x4-x1)
	M_J_point(1,2) = d3[1];			// (y4-y1)
	M_J_point(2,2) = d3[2];			// (z4-z1)

	rhs(0) = point[0]-t[0][0];
	rhs(1) = point[1]-t[0][1];
	rhs(2) = point[2]-t[0][2];

	sol = M_J_point.colPivHouseholderQr().solve(rhs);

	lambda[1] = sol(0);
	lambda[2] = sol(1);
	lambda[3] = sol(2);

	// Partition of unity property
	lambda[0]=1-lambda[1]-lambda[2]-lambda[3];

	return lambda;
}

template <UInt NNODES>
bool Element<NNODES, 3, 3>::isPointInside(const Point & point) const
{
	Real eps       = 2.2204e-016;
	Real tolerance = 10 * eps;

	Element<NNODES, 3, 3> t = *this;
	Eigen::Matrix<Real, 4, 1> bary_coeff = t.getBaryCoordinates(point);
	return (-tolerance <= bary_coeff[0]) && (-tolerance <= bary_coeff[1]) && (-tolerance <= bary_coeff[2]) && (-tolerance <= bary_coeff[3]);
}

template <UInt NNODES>
void Element<NNODES, 3, 3>::print(std::ostream & out) const
{
	out << "Tetrahedron - id " << id_ << "- nodes: ";
	for (UInt i = 0; i < NNODES; ++i)
		out << points_[i].getId() << "  ";
	out << std::endl;
}

#endif
