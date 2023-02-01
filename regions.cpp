#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

#include <eigen3/Eigen/Dense>

void define_cube_as_constraints(std::ofstream &constraints, const int &obs_idx, const double &fx, const double &tx, const double &fy, const double &ty, const double &fz, const double &tz) {
  // initialise output to write in the file
  Eigen::MatrixXd A(6, 3);
  Eigen::VectorXd B(6);

  // define cube vertices
  Eigen::Vector3d p1(tx, fy, tz);
  Eigen::Vector3d q1(fx, fy, tz);
  Eigen::Vector3d r1(fx, fy, fz);

  Eigen::Vector3d p2(tx, fy, tz);
  Eigen::Vector3d q2(tx, ty, tz);
  Eigen::Vector3d r2(tx, ty, fz);

  Eigen::Vector3d p3(fx, ty, tz);
  Eigen::Vector3d q3(tx, ty, tz);
  Eigen::Vector3d r3(fx, ty, fz);

  // compute vectors defining the faces (constraints)
	//H1
	Eigen::Vector3d v1(q1 - p1);
	Eigen::Vector3d v2(r1 - p1);
	Eigen::Vector3d n(v1.cross(v2));
	A(0,0) = n(0); A(0,1) = n(1); A(0,2) = n(2);
	B(0) = n(0)*p1(0)+n(1)*p1(1)+n(2)*p1(2);

	//H2
	v1 = p2 - q2;
	v2 = r2 - q2;
	n = v1.cross(v2);
	A(1,0) = n(0); A(1,1) = n(1); A(1,2) = n(2);
	B(1) = n(0)*p2(0)+n(1)*p2(1)+n(2)*p2(2);

	//H3
	v1 = q2 - p3;
	v2 = r2 - p3;
	n = v1.cross(v2);
	A(2,0) = n(0); A(2,1) = n(1); A(2,2) = n(2);
	B(2) = n(0)*p3(0)+n(1)*p3(1)+n(2)*p3(2);

	//H4
	v1 = q1 - p3;
	v2 = r1 - p3;
	n = v2.cross(v1);
	A(3,0) = n(0); A(3,1) = n(1); A(3,2) = n(2);
	B(3) = n(0)*p3(0)+n(1)*p3(1)+n(2)*p3(2);

	//H5
	v1 = r1 - r3;
	v2 = r2 - r3;
	n = v2.cross(v1);
	A(4,0) = n(0); A(4,1) = n(1); A(4,2) = n(2);
	B(4) = n(0)*r1(0)+n(1)*r1(1)+n(2)*r1(2);

	//H6
	v1 = q1 - p3;
	v2 = q2 - p3;
	n = v1.cross(v2);
	A(5,0) = n(0); A(5,1) = n(1); A(5,2) = n(2);
	B(5) = n(0)*q1(0)+n(1)*q1(1)+n(2)*q1(2);

  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
  constraints << "A" << A.format(CommaInitFmt) << " A_list_.push_back(A);" << "\n";
  constraints << "B" << B.format(CommaInitFmt) << " B_list_.push_back(B);" << "\n";

}

int main(int argc, char** argv) {

  std::cout << "\ngenerating an octomap of concrete blocks scenario" << std::endl;

  // define map, obstacles, ...
  double x_offset = 0.0;
  double y_offset = 10.0;
  double z_offset = 10.0;
  int obs_x_n = 2; double obs_x_size = 20.0; double inc_x = 16.5;
  int obs_y_n = 1; double obs_y_size = 14.0; double inc_y = 18.0;
  int obs_z_n = 1; double obs_z_size = 10.0; double inc_z = 12.0;

  std::ofstream constraints;
  constraints.open("3d_unicycle_regions.txt", std::ios::out | std::ios::trunc);// | std::ios::app | std::ios::binary);
  constraints << "n_regions_ = " << 2 << ";\n\n";

  	define_cube_as_constraints(constraints, 0, 0.0, 90.0, 50.0, 70.0, 0.0, 10.0);
	// // define_cube_as_constraints(constraints, 1, 0.0, 20.0, 20.0, 55.0, 0.0, 10.0);
	// // define_cube_as_constraints(constraints, 2, 90.0, 100.0, 0.0, 55.0, 0.0, 10.0);
	define_cube_as_constraints(constraints, 3, 10.0, 50.0, -50.0, -40.0, 0.0, 10.0);
	// define_cube_as_constraints(constraints, 2, 70.0, 90.0, -50.0, -40.0, 0.0, 10.0);
  	// define_cube_as_constraints(constraints, 0, 0.0, 20.0, 80.0, 100.0, 0.0, 10.0);
  	// define_cube_as_constraints(constraints, 1, 80.0, 100.0, 80.0, 100.0, 0.0, 10.0);

	// define_cube_as_constraints(constraints, 0, 0.0, 45.0, 50.0, 70.0, 0.0, 10.0);
	// define_cube_as_constraints(constraints, 1, 55.0, 100.0, 50.0, 70.0, 0.0, 10.0);
//   define_cube_as_constraints(constraints, 1, 100.0, 110.0, 0.0, 100.0, 0.0, 10.0);
//   define_cube_as_constraints(constraints, 1, 0.0, 100.0, 0.0, -100.0, 0.0, 10.0);
//   define_cube_as_constraints(constraints, 1, -100.0, 0.0, 0.0, 100.0, 0.0, 10.0);
  	// define_cube_as_constraints(constraints, 1, -10.0, 0.0, -60.0, 100.0, 0.0, 10.0);
	// define_cube_as_constraints(constraints, 1, 100.0, 110.0, -60.0, 100.0, 0.0, 10.0);

  constraints.close();
  std::cout << "Wrote constraints" << std::endl << std::endl;
  return 0;
}
