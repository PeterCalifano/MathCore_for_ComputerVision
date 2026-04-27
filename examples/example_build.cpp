#include <iostream>

#include <mathcore/mathcore.h>

int main()
{
    const Eigen::Vector3d axis{1.0, 0.0, 0.0};
    const auto skew = mathcore::skew_symmetric_matrix(axis);
    const auto quat = mathcore::dcmToQuaternion<Eigen::Vector4d>(Eigen::Matrix3d::Identity());

    std::cout << "Skew matrix:\n" << skew << '\n';
    std::cout << "Identity quaternion: " << quat.transpose() << '\n';
    return 0;
}
