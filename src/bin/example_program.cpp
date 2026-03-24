#include <iostream>

#include <mathcore/mathcore.h>

int main()
{
    const Eigen::Vector3d omega{0.1, 0.2, 0.3};
    const auto skew = mathcore::skew_symmetric_matrix(omega);

    std::cout << "Example skew matrix:\n" << skew << '\n';
    return 0;
}
