#include "example_project.h"

int main()
{
    std::cout << "Consumer project using MathCore_for_SpaceNav through find_package().\n";

    const Eigen::Vector3d axis{1.0, 2.0, 3.0};
    const auto skew = mathcore::skew_symmetric_matrix(axis);

    std::cout << "Skew matrix:\n" << skew << '\n';

    return 0;
}
