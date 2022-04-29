#include <igl/AABB.h>
#include <igl/in_element.h>
#include <igl/readMESH.h>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        std::cout << "Description:\n";
        std::cout << "Check if TetMesh[i] is contained in TetMesh[i+1] for i = 1...N-1\n";
        std::cout << "Usage:\n";
        std::cout
            << "check-tet-mesh-hierarchy.exe TetMesh1.mesh TetMesh2.mesh [TetMesh3...N].mesh\n";
        return 0;
    }

    std::list<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> hierarchy{};
    for (auto i = 1u; i < argc; ++i)
    {
        std::string const file = argv[i];
        Eigen::MatrixXd V;
        Eigen::MatrixXi T, F;
        if (!igl::readMESH(file, V, T, F))
        {
            std::cout << "Could not read mesh " << i << ": " << file << "\n";
            return 0;
        }

        hierarchy.push_back(std::make_pair(V, T));
    }

    auto begin   = hierarchy.begin();
    auto it      = hierarchy.begin();
    auto next    = ++hierarchy.begin();
    auto end     = hierarchy.end();
    bool success = true;
    for (; next != end; ++it, ++next)
    {
        Eigen::MatrixXd const& Vfine = it->first;

        Eigen::MatrixXd const& Vcoarse = next->first;
        Eigen::MatrixXi const& Tcoarse = next->second;

        igl::AABB<Eigen::MatrixXd, 3> tree;
        tree.init(Vcoarse, Tcoarse);
        Eigen::VectorXi I;
        igl::in_element(Vcoarse, Tcoarse, Vfine, tree, I);

        if ((I.array() == -1).any())
        {
            auto const i = std::distance(begin, it);
            std::cout << "Mesh " << i << " not entirely contained in mesh " << i + 1 << "\n";
            success = false;
        }
    }

    if (success)
        std::cout << "All fine layers contained in coarse layers\n";

    return 0;
}