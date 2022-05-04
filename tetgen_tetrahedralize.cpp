#include <filesystem>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/readSTL.h>
#include <igl/writeMESH.h>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc == 1 || argc != 3)
    {
        std::cout << "Usage:\n";
        std::cout << "./tetgen_tetrahedralize.exe input.[ply|obj|stl] output.mesh\n";
        return 0;
    }

    std::filesystem::path const input{argv[1]};
    std::string const input_extension = input.extension().string();
    bool const is_valid_input_file =
        std::filesystem::exists(input) && std::filesystem::is_regular_file(input) &&
        (input_extension == ".ply" || input_extension == ".obj" || input_extension == ".stl");

    if (!is_valid_input_file)
    {
        std::cerr << "Expected input to be .[ply|obj|stl]\n";
        return 1;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    if (input_extension == ".ply")
    {
        igl::readPLY(input.string(), V, F);
    }
    if (input_extension == ".obj")
    {
        igl::readOBJ(input.string(), V, F);
    }
    if (input_extension == ".stl")
    {
        Eigen::MatrixXd N;
        igl::readSTL(input.string(), V, F, N);
    }

    Eigen::MatrixXd TV;
    Eigen::MatrixXi TF, TT;
    igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414", TV, TT, TF);

    std::filesystem::path const output = argv[2];
    std::string const output_extension = output.extension().string();
    if (output_extension == ".mesh")
    {
        igl::writeMESH(output.string(), TV, TT, TF);
    }
    else
    {
        std::cout << "Expected .mesh extension for output\n";
    }

    return 0;
}