#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/readSTL.h>
#include <igl/writePLY.h>
#include <igl/writeOBJ.h>
#include <igl/writeSTL.h>

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>

int main(int argc, char **argv)
{
    if (argc == 1 || argc != 3)
    {
        std::cout << "Usage:\n";
        std::cout << "./convert_triangle_mesh_format.exe input.[ply|obj|stl] output.[ply|obj|stl]\n";
        return 0;
    }

    std::filesystem::path const input{argv[1]};
    std::string const input_extension = input.extension().string();
    bool const is_valid_input_file =
        std::filesystem::exists(input) &&
        std::filesystem::is_regular_file(input) &&
        (input_extension == ".ply" ||
         input_extension == ".obj" ||
         input_extension == ".stl");

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

    std::filesystem::path const output{argv[2]};
    std::string const output_extension = output.extension().string();
    if (output_extension == ".ply")
    {
        igl::writePLY(output.string(), V, F);
    }
    if (output_extension == ".obj")
    {
        igl::writeOBJ(output.string(), V, F);
    }

    return 0;
}
