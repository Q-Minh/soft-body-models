#include <igl/readMESH.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>
#include <igl/writeOFF.h>
#include <igl/boundary_facets.h>

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>

int main(int argc, char **argv)
{
    if (argc == 1 || argc != 4)
    {
        std::cout << "Usage:\n";
        std::cout << "./surface_from_tet.exe input.mesh path/to/output [ply|obj|off]\n";
        return 0;
    }

    std::filesystem::path const input{argv[1]};
    std::string const input_extension = input.extension().string();
    bool const is_valid_input_file =
        std::filesystem::exists(input) &&
        std::filesystem::is_regular_file(input) &&
        input_extension == ".mesh";

    if (!is_valid_input_file)
    {
        std::cerr << "Expected input to be .mesh\n";
        return 1;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi T, F;
    igl::readMESH(input.string(), V, T, F);
    F = igl::boundary_facets<Eigen::MatrixXi, Eigen::MatrixXi>(T);
    F.rowwise().reverseInPlace();

    std::filesystem::path const path_to_output = argv[2];
    std::string const output_file_type = argv[3];
    if (output_file_type == "ply")
    {
        std::filesystem::path const output = path_to_output / (input.stem().string() + ".ply");
        igl::writePLY(output.string(), V, F);
    }
    if (output_file_type == "obj")
    {
        std::filesystem::path const output = path_to_output / (input.stem().string() + ".obj");
        igl::writeOBJ(output.string(), V, F);
    }
    if (output_file_type == "off")
    {
        std::filesystem::path const output = path_to_output / (input.stem().string() + ".off");
        igl::writeOFF(output.string(), V, F);
    }

    return 0;
}