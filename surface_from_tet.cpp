#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <igl/boundary_facets.h>
#include <igl/readMESH.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writePLY.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

bool verify_tet_mesh_topology(Eigen::MatrixXi const& T)
{
    using face_type = std::array<int, 3u>;

    auto const less = [](face_type const& f1, face_type const& f2) {
        face_type f1_ = f1;
        face_type f2_ = f2;
        std::sort(f1_.begin(), f1_.end());
        std::sort(f2_.begin(), f2_.end());
        return f1_ < f2_;
    };

    std::map<face_type, unsigned int, decltype(less)> face_occurrences(less);
    for (auto t = 0u; t < T.rows(); ++t)
    {
        face_type f1{T(t, 0), T(t, 1), T(t, 3)};
        face_type f2{T(t, 1), T(t, 2), T(t, 3)};
        face_type f3{T(t, 2), T(t, 0), T(t, 3)};
        face_type f4{T(t, 0), T(t, 2), T(t, 1)};

        face_occurrences[f1]++;
        face_occurrences[f2]++;
        face_occurrences[f3]++;
        face_occurrences[f4]++;
    }

    bool const faces_referenced_once_or_twice = std::all_of(
        face_occurrences.begin(),
        face_occurrences.end(),
        [](std::pair<face_type const, unsigned int> const& kv) {
            return kv.second == 1u || kv.second == 2u;
        });

    return faces_referenced_once_or_twice;
}

bool verify_surface(Eigen::MatrixXi const& F)
{
    using edge_type = std::array<int, 2u>;

    auto const less = [](edge_type const& e1, edge_type const& e2) {
        edge_type e1_ = e1;
        edge_type e2_ = e2;
        std::sort(e1_.begin(), e1_.end());
        std::sort(e2_.begin(), e2_.end());
        return e1_ < e2_;
    };

    std::map<edge_type, unsigned int, decltype(less)> edge_occurrences(less);
    for (auto f = 0u; f < F.rows(); ++f)
    {
        edge_type e1{F(f, 0), F(f, 1)};
        edge_type e2{F(f, 1), F(f, 2)};
        edge_type e3{F(f, 2), F(f, 0)};

        edge_occurrences[e1]++;
        edge_occurrences[e2]++;
        edge_occurrences[e3]++;
    }

    bool const edge_referenced_once_or_twice = std::all_of(
        edge_occurrences.begin(),
        edge_occurrences.end(),
        [](std::pair<edge_type const, unsigned int> const& kv) {
            return kv.second == 1u || kv.second == 2u;
        });

    return edge_referenced_once_or_twice;
}

std::vector<int> unreferenced_vertices(Eigen::MatrixXd const& V, Eigen::MatrixXi const& F)
{
    std::vector<unsigned int> vertex_occurrences{};
    vertex_occurrences.resize(V.rows(), 0u);
    for (auto f = 0u; f < F.rows(); ++f)
    {
        for (auto j = 0u; j < 3u; ++j)
        {
            vertex_occurrences[F(f, j)]++;
        }
    }
    std::vector<int> unreferencedV{};
    for (auto i = 0u; i < vertex_occurrences.size(); ++i)
    {
        if (vertex_occurrences[i] > 0u)
            continue;

        unreferencedV.push_back(i);
    }
    return unreferencedV;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> remove_unreferenced_vertices(
    Eigen::MatrixXd const& V,
    Eigen::MatrixXi const& F,
    std::vector<int> const& unreferencedV)
{
    std::unordered_set<int> unrefV(unreferencedV.begin(), unreferencedV.end());
    std::vector<int> prefix_sum{};
    prefix_sum.reserve(V.rows());
    int sum = 0;
    for (auto i = 0u; i < V.rows(); ++i)
    {
        auto it = unrefV.find(i);
        prefix_sum.push_back(sum);
        if (it == unrefV.end())
            ++sum;
    }
    prefix_sum.push_back(sum);

    Eigen::MatrixXd VF;
    Eigen::MatrixXi FF;
    VF.resize(prefix_sum.back(), 3u);
    FF.resize(F.rows(), F.cols());

    auto const is_removed = [&](int i) {
        return prefix_sum[i] == prefix_sum[i + 1];
    };
    for (auto i = 0u; i < V.rows(); ++i)
    {
        if (is_removed(i))
            continue;

        VF.row(prefix_sum[i]) = V.row(i);
    }
    for (auto f = 0u; f < F.rows(); ++f)
    {
        for (auto j = 0u; j < 3u; ++j)
        {
            FF(f, j) = prefix_sum[F(f, j)];
        }
    }

    return {VF, FF};
}

int main(int argc, char** argv)
{
    if (argc == 1 || argc != 4)
    {
        std::cout << "Usage:\n";
        std::cout << "./surface_from_tet.exe input.mesh path/to/output [ply|obj|off]\n";
        return 0;
    }

    std::filesystem::path const input{argv[1]};
    std::string const input_extension = input.extension().string();
    bool const is_valid_input_file    = std::filesystem::exists(input) &&
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
    std::vector<int> const unreferencedV = unreferenced_vertices(V, F);
    std::tie(V, F)                       = remove_unreferenced_vertices(V, F, unreferencedV);
    F                                    = F.rowwise().reverse().eval();

    bool const does_tet_mesh_have_correct_topology  = verify_tet_mesh_topology(T);
    bool const does_face_mesh_have_correct_topology = verify_surface(F);

    std::cout << std::boolalpha << "Tet mesh topology: " << does_tet_mesh_have_correct_topology
              << "\n";
    std::cout << std::boolalpha << "Surface mesh topology: " << does_face_mesh_have_correct_topology
              << "\n";
    std::cout << "Unreferenced vertices: " << unreferencedV.size() << "\n";

    std::filesystem::path const path_to_output = argv[2];
    std::string const output_file_type         = argv[3];
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
