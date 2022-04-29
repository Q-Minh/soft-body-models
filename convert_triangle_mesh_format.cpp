#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/readSTL.h>
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>
#include <igl/boundary_loop.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>

/**
 * @brief
 * Taken from https://github.com/libigl/libigl/blob/main/include/igl/copyleft/cgal/mesh_to_polyhedron.cpp .
 * We add the verbose flag to the polyhedron incremental builder to get some useful output.
 *
 * @tparam DerivedV
 * @tparam DerivedF
 * @tparam Polyhedron
 * @param V
 * @param F
 * @param poly
 * @return true
 * @return false
 */
template <
    typename DerivedV,
    typename DerivedF,
    typename Polyhedron>
bool mesh_to_polyhedron(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    Polyhedron &poly)
{
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(poly.hds(), true);
    B.begin_surface(V.rows(), F.rows());
    typedef typename HalfedgeDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
    assert(V.cols() == 3 && "V must be #V by 3");
    for (int v = 0; v < V.rows(); v++)
    {
        B.add_vertex(Point(V(v, 0), V(v, 1), V(v, 2)));
    }
    assert(F.cols() == 3 && "F must be #F by 3");
    for (int f = 0; f < F.rows(); f++)
    {
        B.begin_facet();
        for (int c = 0; c < 3; c++)
        {
            B.add_vertex_to_facet(F(f, c));
        }
        B.end_facet();
    }
    if (B.error())
    {
        B.rollback();
        return false;
    }
    B.end_surface();
    return poly.is_valid();
}

int main(int argc, char **argv)
{
    if (argc == 1 || argc != 2)
    {
        std::cout << "Usage:\n";
        std::cout << "./check_mesh.exe input.[ply|obj]\n";
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
        igl::readSTL(input.string(), V, F);
    }

    Eigen::MatrixXi B{};
    igl::is_vertex_manifold(F, B);
    bool const are_edges_manifold = igl::is_edge_manifold(F);
    std::vector<int> vertex_loop{};
    igl::boundary_loop(F, vertex_loop);

    std::cout << "Vertices: " << V.rows() << "\n";
    std::cout << "Faces: " << F.rows() << "\n";
    std::cout << "Manifold vertices: " << B.array().sum() << "\n";
    std::cout << "Manifold edges: " << std::boolalpha << are_edges_manifold << "\n";
    std::cout << "Watertight: " << std::boolalpha << vertex_loop.empty() << "\n";

    using Kernel = CGAL::Simple_cartesian<double>;
    using SurfaceMesh = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;

    SurfaceMesh M{};
    bool const mesh_to_polyhedron_success = mesh_to_polyhedron(V, F, M);

    std::cout << "Mesh to CGAL polyhedron test: " << mesh_to_polyhedron_success << "\n";
    std::cout << "CGAL Polyhedron vertices: " << M.size_of_vertices() << "\n";
    std::cout << "CGAL Polyhedron facets: " << M.size_of_facets() << "\n";

    return 0;
}