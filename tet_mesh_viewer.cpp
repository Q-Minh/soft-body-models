#include <filesystem>
#include <igl/boundary_facets.h>
#include <igl/colormap.h>
#include <igl/file_dialog_open.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readMESH.h>
#include <list>
#include <string>

int main(int argc, char** argv)
{
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    struct TetMesh
    {
        TetMesh(Eigen::MatrixXd const& V, Eigen::MatrixXi const& T, Eigen::MatrixXi const& F)
            : V(V), T(T), F(F)
        {
        }

        Eigen::MatrixXd V;
        Eigen::MatrixXi T, F;
    };

    std::list<TetMesh> hierarchy{};

    menu.callback_draw_viewer_window = [&]() {
        ImGui::SetNextWindowSize(ImVec2(300.0f, 480.0f), ImGuiSetCond_FirstUseEver);
        ImGui::Begin("Tet Mesh Viewer");

        float const w = ImGui::GetContentRegionAvailWidth();
        float const p = ImGui::GetStyle().FramePadding.x;

        if (ImGui::Button("Load tet mesh", ImVec2((w - p) / 2.f, 0)))
        {
            std::string const filename = igl::file_dialog_open();
            std::filesystem::path const mesh{filename};
            if (std::filesystem::exists(mesh) && std::filesystem::is_regular_file(mesh) &&
                mesh.extension().string() == ".mesh")
            {
                Eigen::MatrixXd V;
                Eigen::MatrixXi T, F;
                if (igl::readMESH(mesh.string(), V, T, F))
                {
                    igl::boundary_facets(T, F);
                    F = F.rowwise().reverse().eval();
                    hierarchy.push_back(TetMesh(V, T, F));
                }
            }
        }

        viewer.data().clear();

        Eigen::MatrixXd V, C;
        Eigen::MatrixXi F;
        V.resize(0, 3);
        C.resize(0, 3);
        F.resize(0, 3);
        int i            = 0;
        auto const cmap  = igl::ColorMapType::COLOR_MAP_TYPE_VIRIDIS;
        double const max = static_cast<double>(hierarchy.size());
        double const min = 0.;
        for (auto const& mesh : hierarchy)
        {
            auto const Voffset = V.rows();
            auto const Foffset = F.rows();
            auto const Coffset = C.rows();

            V.conservativeResize(V.rows() + mesh.V.rows(), 3u);
            C.conservativeResize(C.rows() + mesh.V.rows(), 3u);
            F.conservativeResize(F.rows() + mesh.F.rows(), 3u);

            Eigen::VectorXi Z{};
            Z.resize(mesh.V.rows());
            Z.setConstant(static_cast<double>(i));

            Eigen::MatrixXd Cmesh{};
            igl::colormap(cmap, Z, min, max, Cmesh);

            V.block(Voffset, 0, mesh.V.rows(), 3) = mesh.V;
            C.block(Coffset, 0, mesh.V.rows(), 3) = Cmesh;
            F.block(Foffset, 0, mesh.F.rows(), 3) = mesh.F;
            F.block(Foffset, 0, mesh.F.rows(), 3).array() += Voffset;

            ++i;
            ImGui::Text("Mesh %d", i);
            ImGui::Text("Vertices: %d", mesh.V.rows());
            ImGui::Text("Tets: %d", mesh.T.rows());
            ImGui::Text("Boundary faces: %d", mesh.F.rows());
            ImGui::Text("\n");
        }
        viewer.data().set_mesh(V, F);
        viewer.data().set_colors(C);

        ImGui::End();
    };

    viewer.launch();
    return 0;
}