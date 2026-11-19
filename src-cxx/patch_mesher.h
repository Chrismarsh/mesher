// Local constrained remeshing with the same Criteria as ordinary Mesher.
// Included after read_json() in mesher.cpp.
#pragma once
#include <limits>
#include <memory>
#include <cstdio>
#include <chrono>
#include <CGAL/Constrained_triangulation_plus_2.h>

using PatchCDT = CGAL::Constrained_triangulation_plus_2<CDT>;
using PatchMesher = CGAL::Delaunay_mesher_2<PatchCDT, Criteria>;
using PatchClock = std::chrono::steady_clock;

double patch_seconds(PatchClock::time_point start) {
    return std::chrono::duration<double>(PatchClock::now() - start).count();
}

struct PatchStatistics {
    size_t refinement_steps = 0;
    double refinement_seconds = 0;
    bool refining = false;
    PatchClock::time_point refinement_started;
};
static PatchStatistics patch_stats;

void write_patch_statistics(const std::string& output) {
    if (output.empty()) return;
    const double refinement_seconds = patch_stats.refining
        ? patch_seconds(patch_stats.refinement_started) : patch_stats.refinement_seconds;
    std::ofstream out(output + ".stats.json");
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "{\"refinement_steps\":" << patch_stats.refinement_steps
        << ",\"refinement_seconds\":" << refinement_seconds
        << ",\"raster_evaluations\":" << mesher_raster_stats.evaluations
        << ",\"raster_seconds\":" << mesher_raster_stats.seconds << "}\n";
}

Point patch_point(const pt::ptree& coordinate) {
    if (coordinate.size() < 2) throw std::invalid_argument("Patch coordinate needs x and y");
    auto it = coordinate.begin();
    const double x = it++->second.get_value<double>(), y = it->second.get_value<double>();
    if (!std::isfinite(x) || !std::isfinite(y)) throw std::invalid_argument("Non-finite patch coordinate");
    return Point(x, y);
}

template<class Triangulation>
void mark_patch_domain(Triangulation& cdt, OGRGeometry* domain) {
    // A single rounded centroid near a raster-stair boundary can falsely mark
    // an exterior sliver as inside. Mark whole regions separated by constrained
    // edges instead. The region touching infinity is ALWAYS outside; classify
    // other regions using their largest face to avoid near-collinear samples.
    // This also distinguishes holes from regions enclosed by physical streams.
    for (auto f = cdt.all_faces_begin(); f != cdt.all_faces_end(); ++f) f->id = -1;
    int region = 0;
    for (auto seed = cdt.all_faces_begin(); seed != cdt.all_faces_end(); ++seed) {
        if (seed->id >= 0) continue;
        std::vector<Face_handle> faces{seed};
        seed->id = region++;
        bool exterior = false;
        double largest_area = 0;
        Face_handle representative;
        for (size_t i = 0; i < faces.size(); ++i) {
            auto face = faces[i];
            if (cdt.is_infinite(face)) exterior = true;
            else {
                const double area = CGAL::to_double(CGAL::area(face->vertex(0)->point(),
                    face->vertex(1)->point(), face->vertex(2)->point()));
                if (area > largest_area) { largest_area = area; representative = face; }
            }
            for (int side = 0; side < 3; ++side) {
                auto neighbour = face->neighbor(side);
                if (!face->is_constrained(side) && neighbour->id < 0) {
                    neighbour->id = seed->id;
                    faces.push_back(neighbour);
                }
            }
        }
        bool inside = false;
        if (!exterior && largest_area > 0) {
            const auto p = CGAL::centroid(representative->vertex(0)->point(),
                representative->vertex(1)->point(), representative->vertex(2)->point());
            OGRPoint center(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
            inside = domain->Contains(&center);
        }
        for (auto face : faces) face->set_in_domain(inside);
    }
}


int mesh_patch(const std::string& domain_file, const std::string& constraints_file,
               const std::string& output, const Criteria& criteria, size_t max_vertices,
               const std::string& interface_file = "")
{
    if (output.empty() || constraints_file.empty() || max_vertices < 3)
        throw std::invalid_argument("Patch output, constraints and a vertex budget >= 3 are required");
    std::ifstream input(domain_file);
    std::string json((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    std::unique_ptr<OGRGeometry> domain(OGRGeometryFactory::createFromGeoJson(json.c_str()));
    if (!domain || domain->IsEmpty() || !domain->IsValid())
        throw std::invalid_argument("Patch domain must be nonempty valid polygon geometry");
    auto kind = wkbFlatten(domain->getGeometryType());
    if (kind != wkbPolygon && kind != wkbMultiPolygon)
        throw std::invalid_argument("Patch domain must be Polygon or MultiPolygon");

    const auto started = PatchClock::now();
    PatchCDT cdt;
    // Python supplies all boundary and physical constraints, noded together.
    // The old MPI seam is deliberately absent from this graph.
    const auto constraints = read_json(constraints_file);
    for (const auto& feature : constraints.get_child("features")) {
        if (feature.second.get<std::string>("geometry.type") != "LineString")
            throw std::invalid_argument("Patch constraints must be noded LineStrings");
        std::vector<Vertex_handle> vertices;
        for (const auto& coordinate : feature.second.get_child("geometry.coordinates")) {
            if (coordinate.second.size() < 2)
                throw std::invalid_argument("Patch coordinates need x and y");
            auto it = coordinate.second.begin();
            double x = it++->second.get_value<double>();
            double y = it->second.get_value<double>();
            if (!std::isfinite(x) || !std::isfinite(y))
                throw std::invalid_argument("Non-finite patch coordinate");
            vertices.push_back(cdt.insert(Point(x, y)));
        }
        for (size_t i = 1; i < vertices.size(); ++i)
            if (vertices[i - 1] != vertices[i]) cdt.insert_constraint(vertices[i - 1], vertices[i]);
    }
    // Track each caller-owned interface through CGAL's constraint hierarchy.
    // Its ordered children are exported by ID, never rediscovered by proximity.
    std::vector<PatchCDT::Constraint_id> interfaces;
    if (!interface_file.empty()) {
        const auto records = read_json(interface_file);
        for (const auto& record : records.get_child("edges")) {
            std::vector<Point> points;
            for (const auto& point : record.second) points.push_back(patch_point(point.second));
            if (points.size() < 2 || points.front() == points.back())
                throw std::invalid_argument("Invalid retained interface");
            interfaces.push_back(cdt.insert_constraint(points.begin(), points.end()));
        }
    }
    if (cdt.dimension() != 2)
        throw std::runtime_error("Patch constraint graph is not two-dimensional");
    const size_t input_vertices = cdt.number_of_vertices();
    size_t steps = 0;
    auto budget_error = [&]() {
        return std::runtime_error("Patch vertex budget exceeded: vertices=" + std::to_string(cdt.number_of_vertices()) +
            ", limit=" + std::to_string(max_vertices) + ", input_vertices=" + std::to_string(input_vertices) +
            ", refinement_steps=" + std::to_string(steps) +
            ". Check boundary detail and area/tolerance settings before increasing mpi_remesh_max_vertices.");
    };
    if (input_vertices > max_vertices) throw budget_error();
    std::cout << "Patch input: vertices=" << input_vertices << ", budget=" << max_vertices << std::endl;

    mark_patch_domain(cdt, domain.get());
    PatchMesher mesher(cdt, criteria);
    const auto refine_started = PatchClock::now();
    auto progress = refine_started;
    patch_stats.refinement_started = refine_started;
    patch_stats.refining = true;
    mesher.init(true);
    while (!mesher.is_refinement_done()) {
        if (cdt.number_of_vertices() >= max_vertices)
            throw budget_error();
        mesher.step_by_step_refine_mesh();
        ++steps;
        ++patch_stats.refinement_steps;
        if ((steps % 128) == 0 && patch_seconds(progress) >= 10) {
            std::cout << "Patch progress: vertices=" << cdt.number_of_vertices()
                      << ", refinement_steps=" << steps << ", elapsed_seconds=" << patch_seconds(started)
                      << ", raster_evaluations=" << mesher_raster_stats.evaluations << std::endl;
            progress = PatchClock::now();
        }
    }
    if (cdt.number_of_vertices() > max_vertices) throw budget_error();

    const double refinement_seconds = patch_seconds(refine_started);
    patch_stats.refinement_seconds = refinement_seconds;
    patch_stats.refining = false;
    size_t bad = 0, faces = 0;
    size_t bad_area = 0, bad_angle = 0, bad_raster = 0;
    std::ofstream quality_report;
    const std::string report_path = output + ".quality.json";
    auto is_bad = criteria.is_bad_object();
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        if (!f->is_in_domain()) continue;
        ++faces;
        Criteria::Quality quality;
        if (is_bad(f, quality) == CGAL::Mesh_2::NOT_BAD) continue;
        if (!bad) {
            quality_report.open(report_path + ".tmp");
            quality_report << std::setprecision(std::numeric_limits<double>::max_digits10)
                           << "{\"faces\":[";
        } else quality_report << ',';
        ++bad;
        const char* reason = "raster";
        if (quality.area() > criteria.maxarea()) { reason = "area"; ++bad_area; }
        else if (quality.sine() < criteria.bound()) { reason = "angle"; ++bad_angle; }
        else ++bad_raster;
        const double angle = std::asin(std::sqrt(std::max(0.0, std::min(1.0, quality.sine()))))
                           * 180 / std::acos(-1.0);
        quality_report << "{\"reason\":\"" << reason << "\",\"area\":" << quality.area()
                       << ",\"min_angle_degrees\":" << angle << ",\"vertices\":[";
        for (int i = 0; i < 3; ++i) {
            if (i) quality_report << ',';
            quality_report << '[' << f->vertex(i)->point().x() << ',' << f->vertex(i)->point().y() << ']';
        }
        quality_report << "],\"constrained_edges\":[";
        for (int i = 0; i < 3; ++i) {
            if (i) quality_report << ',';
            quality_report << (f->is_constrained(i) ? "true" : "false");
        }
        quality_report << "]}";
        if (bad <= 10)
            std::cout << "Patch quality failure: reason=" << reason << ", area=" << quality.area()
                      << ", min_angle_degrees=" << angle << std::endl;
    }
    if (bad) {
        quality_report << "],\"failed_faces\":" << bad << ",\"area_failures\":" << bad_area
                       << ",\"angle_failures\":" << bad_angle << ",\"raster_failures\":" << bad_raster << "}\n";
        quality_report.close();
        if (!quality_report || std::rename((report_path + ".tmp").c_str(), report_path.c_str()) != 0)
            throw std::runtime_error("Cannot write patch quality diagnostics");
    }
    // Never silently accept quality failures next to a frozen interface. The
    // caller can enlarge the cavity; small physical angles may still fail.
    if (!faces || bad) throw std::runtime_error("Patch has " + std::to_string(bad) +
        " triangles failing Mesher criteria (or an empty domain)");

    const std::string temporary = output + ".tmp";
    std::ofstream out(temporary);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "{\"vertices\":[";
    size_t id = 0;
    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
        v->info() = id;
        if (id++) out << ',';
        out << '[' << v->point().x() << ',' << v->point().y() << ",0]";
    }
    out << "],\"triangles\":[";
    bool first = true;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        if (!f->is_in_domain()) continue;
        if (!first) out << ',';
        first = false;
        out << '[' << f->vertex(0)->info() << ',' << f->vertex(1)->info() << ',' << f->vertex(2)->info() << ']';
    }
    out << "],\"interfaces\":[";
    for (size_t i = 0; i < interfaces.size(); ++i) {
        if (i) out << ',';
        out << "{\"id\":" << i << ",\"vertices\":[";
        bool first_vertex = true;
        for (auto v = cdt.vertices_in_constraint_begin(interfaces[i]);
             v != cdt.vertices_in_constraint_end(interfaces[i]); ++v) {
            if (!first_vertex) out << ',';
            first_vertex = false;
            out << (*v)->info();
        }
        out << "]}";
    }
    out << "],\"criteria_failures\":0,\"stats\":{\"refinement_steps\":" << steps
        << ",\"refinement_seconds\":" << refinement_seconds
        << ",\"raster_evaluations\":" << mesher_raster_stats.evaluations
        << ",\"raster_seconds\":" << mesher_raster_stats.seconds << "}}\n";
    out.close();
    if (!out || std::rename(temporary.c_str(), output.c_str()) != 0)
        throw std::runtime_error("Cannot publish patch output");
    std::cout << "Patch: " << faces << " triangles, all criteria passed" << std::endl;
    return 0;
}
