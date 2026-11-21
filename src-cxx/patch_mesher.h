// Local constrained remeshing with the same Criteria as ordinary Mesher.
// Included after read_json() in mesher.cpp.
#pragma once
#include <array>
#include <limits>
#include <map>
#include <set>
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
    size_t quality_repair_vertices = 0;
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
        << ",\"quality_repair_vertices\":" << patch_stats.quality_repair_vertices
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


// The input sector containing this face is bounded by the first constrained
// ray on each side of the vertex. Traversing the fan distinguishes an acute
// sector from its reflex complement, including at an interior stream junction.
bool acute_constraint_sector(const PatchCDT& cdt, Face_handle face, int corner, double bound,
                             std::array<Vertex_handle, 2>* rays = nullptr) {
    const auto vertex = face->vertex(corner);
    auto ray = [&](int direction) -> Vertex_handle {
        auto current = face;
        do {
            const int i = current->index(vertex);
            const int opposite = (i + 3 - direction) % 3;
            if (current->is_constrained(opposite))
                return current->vertex((i + direction) % 3);
            current = current->neighbor(opposite);
        } while (current != face && !cdt.is_infinite(current));
        return Vertex_handle();
    };
    const auto cw = ray(1), ccw = ray(2);
    if (cw == Vertex_handle() || ccw == Vertex_handle() || cw == ccw) return false;
    const auto& p = vertex->point();
    const auto& a = cw->point();
    const auto& b = ccw->point();
    if (CGAL::orientation(p, a, b) != CGAL::LEFT_TURN || CGAL::angle(a, p, b) != CGAL::ACUTE)
        return false;
    const double twice_area = 2 * CGAL::to_double(CGAL::area(p, a, b));
    const bool blocked = twice_area * twice_area < bound * CGAL::to_double(
        CGAL::squared_distance(p, a) * CGAL::squared_distance(p, b));
    if (blocked && rays) *rays = {cw, ccw};
    return blocked;
}

struct PatchFailure {
    Face_handle face;
    double area, angle;
    const char* reason;
    int constrained_corner;
    std::set<size_t> input_interfaces;
    double input_sector_degrees = -1.;
    std::array<Vertex_handle, 2> input_rays{};
};

// Bisect only existing, unconstrained edges. The midpoint must still locate
// in one of its two incident faces after rounding to working coordinates.
// Vertex handles stay valid across insertions. Face handles do not.
size_t bisect_patch_edges(PatchCDT& cdt, const std::vector<PatchFailure>& failures, size_t budget) {
    std::vector<std::pair<Vertex_handle, Vertex_handle>> edges;
    for (const auto& failure : failures) {
        const auto f = failure.face;
        int longest = -1;
        double length = 0;
        for (int i = 0; i < 3; ++i) {
            if (f->is_constrained(i)) continue;
            const double d = CGAL::to_double(cdt.segment(f, i).squared_length());
            if (d > length) { longest = i; length = d; }
        }
        if (longest >= 0) edges.emplace_back(f->vertex((longest+1)%3), f->vertex((longest+2)%3));
    }
    const auto before = cdt.number_of_vertices();
    for (const auto& edge : edges) {
        if (cdt.number_of_vertices() >= budget) break;
        Face_handle f;
        int side;
        if (!cdt.is_edge(edge.first, edge.second, f, side) || f->is_constrained(side)) continue;
        const auto p = CGAL::midpoint(edge.first->point(), edge.second->point());
        PatchCDT::Locate_type location;
        int index;
        auto located = cdt.locate(p, location, index, f);
        if (located != f && located != f->neighbor(side)) continue;
        if (location == PatchCDT::EDGE && located->neighbor(index) !=
                (located == f ? f->neighbor(side) : f)) continue;
        if (location == PatchCDT::EDGE || location == PatchCDT::FACE)
            cdt.insert(p, location, located, index);
    }
    return cdt.number_of_vertices() - before;
}

struct PatchAudit {
    std::vector<PatchFailure> failures;
    size_t faces = 0, bad_area = 0, bad_raster = 0;
    bool blocked = false;
    double worst_angle = 180;
};

// Walk each in-domain sector once, starting on its clockwise constrained
// ray. Splitting either ray or adding interior vertices cannot increase the
// angle between its prescribed directions. Record the actual constraint
// hierarchy IDs so recovery can release the faces that imposed the angle.
PatchAudit audit_input_sectors(const PatchCDT& cdt, double bound,
                              const std::vector<PatchCDT::Constraint_id>& interfaces) {
    std::map<PatchCDT::Constraint_id, size_t> ids;
    for (size_t i = 0; i < interfaces.size(); ++i) ids.emplace(interfaces[i], i);
    PatchAudit audit;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        if (!f->is_in_domain()) continue;
        ++audit.faces;
        PatchFailure failure{f, 0., 180., "angle", -1};
        for (int i = 0; i < 3; ++i) {
            if (!f->is_constrained((i+2)%3)) continue;
            std::array<Vertex_handle, 2> rays;
            if (!acute_constraint_sector(cdt, f, i, bound, &rays)) continue;
            const auto vertex = f->vertex(i);
            const auto u = rays[0]->point() - vertex->point();
            const auto v = rays[1]->point() - vertex->point();
            const double sector = std::atan2(CGAL::to_double(u.x()*v.y()-u.y()*v.x()),
                                             CGAL::to_double(u*v)) * 180. / std::acos(-1.);
            if (failure.constrained_corner < 0 || sector < failure.input_sector_degrees) {
                failure.input_sector_degrees = sector;
                failure.constrained_corner = i;
                failure.input_rays = rays;
            }
            for (auto ray : rays) {
                auto end = cdt.contexts_end(vertex, ray);
                for (auto it = cdt.contexts_begin(vertex, ray); it != end; ++it) {
                    const auto id = ids.find(it->id());
                    if (id != ids.end()) failure.input_interfaces.insert(id->second);
                }
            }
        }
        if (failure.constrained_corner < 0) continue;
        const auto& a = f->vertex(0)->point();
        const auto& b = f->vertex(1)->point();
        const auto& c = f->vertex(2)->point();
        failure.area = CGAL::to_double(CGAL::area(a, b, c));
        for (int i = 0; i < 3; ++i) {
            const auto p = f->vertex(i)->point();
            const auto u = f->vertex((i+1)%3)->point() - p;
            const auto v = f->vertex((i+2)%3)->point() - p;
            const double angle = std::atan2(2*failure.area, CGAL::to_double(u*v)) * 180. / std::acos(-1.);
            failure.angle = std::min(failure.angle, angle);
        }
        audit.failures.push_back(std::move(failure));
    }
    return audit;
}

PatchAudit audit_patch_faces(PatchCDT& cdt, const Criteria& criteria) {
    PatchAudit audit;
    auto is_bad = criteria.is_bad_object();
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        if (!f->is_in_domain()) continue;
        ++audit.faces;
        Criteria::Quality quality;
        if (is_bad(f, quality) == CGAL::Mesh_2::NOT_BAD) continue;
        const char* reason = "raster";
        int corner = -1;
        if (quality.area() > criteria.maxarea()) { reason = "area"; ++audit.bad_area; }
        else if (quality.sine() < criteria.bound()) {
            reason = "angle";
            for (int i = 0; i < 3; ++i) {
                if (acute_constraint_sector(cdt, f, i, criteria.bound())) { corner = i; break; }
            }
            audit.blocked |= corner >= 0;
        } else ++audit.bad_raster;
        const double angle = quality.min_angle_degrees();
        audit.worst_angle = std::min(audit.worst_angle, angle);
        audit.failures.push_back({f, quality.area(), angle, reason, corner});
    }
    return audit;
}

void require_patch_quality(const PatchAudit& audit, const std::string& output) {
    const size_t bad = audit.failures.size();
    if (!bad && audit.faces) return;
    std::ofstream quality_report;
    const std::string report_path = output + ".quality.json";
    if (bad) {
        quality_report.open(report_path + ".tmp");
        quality_report << std::setprecision(std::numeric_limits<double>::max_digits10) << "{\"faces\":[";
    }
    for (size_t j = 0; j < bad; ++j) {
        const auto& failure = audit.failures[j];
        const auto f = failure.face;
        if (j) quality_report << ',';
        quality_report << "{\"reason\":\"" << failure.reason << "\",\"area\":" << failure.area
                       << ",\"min_angle_degrees\":" << failure.angle;
        if (failure.constrained_corner >= 0)
            quality_report << ",\"constrained_corner\":" << failure.constrained_corner;
        if (failure.input_sector_degrees >= 0.) {
            quality_report << ",\"input_sector_degrees\":" << failure.input_sector_degrees
                           << ",\"input_rays\":[";
            for (size_t i = 0; i < failure.input_rays.size(); ++i) {
                const auto& p = failure.input_rays[i]->point();
                if (i) quality_report << ',';
                quality_report << '[' << p.x() << ',' << p.y() << ']';
            }
            quality_report << "],\"input_interface_ids\":[";
            bool first = true;
            for (auto id : failure.input_interfaces) {
                if (!first) quality_report << ',';
                quality_report << id;
                first = false;
            }
            quality_report << ']';
        }
        quality_report << ",\"vertices\":[";
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
        if (j < 10)
            std::cout << "Patch quality failure: reason=" << failure.reason << ", area=" << failure.area
                      << ", min_angle_degrees=" << failure.angle << std::endl;
    }
    if (bad) {
        quality_report << "],\"failed_faces\":" << bad << ",\"area_failures\":" << audit.bad_area
                       << ",\"angle_failures\":" << (bad - audit.bad_area - audit.bad_raster) << ",\"raster_failures\":" << audit.bad_raster << "}\n";
        quality_report.close();
        if (!quality_report || std::rename((report_path + ".tmp").c_str(), report_path.c_str()) != 0)
            throw std::runtime_error("Cannot write patch quality diagnostics");
    }
    // Never silently accept quality failures next to a frozen interface. The
    // caller can enlarge the cavity; small physical angles may still fail.
    throw std::runtime_error("Patch has " + std::to_string(bad) +
        " triangles failing Mesher criteria (or an empty domain)");
}

void record_failed_patch_edge(const std::string& output, const Point& a, const Point& b) {
    std::ofstream report(output + ".edge.json");
    report << std::setprecision(std::numeric_limits<double>::max_digits10)
           << "{\"edge\":[[" << a.x() << ',' << a.y() << "],[" << b.x() << ',' << b.y() << "]]}\n";
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
    bool export_constrained_edges = false;
    if (!interface_file.empty()) {
        const auto records = read_json(interface_file);
        export_constrained_edges = records.get<bool>("export_constrained_edges", false);
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
    const auto input_audit = audit_input_sectors(cdt, criteria.bound(), interfaces);
    if (!input_audit.failures.empty()) require_patch_quality(input_audit, output);

    const auto refine_started = PatchClock::now();
    auto progress = refine_started;
    patch_stats.refinement_started = refine_started;
    patch_stats.refining = true;
    PatchAudit audit;
    size_t previous_count = std::numeric_limits<size_t>::max();
    double previous_angle = 0;
    while (true) {
        PatchMesher mesher(cdt, criteria);
        mesher.init(true);
        while (!mesher.is_refinement_done()) {
            if (cdt.number_of_vertices() >= max_vertices) throw budget_error();
            const auto before = cdt.number_of_vertices();
            if (mesher.is_edges_refinement_done()) {
                mesher.try_one_step_refine_mesh();
            } else {
                const auto edge = mesher.next_encroached_edge();
                // Copy the coordinates: a failed insertion can invalidate
                // face handles. Never call next_refinement_point() merely
                // to preview a split, since it updates CGAL cluster state.
                const auto a = edge.first->vertex((edge.second+1)%3)->point();
                const auto b = edge.first->vertex((edge.second+2)%3)->point();
                if (std::nextafter(a.x(), b.x()) == b.x() && std::nextafter(a.y(), b.y()) == b.y()) {
                    record_failed_patch_edge(output, a, b);
                    throw std::runtime_error("Required constraint split has no representable interior point");
                }
                try {
                    mesher.try_one_step_refine_mesh();
                } catch (const CGAL::Assertion_exception&) {
                    record_failed_patch_edge(output, a, b);
                    throw;
                }
            }
            if (cdt.number_of_vertices() == before) continue;
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

        audit = audit_patch_faces(cdt, criteria);
        // Mesh_2 may drop ordinary BAD faces near protected constraint
        // clusters. Repair only residual shape failures and re-audit the entire
        // patch. Stop when the count and worst angle cease to improve.
        const auto count = audit.failures.size();
        if (!count || audit.blocked || audit.bad_area || audit.bad_raster ||
            count > previous_count || (count == previous_count && audit.worst_angle <= previous_angle)) break;
        previous_count = count;
        previous_angle = audit.worst_angle;
        if (cdt.number_of_vertices() >= max_vertices) throw budget_error();
        const auto inserted = bisect_patch_edges(cdt, audit.failures, max_vertices);
        if (!inserted) break;
        patch_stats.quality_repair_vertices += inserted;
        mark_patch_domain(cdt, domain.get());
    }
    const double refinement_seconds = patch_seconds(refine_started);
    patch_stats.refinement_seconds = refinement_seconds;
    patch_stats.refining = false;
    require_patch_quality(audit, output);

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
    if (export_constrained_edges) {
        out << "],\"constrained_edges\":[";
        first = true;
        for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
            auto f = edge->first;
            const int side = edge->second;
            if (!f->is_constrained(side) ||
                (!f->is_in_domain() && !f->neighbor(side)->is_in_domain())) continue;
            if (!first) out << ',';
            first = false;
            out << '[' << f->vertex((side+1)%3)->info() << ',' << f->vertex((side+2)%3)->info() << ']';
        }
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
        << ",\"quality_repair_vertices\":" << patch_stats.quality_repair_vertices
        << ",\"refinement_seconds\":" << refinement_seconds
        << ",\"raster_evaluations\":" << mesher_raster_stats.evaluations
        << ",\"raster_seconds\":" << mesher_raster_stats.seconds << "}}\n";
    out.close();
    if (!out || std::rename(temporary.c_str(), output.c_str()) != 0)
        throw std::runtime_error("Cannot publish patch output");
    std::cout << "Patch: " << audit.faces << " triangles, all criteria passed" << std::endl;
    return 0;
}
