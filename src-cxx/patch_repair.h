// Batched conforming subdivisions. Never add Steiner points or weaken Criteria.
#pragma once

int repair_patch(const std::string& input, const std::string& output,
                 const Criteria& criteria, size_t max_vertices)
{
    if (output.empty() || max_vertices < 3)
        throw std::invalid_argument("Repair requires output and a vertex budget >= 3");
    const auto records = read_json(input);
    const auto started = PatchClock::now();
    const std::string temporary = output + ".tmp";
    std::ofstream out(temporary);
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << "{\"repairs\":[";
    size_t input_vertices = 0, accepted = 0, rejected = 0;
    bool first = true;
    for (const auto& record : records.get_child("repairs")) {
        CDT cdt;
        std::vector<Vertex_handle> boundary;
        OGRLinearRing ring;
        for (const auto& coordinate : record.second.get_child("boundary")) {
            const auto point = patch_point(coordinate.second);
            boundary.push_back(cdt.insert(point));
            ring.addPoint(point.x(), point.y());
        }
        ring.closeRings();
        OGRPolygon domain;
        domain.addRing(&ring);
        if (boundary.size() < 3 || !domain.IsValid() || domain.get_Area() <= 0)
            throw std::invalid_argument("Invalid subdivision boundary");
        for (size_t i = 0; i < boundary.size(); ++i)
            cdt.insert_constraint(boundary[i], boundary[(i + 1) % boundary.size()]);
        input_vertices += cdt.number_of_vertices();
        if (input_vertices > max_vertices)
            throw std::runtime_error("Patch vertex budget exceeded in subdivision batch");
        mark_patch_domain(cdt, &domain);
        auto is_bad = criteria.is_bad_object();
        bool valid = true;
        size_t faces = 0;
        for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
            if (!f->is_in_domain()) continue;
            ++faces;
            Criteria::Quality quality;
            if (is_bad(f, quality) != CGAL::Mesh_2::NOT_BAD) { valid = false; break; }
        }
        valid = valid && faces;
        if (!first) out << ',';
        first = false;
        out << "{\"id\":" << record.second.get<size_t>("id")
            << ",\"accepted\":" << (valid ? "true" : "false") << ",\"triangles\":[";
        if (valid) {
            ++accepted;
            bool first_face = true;
            for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
                if (!f->is_in_domain()) continue;
                if (!first_face) out << ',';
                first_face = false;
                out << '[';
                for (int i = 0; i < 3; ++i) {
                    if (i) out << ',';
                    out << '[' << f->vertex(i)->point().x() << ',' << f->vertex(i)->point().y() << ']';
                }
                out << ']';
            }
        } else ++rejected;
        out << "]}";
    }
    out << "],\"stats\":{\"subdivision_seconds\":" << patch_seconds(started)
        << ",\"raster_evaluations\":" << mesher_raster_stats.evaluations
        << ",\"raster_seconds\":" << mesher_raster_stats.seconds << "}}\n";
    out.close();
    if (!out || std::rename(temporary.c_str(), output.c_str()) != 0)
        throw std::runtime_error("Cannot publish subdivision output");
    std::cout << "Subdivision: accepted=" << accepted << ", rejected=" << rejected << std::endl;
    return 0;
}
