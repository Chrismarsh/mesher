#include "mesh.h"
#include <iostream>
#include <stdexcept>
#include <vector>

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void mandatory_raster_cases()
{
    using CGAL::Mesh_2::BAD;
    using CGAL::Mesh_2::NOT_BAD;
    using CGAL::Mesh_2::IMPERATIVELY_BAD;
    std::vector<std::tuple<boost::shared_ptr<raster>, double, double>> rasters, categories;
    const std::string metric = "rmse";
    Criteria criteria(.125, 100, 1, rasters, categories, metric);
    auto bad = criteria.is_bad_object();
    Criteria::Quality q;
    q._area = 20; q._size = .2; q._sine = .7;
    require(bad(q) == NOT_BAD, "Well-shaped geometry passes");
    q._tolerance.emplace_back(50, 30);
    require(bad(q) == IMPERATIVELY_BAD, "Numeric raster error must not be dropped at protected clusters");
    q._tolerance = {{0, -1}};
    require(bad(q) == NOT_BAD, "A disabled raster must not fail a cached check");
    q._category_tol.emplace_back(.4, .6);
    require(bad(q) == IMPERATIVELY_BAD, "Categorical error must request mandatory refinement");
    q._category_tol.clear(); q._weighted_raster_bad = true;
    require(bad(q) == IMPERATIVELY_BAD, "Preserve weighted failure when CGAL reads cached Quality");
    q._area = 1;
    require(bad(q) == NOT_BAD, "Minimum area still stops raster refinement");
    q._sine = .01;
    require(bad(q) == BAD, "Shape-only refinement remains nonmandatory even below minimum area");
    q._area = 101;
    require(bad(q) == IMPERATIVELY_BAD, "Maximum area remains mandatory");

    // Exercise both predicate overloads on an actual raster and CDT face.
    auto source = boost::make_shared<raster>();
    auto ds = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 8, 8, 1, GDT_Float32, nullptr);
    double gt[] = {0, 1, 0, 8, 0, -1};
    require(ds && ds->SetGeoTransform(gt) == CE_None, "Create raster fixture");
    require(ds->GetRasterBand(1)->Fill(0) == CE_None, "Fill raster fixture");
    float value = 100;
    require(ds->GetRasterBand(1)->RasterIO(GF_Write, 3, 5, 1, 1, &value, 1, 1, GDT_Float32, 0, 0) == CE_None,
            "Write interior raster anomaly");
    source->setDs(ds);
    rasters.emplace_back(source, 1, 1);
    CDT mesh;
    mesh.insert(Point(1, 1)); mesh.insert(Point(6, 1)); mesh.insert(Point(3, 6));
    auto face = mesh.finite_faces_begin();
    for (bool weighted : {false, true}) {
        Criteria check(.125, 100, 1, rasters, categories, metric, false, weighted, 1);
        auto evaluate = check.is_bad_object();
        Criteria::Quality quality;
        require(evaluate(face, quality) == IMPERATIVELY_BAD, "Raster anomaly requires refinement");
        require(evaluate(quality) == IMPERATIVELY_BAD, "Cached and face checks must agree");
        if (weighted) {
            Criteria permissive(.125, 100, 1, rasters, categories, metric, false, true, 0);
            Criteria::Quality passed;
            auto accept = permissive.is_bad_object();
            require(accept(face, passed) == NOT_BAD && accept(passed) == NOT_BAD,
                    "A satisfied weighted score must remain accepted in both checks");
        }
    }
}

int main()
{
    GDALAllRegister();
    try {
        mandatory_raster_cases();
        std::vector<Criteria::Quality> values;
        for (double size : {.01, .5, 1., 2., 100.}) {
            for (double sine : {.001, .125, .5, .75}) {
                Criteria::Quality q;
                q._area = 6250000 * size;
                q._size = size;
                q._sine = sine;
                values.push_back(q);
                q._tolerance.emplace_back(50, 30);
                q._category_tol.emplace_back(.4, .6);
                values.push_back(q);
            }
        }
        auto equivalent = [](const Criteria::Quality& a, const Criteria::Quality& b) {
            return !(a < b) && !(b < a);
        };
        for (auto a : values) {
            require(!(a < a), "Queue order must be irreflexive");
            for (auto b : values) {
                require(!(a < b && b < a), "Queue order must be asymmetric");
                if (a._size > 1 && b._size <= 1)
                    require(a < b, "Oversized faces must be refined first");
                if (a._size > 1 && b._size > 1 && a._size > b._size)
                    require(a < b, "Largest area violation must be refined first");
                if (a._size <= 1 && b._size <= 1 && a._sine < b._sine)
                    require(a < b, "Smallest angle must be refined first below the area limit");
                for (auto c : values) {
                    require(!(a < b && b < c) || a < c, "Queue order must be transitive");
                    require(!(equivalent(a, b) && equivalent(b, c)) || equivalent(a, c),
                            "Queue equivalence must be transitive");
                }
            }
        }
        std::cout << "Mesh criteria ordering regressions passed" << std::endl;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
