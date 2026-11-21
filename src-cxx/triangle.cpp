// Mesher
// Copyright (C) 2017 Christopher Marsh

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "triangle.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

triangle::triangle() : is_nan(false) {}
triangle::~triangle() = default;

void triangle::make_rasterized(vertex v0_in, vertex v1_in, vertex v2_in, const raster& r)
{
    is_nan = false;
    rasterized_triangle.reset();
    const auto gt = r.getGt();
    const auto source = r.getDs();
    if (!source || !r.getBand())
        throw std::runtime_error("Cannot sample a triangle from an unopened raster");

    OGRLinearRing ring;
    double xmin = std::numeric_limits<double>::infinity(), ymin = xmin;
    double xmax = -xmin, ymax = -xmin;
    for (const auto v : {v0_in, v1_in, v2_in}) {
        auto pixel = r.pixel_coordinates(v[0], v[1]);
        if (!std::isfinite(pixel.first) || !std::isfinite(pixel.second))
            throw std::runtime_error("Non-finite triangle pixel coordinates in " + std::string(source->GetDescription()));
        xmin = std::min(xmin, pixel.first); xmax = std::max(xmax, pixel.first);
        ymin = std::min(ymin, pixel.second); ymax = std::max(ymax, pixel.second);
        ring.addPoint(v[0], v[1]);
    }
    ring.closeRings();
    OGRPolygon poly;
    poly.addRing(&ring);

    // Intersect in floating-point pixel space BEFORE converting to integers or
    // allocating. A disjoint triangle can otherwise create negative dimensions
    // (converted to huge unsigned sizes by GDAL), or overflow the integer cast.
    const int width = source->GetRasterXSize(), height = source->GetRasterYSize();
    if (xmax <= 0 || ymax <= 0 || xmin >= width || ymin >= height) {
        is_nan = true; // No samples available, as for an all-NoData triangle.
        return;
    }
    auto clipped_start = [](double value, int limit) {
        return static_cast<int>(std::max(0.0, std::min(double(limit), std::floor(value))));
    };
    auto clipped_end = [](double value, int limit) {
        return static_cast<int>(std::max(0.0, std::min(double(limit), std::floor(value) + 1.0)));
    };
    const int x1 = clipped_start(xmin, width), y1 = clipped_start(ymin, height);
    const int x2 = clipped_end(xmax, width), y2 = clipped_end(ymax, height);
    const int xsize = x2 - x1, ysize = y2 - y1;
    if (xsize <= 0 || ysize <= 0) {
        is_nan = true;
        return;
    }
    if (size_t(xsize) > std::numeric_limits<size_t>::max() / size_t(ysize) / sizeof(float))
        throw std::runtime_error("Triangle raster window is too large to allocate");
    const size_t count = size_t(xsize) * size_t(ysize);
    std::ostringstream context;
    context << "Triangle raster window (" << x1 << ',' << y1 << ") "
            << xsize << 'x' << ysize << " in " << source->GetDescription();
    auto check = [&](CPLErr error, const char* action) {
        if (error != CE_None)
            throw std::runtime_error(context.str() + ": " + action + ": " + CPLGetLastErrorMsg());
    };

    double new_gt[6];
    std::copy(gt, gt + 6, new_gt);
    GDALApplyGeoTransform(gt, x1, y1, &new_gt[0], &new_gt[3]);
    auto mem_drv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!mem_drv) throw std::runtime_error("GDAL MEM driver is unavailable");
    auto rvds = mem_drv->Create("", xsize, ysize, 1, GDT_Float32, nullptr);
    if (!rvds) throw std::runtime_error(context.str() + ": MEM allocation failed: " + CPLGetLastErrorMsg());
    // Take ownership immediately so failed RasterIO/rasterization also closes it.
    std::unique_ptr<GDALDataset, decltype(&GDALClose)> owner(rvds, GDALClose);
    rasterized_triangle = boost::make_shared<raster>();
    check(rvds->SetGeoTransform(new_gt), "setting geotransform");
    rasterized_triangle->setDs(owner.release());
    check(rvds->SetProjection(source->GetProjectionRef()), "setting projection");

    int band_number = 1;
    double burn = 1;
    OGRGeometryH geometry = reinterpret_cast<OGRGeometryH>(&poly);
    const char* options[] = {"ALL_TOUCHED=TRUE", nullptr};
    check(GDALRasterizeGeometries(rvds, 1, &band_number, 1, &geometry,
          nullptr, nullptr, &burn, options, nullptr, nullptr), "rasterizing triangle");

    std::unique_ptr<float[]> mask(new float[count]);
    check(rvds->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, xsize, ysize,
          mask.get(), xsize, ysize, GDT_Float32, 0, 0), "reading triangle mask");
    rasterized_triangle->setMask(mask.release());
    std::unique_ptr<float[]> data(new float[count]);
    check(r.getBand()->RasterIO(GF_Read, x1, y1, xsize, ysize,
          data.get(), xsize, ysize, GDT_Float32, 0, 0), "reading source raster");
    int has_nodata = false;
    double nodata = r.getBand()->GetNoDataValue(&has_nodata);
    if (has_nodata) check(rvds->GetRasterBand(1)->SetNoDataValue(nodata), "copying NoData");
    rasterized_triangle->setBand(data.get(), xsize, ysize);
    data.release();

    v0[0] = v0_in[0]; v0[1] = v0_in[1];
    v1[0] = v1_in[0]; v1[1] = v1_in[1];
    v2[0] = v2_in[0]; v2[1] = v2_in[1];
    v0[2] = rasterized_triangle->getXY(v0[0], v0[1]);
    v1[2] = rasterized_triangle->getXY(v1[0], v1[1]);
    v2[2] = rasterized_triangle->getXY(v2[0], v2[1]);
    is_nan = std::isnan(v0[2]) && std::isnan(v1[2]) && std::isnan(v2[2]);
}
