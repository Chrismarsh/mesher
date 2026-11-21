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

#include "raster.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

void raster::open(std::string path)
{
    auto opened = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!opened) throw std::runtime_error("Unable to open raster " + path);
    setDs(opened);
}

raster::raster() : gt(nullptr), ds(nullptr), band(nullptr), mask(nullptr), data(nullptr) {}

raster::~raster()
{
    delete[] gt;
    delete[] mask;
    delete[] data;
    if (ds) GDALClose(ds);
}

double raster::getXY(double x, double y)
{
    auto pixel = xy_to_pxpy(x, y);
    // Vertices on the triangle boundary still sample the source pixel even if
    // GDAL's polygon burn excludes that edge. The mask applies to area samples.
    return getpXpY(pixel.first, pixel.second, false);
}

std::pair<double,double> raster::pixel_coordinates(double x, double y) const
{
    if (!ds || !gt) throw std::runtime_error("Raster geotransform is unavailable");
    // Subtract the world origin first to avoid cancellation at large UTM values.
    x -= gt[0]; y -= gt[3];
    return {inverse_gt[1] * x + inverse_gt[2] * y,
            inverse_gt[4] * x + inverse_gt[5] * y};
}

std::pair<int,int> raster::xy_to_pxpy(double x, double y)
{
    auto pixel = pixel_coordinates(x, y);
    const double precision = 64 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(x), std::abs(y), std::abs(gt[0]), std::abs(gt[3])});
    auto index = [](double value, int size, double precision) {
        // Permit coordinate roundoff on the outer edge, not a whole pixel of
        // extrapolation. Out-of-range values use a sentinel checked below.
        const double epsilon = std::max(precision, 64 * std::numeric_limits<double>::epsilon() * std::max(1, size));
        if (!std::isfinite(value) || value < -epsilon || value > size + epsilon)
            return -1;
        value = std::max(0.0, std::min(double(size - 1), std::floor(value)));
        return static_cast<int>(value);
    };
    return {index(pixel.first, ds->GetRasterXSize(), precision * (std::abs(inverse_gt[1]) + std::abs(inverse_gt[2]))),
            index(pixel.second, ds->GetRasterYSize(), precision * (std::abs(inverse_gt[4]) + std::abs(inverse_gt[5])))};
}

double *raster::getGt() const { return gt; }
GDALDataset* raster::getDs() const { return ds; }
GDALRasterBand* raster::getBand() const { return band; }

void raster::setDs(GDALDataset* source)
{
    if (ds && ds != source) GDALClose(ds);
    ds = source; // This object owns the dataset, including on validation failure.
    delete[] data; data = nullptr;
    delete[] mask; mask = nullptr;
    delete[] gt; gt = new double[6];
    if (!ds || ds->GetGeoTransform(gt) != CE_None || !GDALInvGeoTransform(gt, inverse_gt))
        throw std::runtime_error("Raster must have an invertible geotransform");
    band = ds->GetRasterBand(1);
    if (!band) throw std::runtime_error("Raster has no data band");
}

void raster::setMask(float *values)
{
    delete[] mask;
    mask = values;
}

double raster::getpXpY(int px, int py, bool apply_mask)
{
    const double missing = std::numeric_limits<double>::quiet_NaN();
    if (!ds || !band || px < 0 || py < 0 || px >= ds->GetRasterXSize() || py >= ds->GetRasterYSize())
        return missing;
    const size_t index = size_t(py) * size_t(ds->GetRasterXSize()) + size_t(px);
    if (apply_mask && mask && mask[index] == 0) return missing;
    float value;
    if (data) value = data[index];
    else if (band->RasterIO(GF_Read, px, py, 1, 1, &value, 1, 1, GDT_Float32, 0, 0) != CE_None)
        throw std::runtime_error("Cannot read raster sample: " + std::string(CPLGetLastErrorMsg()));
    int has_nodata = false;
    const double nodata = band->GetNoDataValue(&has_nodata);
    if (!std::isfinite(value) || (has_nodata && value == static_cast<float>(nodata))) return missing;
    return value;
}

void raster::setBand(float *values, int xsize, int ysize)
{
    if (band->RasterIO(GF_Write, 0, 0, xsize, ysize, values, xsize, ysize, GDT_Float32, 0, 0) != CE_None)
        throw std::runtime_error("Cannot update triangle raster band: " + std::string(CPLGetLastErrorMsg()));
    delete[] data;
    data = values;
}
