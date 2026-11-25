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

// Open a source for read-only sampling and transfer the GDAL handle to this wrapper, which also
// validates the transform and first band.
void raster::open(std::string path)
{
    auto opened = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!opened) throw std::runtime_error("Unable to open raster " + path);
    setDs(opened);
}

// Start without a dataset or buffers so destruction is safe even when opening a raster fails
// partway.
raster::raster() : gt(nullptr), ds(nullptr), band(nullptr), mask(nullptr), data(nullptr) {}

// The wrapper owns the sample and mask buffers and the dataset. GDAL owns the band as part of
// that dataset.
raster::~raster()
{
    delete[] gt;
    delete[] mask;
    delete[] data;
    if (ds) GDALClose(ds);
}

// Sample a terrain vertex independently of the polygon burn. GDAL can exclude an edge pixel
// from the area mask even when the vertex lies inside the raster.
double raster::getXY(double x, double y)
{
    auto pixel = xy_to_pxpy(x, y);
    // Vertices on the triangle boundary still sample the source pixel even if
    // GDAL's polygon burn excludes that edge. The mask applies to area samples.
    return getpXpY(pixel.first, pixel.second, false);
}

// Return fractional column and row using the inverse affine transform. Window bounds and point
// sampling both need the original grid's rotation and shear.
std::pair<double,double> raster::pixel_coordinates(double x, double y) const
{
    if (!ds || !gt) throw std::runtime_error("Raster geotransform is unavailable");
    // Subtract the world origin first to avoid cancellation at large UTM values.
    x -= gt[0]; y -= gt[3];
    return {inverse_gt[1] * x + inverse_gt[2] * y,
            inverse_gt[4] * x + inverse_gt[5] * y};
}

// Return the containing pixel index on each axis, admitting roundoff at the outer edge. An
// out-of-range index makes the sample missing.
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

// Borrow the six affine coefficients so a triangle window can preserve the source grid when it
// gets its own origin.
double *raster::getGt() const { return gt; }
// Borrow the dataset for dimensions and projection metadata. Its lifetime is tied to this
// wrapper.
GDALDataset* raster::getDs() const { return ds; }
// Borrow the source band for a window read or a NoData lookup.
GDALRasterBand* raster::getBand() const { return band; }

// Take ownership of a dataset and discard cached samples and georeferencing from the previous
// one. The new band and inverse transform must describe the same grid.
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

// Take ownership of the polygon burn buffer. Its zero entries exclude pixels from area
// statistics.
void raster::setMask(float *values)
{
    delete[] mask;
    mask = values;
}

// Use NaN for missing or excluded samples. A failed GDAL read raises an error so a broken
// raster cannot pass a terrain tolerance check.
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

// Write the samples to the triangle band before accepting ownership. On success pixel reads
// use the same buffer directly. If the write fails the caller still owns it.
void raster::setBand(float *values, int xsize, int ysize)
{
    if (band->RasterIO(GF_Write, 0, 0, xsize, ysize, values, xsize, ysize, GDT_Float32, 0, 0) != CE_None)
        throw std::runtime_error("Cannot update triangle raster band: " + std::string(CPLGetLastErrorMsg()));
    delete[] data;
    data = values;
}
