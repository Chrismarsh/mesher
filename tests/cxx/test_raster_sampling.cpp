#include "triangle.h"
#include <array>
#include <climits>
#include <iostream>
#include <stdexcept>
#include <string>

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

void sampling_cases(std::array<double, 6> gt)
{
    auto ds = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 10, 8, 1, GDT_Float32, nullptr);
    require(ds != nullptr, "Create test raster");
    require(ds->SetGeoTransform(gt.data()) == CE_None, "Set transform");
    require(ds->GetRasterBand(1)->Fill(0) == CE_None, "Fill zero-valued raster");
    raster r;
    r.setDs(ds);
    auto make = [&](double x0, double y0, double x1, double y1, double x2, double y2) {
        vertex a, b, c;
        GDALApplyGeoTransform(gt.data(), x0, y0, &a[0], &a[1]);
        GDALApplyGeoTransform(gt.data(), x1, y1, &b[0], &b[1]);
        GDALApplyGeoTransform(gt.data(), x2, y2, &c[0], &c[1]);
        triangle t;
        CPLErrorReset();
        t.make_rasterized(a, b, c, r);
        require(CPLGetLastErrorType() < CE_Failure, "No invalid GDAL reads or allocations");
        return t;
    };
    // Partial overlap at every edge and corner, with a genuine zero elevation.
    for (auto coords : {std::array<double,6>{-2, 2, 2, 2, 2, 6}, {8, 2, 12, 2, 8, 6},
                        {2, -2, 2, 2, 6, 2}, {2, 6, 2, 10, 6, 6},
                        {-2, -2, 3, 2, 2, 3}, {8, 6, 12, 10, 7, 6}}) {
        auto t = make(coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
        require(!t.is_nan && t.rasterized_triangle, "Retain samples from partially overlapping triangle");
        require(std::isnan(t.v0[2]) || t.v0[2] == 0, "Finite vertex must equal source zero");
        auto sampled = t.rasterized_triangle;
        require(sampled->getDs()->GetRasterXSize() <= 10 && sampled->getDs()->GetRasterYSize() <= 8,
                "Window remains within source size");
        require(std::isnan(sampled->getpXpY(-1, 0)), "Reject negative pixel index");
        require(std::isnan(sampled->getpXpY(INT_MAX, INT_MAX)), "Reject huge pixel index");
    }
    // Entirely outside: no negative dimensions, allocation or integer overflow.
    for (auto coords : {std::array<double,6>{-15, 2, -6, 2, -6, 7}, {11, 2, 15, 2, 11, 6},
                        {2, -8, 2, -2, 6, -2}, {2, 9, 2, 15, 6, 9},
                        {-4, -4, -2, -4, -4, -2}, {12, 10, 15, 10, 12, 13},
                        {1e100, 2, 1e100, 6, 2e100, 2}}) {
        auto t = make(coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
        require(t.is_nan && !t.rasterized_triangle, "Disjoint triangle has no raster samples");
    }
    auto inside = make(1, 1, 8, 1, 1, 6);
    require(!inside.is_nan && inside.v0[2] == 0 && inside.v1[2] == 0 && inside.v2[2] == 0,
            "Zero is valid without a NoData tag");
    require(ds->GetRasterBand(1)->SetNoDataValue(-9999) == CE_None, "Set NoData");
    require(ds->GetRasterBand(1)->Fill(-9999) == CE_None, "Fill NoData");
    auto missing = make(1, 1, 8, 1, 1, 6);
    require(missing.is_nan, "Preserve source NoData through the temporary raster");
}

void failed_read_case()
{
    const char* vrt = "<VRTDataset rasterXSize='10' rasterYSize='10'>"
        "<GeoTransform>0,1,0,10,0,-1</GeoTransform><VRTRasterBand dataType='Float32' band='1'>"
        "<SimpleSource><SourceFilename>/nonexistent/mesher-test-raster.tif</SourceFilename>"
        "<SourceBand>1</SourceBand><SrcRect xOff='0' yOff='0' xSize='10' ySize='10'/>"
        "<DstRect xOff='0' yOff='0' xSize='10' ySize='10'/></SimpleSource>"
        "</VRTRasterBand></VRTDataset>";
    raster r;
    r.setDs(static_cast<GDALDataset*>(GDALOpen(vrt, GA_ReadOnly)));
    vertex a = {1,9,0}, b = {8,9,0}, c = {1,2,0};
    bool rejected = false;
    CPLPushErrorHandler(CPLQuietErrorHandler);
    try {
        triangle t;
        t.make_rasterized(a, b, c, r);
    } catch (const std::runtime_error& error) {
        rejected = std::string(error.what()).find("reading source raster") != std::string::npos;
    }
    CPLPopErrorHandler();
    require(rejected, "Failed RasterIO must raise, never evaluate uninitialized data");
}

int main()
{
    GDALAllRegister();
    try {
        sampling_cases({0, 1, 0, 8, 0, -1});
        sampling_cases({500000, 2, .5, 7600000, .25, 3});
        failed_read_case();
        std::cout << "Raster sampling regressions passed" << std::endl;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
