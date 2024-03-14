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
#pragma once
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

#include <iostream>
#include <string>
#include <utility>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

class raster
{
public:
    raster();

    virtual ~raster();

    void open(std::string path);
    void create_memory_raster(double xsize, double ysize, GDALDataType type);

    double getXY(double X, double Y);
    double getpXpY(int px, int py);
    GDALRasterBand *getBand() const;

    GDALDataset *getDs() const;
    void setDs(GDALDataset* ds);
    std::pair<int,int> xy_to_pxpy(double pX, double pY);
    void setMask(float *mask);
    void setBand(float *data, int xsize, int ysize);
    double *getGt() const;

private:
    double* gt;
    GDALDataset* ds;
    GDALRasterBand* band;
    float* mask;
    float* data;



};
