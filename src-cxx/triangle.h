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

#include <utility>
#include <boost/shared_ptr.hpp>

#include <ogr_spatialref.h>
#include <ogr_geometry.h>
#include <gdal_alg.h>
#include <ogrsf_frmts.h>

#include <cmath>

#include "raster.h"

typedef double vertex[3];

class triangle
{
public:
    triangle();
    void make_rasterized(vertex v0_in, vertex v1_in, vertex v2_in, const raster& r);


    boost::shared_ptr<raster> rasterized_triangle;
    virtual ~triangle();

    bool is_nan;

    //vertexes transformed to rasterized pixel offsets
    //x,y,z
    double v0[3];
    double v1[3];
    double v2[3];
};

