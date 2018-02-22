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


//Defining this forces a custom exception to be raised within the CGAL code
// Constrained_triangulat_2.h l 846
// this allows for catching constraints that intersect the domain boundaries and are invalid
#define CGAL_CT2_WANTS_TO_HAVE_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS
#define CGAL_CDT2_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS throw std::runtime_error("Point results in intersection contraints");


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include "mesh_2_criteria_area.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//Set up a vertex that can hold extra information
template < typename Info_, typename GT,
        typename Vb = CGAL::Delaunay_mesh_vertex_base_2<GT> >
class Delaunay_mesh_vertex_base_with_info_2
        : public Vb
{
    Info_ _info;

public:
    typedef typename Vb::Face_handle                   Face_handle;
    typedef typename Vb::Point                         Point;
    typedef Info_                                      Info;

    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
        typedef Delaunay_mesh_vertex_base_with_info_2<Info, GT, Vb2>   Other;
    };

    Delaunay_mesh_vertex_base_with_info_2()
            : Vb() {}

    Delaunay_mesh_vertex_base_with_info_2(const Point & p)
            : Vb(p) {}

    Delaunay_mesh_vertex_base_with_info_2(const Point & p, Face_handle c)
            : Vb(p, c) {}

    Delaunay_mesh_vertex_base_with_info_2(Face_handle c)
            : Vb(c) {}

    const Info& info() const { return _info; }
    Info&       info()       { return _info; }
};

typedef Delaunay_mesh_vertex_base_with_info_2<size_t, K> Vb;

//setup a mesh face class that can hold additional information
template<class Gt, class Fb = CGAL::Delaunay_mesh_face_base_2<Gt> >
class Delaunay_mesh_face_base_info_2 : public Fb
{
public:
    typedef Gt Geom_traits;
    typedef typename Fb::Vertex_handle Vertex_handle;
    typedef typename Fb::Face_handle Face_handle;

    int id;

    template<typename TDS2>
    struct Rebind_TDS
    {
        typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef Delaunay_mesh_face_base_info_2<Gt, Fb2> Other;
    };

    Delaunay_mesh_face_base_info_2() : Fb()
    {
        id = 0;
    }

    Delaunay_mesh_face_base_info_2(Vertex_handle v0,
                                   Vertex_handle v1,
                                   Vertex_handle v2)
            : Fb(v0, v1, v2)
    {
        id = 0;
    }


    Delaunay_mesh_face_base_info_2(Vertex_handle v0,
                                   Vertex_handle v1,
                                   Vertex_handle v2,
                                   Face_handle n0,
                                   Face_handle n1,
                                   Face_handle n2)
            : Fb(v0, v1, v2, n0, n1, n2)
    {
        id = 0;
    }


};

typedef Delaunay_mesh_face_base_info_2<K> Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

typedef mesh_2_criterion_area<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef CDT::Face_circulator Face_circulator;
typedef CDT::Point Point;
