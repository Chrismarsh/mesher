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

#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>
#include <algorithm>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <vector>
#include <cmath>
#include "triangle.h"

    template<class CDT>
    class mesh_2_criterion_area
    {
    protected:
        typedef typename CDT::Geom_traits Geom_traits;
        double max_area;
        double min_area;
        double B;
        Geom_traits traits;
        const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& r;
        const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& category_rasters;
        const std::string& error_metric;
        const bool is_geographic;
        const bool use_weights;
        const double weight_threshold;

        OGRCoordinateTransformation* prj_trans;
    public:

        mesh_2_criterion_area(const double aspect_bound = 0.125,
                              const double max_area = 0,
                              const double min_area = 1,
                              const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& rasters = std::vector< std::tuple< boost::shared_ptr<raster>,double, double>> (),
                              const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& category_rasters = std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>(),
                              const std::string& error_metric = std::string(),
                              const bool is_geographic = false,
                              const bool use_weights = false,
                              const double weight_threshold = 0,
                              const Geom_traits &traits = Geom_traits())
        : r(rasters),category_rasters(category_rasters),error_metric(error_metric),is_geographic(is_geographic),use_weights(use_weights),weight_threshold(weight_threshold)

        {
            this->max_area = max_area;
            this->min_area = min_area;
            B = aspect_bound;
            this->traits=traits;
            prj_trans = nullptr;

            if(is_geographic)
            {
                const char* wkt =  std::get<0>(r.at(0))->getDs()->GetProjectionRef();
                OGRSpatialReference srs;

#if GDAL_VERSION_MAJOR < 3
              char* srs_wkt_nonconst = const_cast<char*>(wkt);
              srs.importFromWkt(&srs_wkt_nonconst);
#else
              srs.importFromWkt(&wkt);
#endif

                const char* out_wkt = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ";
                OGRSpatialReference  srs_out;
                srs_out.importFromProj4(out_wkt);

                prj_trans = OGRCreateCoordinateTransformation( &srs,
                                                               &srs_out );

                if(!prj_trans)
                {
                    std::cout << "Unable to create geographic transform" << std::endl;
                    exit(1);
                }
            }
        }

        inline
        double maxarea() const
        { return max_area; }

        inline
        void set_area_bound(const double ab)
        { max_area = ab; }

        inline
        double bound() const
        { return B; }

        inline
        void set_bound(const double bound)
        { B = bound; }


        class Quality
        {
        public:

            Quality()
            {
                _sine = 0;
                _area = 0;
            };


            const double &area() const
            { return _area; }

            const double &sine() const
            { return _sine; }


        // q1<q2 means q1 is prioritised over q2
        // ( q1 == *this, q2 == q )
            bool operator<(const Quality &q) const
            {
                if ( this->_area < q._area)
                    return true;
                else
                {
                    if( this->_sine < q._sine)
                        return true;
                    else
                    {
                        if(_tolerance.size() != 0 && q._tolerance.size() == 0)
                            return false; //they're empty, we aren't, we should check their tolerance
                        if(_tolerance.size() == 0 && q._tolerance.size() != 0)
                            return true; //we're empty, they aren't, we should check our tolerance

                        if(_tolerance.size() != 0 && q._tolerance.size() != 0)
                        {
                            //compare who has the the greatest current error
                            auto my_max_tol = std::max_element(_tolerance.begin(), _tolerance.end(),
                                                               [] (std::pair< double,double> const& lhs, std::pair< double,double> const& rhs)
                                                               {return lhs.first < rhs.first;});

                            auto their_max_tol = std::max_element(q._tolerance.begin(), q._tolerance.end(),
                                                                  [] (std::pair< double,double> const& lhs, std::pair< double,double> const& rhs)
                                                                  {return lhs.first < rhs.first;});

                            if(my_max_tol->first > their_max_tol->first )
                                return true; //remember we want to be preferential if our max tolerance is worse
                        }

                        return (sine() < q.sine()); //check angles as tie breaker
                    }
                }
            }

            double _area;
            double _sine;
            std::vector< std::pair<double,double> > _tolerance; // first = tol, seond = max tol

            std::vector<std::pair<double,double>> _category_tol; //first = current percent, second = threshold percent

        };

        class Is_bad
        {

        protected:
            const double B;
            const double max_area;
            const double min_area;
            const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& r;
            const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& category_rasters;
            const std::string& error_metric;
            const bool is_geographic;
            const bool use_weights;
            const double weight_threshold;
            OGRCoordinateTransformation* prj_trans;
            const Geom_traits &traits;

            boost::function<double(const typename CDT::Face_handle&, const raster&)> error_fn;
        public:

            typedef typename CDT::Point Point_2;

            Is_bad(const double aspect_bound,
                   const double area_bound,
                   const double min_area,
                   const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& r,
                   const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& category_rasters,
                   const std::string& error_metric,
                   const bool is_geographic=false,
                   const bool use_weights=false,
                   const double weight_threshold=0,
                   OGRCoordinateTransformation* prj_trans=nullptr,
                   const Geom_traits &traits = Geom_traits() )
                    : B(aspect_bound), max_area(area_bound), min_area(min_area),
                      r(r), category_rasters(category_rasters),error_metric(error_metric),is_geographic(is_geographic),
                      use_weights(use_weights),weight_threshold(weight_threshold),traits(traits)
            {
                this->prj_trans = prj_trans;
                if(!prj_trans && is_geographic)
                {
                    std::cout << "NULL prj_trans for geographic data!" << std::endl;
                    exit(1);
                }
                if(error_metric == "rmse" || error_metric == "")
                    error_fn = boost::bind(&Is_bad::rmse_tolerance,this,_1,_2);
                else if(error_metric == "mean_tol")
                    error_fn = boost::bind(&Is_bad::mean_tolerance,this,_1,_2);
                else if(error_metric == "max_tol")
                    error_fn = boost::bind(&Is_bad::max_diff,this,_1,_2);
                else
                {
                    std::cout << "Unknown error function selected" << std::endl;
                    exit(1);
                }
            }

            double rmse_tolerance(const typename CDT::Face_handle &fh, const raster &r) const
            {
                triangle t;
                vertex _v0;
                vertex _v1;
                vertex _v2;

                _v0[0] = fh->vertex(0)->point().x();
                _v0[1] = fh->vertex(0)->point().y();

                _v1[0] = fh->vertex(1)->point().x();
                _v1[1] = fh->vertex(1)->point().y();

                _v2[0] = fh->vertex(2)->point().x();
                _v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(_v0, _v1, _v2, r);
                if(t.is_nan)
                    return 0; //bail

                auto pxpy = t.rasterized_triangle->xy_to_pxpy(t.v0[0],t.v0[1]);
                t.v0[0] = pxpy.first;
                t.v0[1] = pxpy.second;

                pxpy = t.rasterized_triangle->xy_to_pxpy(t.v1[0],t.v1[1]);
                t.v1[0] = pxpy.first;
                t.v1[1] = pxpy.second;


                pxpy = t.rasterized_triangle->xy_to_pxpy(t.v2[0],t.v2[1]);
                t.v2[0] = pxpy.first;
                t.v2[1] = pxpy.second;


                //Xs or Y s are collinear. This seems to happen if we've set a min triangle area less than the pixel resolution of our dem
                //when this happens, all the points end up collinear in 1 axis, and anything that depends on this
                if(  (t.v0[0] == 0 && t.v1[0] == 0 && t.v2[0] == 0) || (t.v0[1] == 0 && t.v1[1] == 0 && t.v2[1] == 0) )
                {
                    return 0; //bail out nothing we can do
                }

                //create the vectors veco0 and
                double u1,u2,u3;
                double v1,v2,v3;
                double o1,o2,o3; //origin of tri

                o1 = t.v0[0];
                o2 = t.v0[1];
                o3 = t.v0[2];

                //if we have only nan z values, bail.
                if(std::isnan(t.v0[2]) || std::isnan(t.v1[2]) || std::isnan(t.v2[2]))
                {
                    t.is_nan=true;
                    return 0;
                }

                //following http://www.had2know.com/academics/equation-plane-through-3-points.html
                //create the two vectors
                u1 = t.v1[0] - o1;
                u2 = t.v1[1] - o2;
                u3 = t.v1[2] - o3;

                v1 = t.v2[0] - o1;
                v2 = t.v2[1] - o2;
                v3 = t.v2[2] - o3;

                //calculate the normal vector via cross product
                double a, b, c;
                a = u2 * v3 - v2 * u3;
                b = v1 * u3 - u1 * v3;
                c = u1 * v2 - v1 * u2;


                //solve for d
                double d =a*o1 + b*o2 + c*o3;

                double rmse = 0;
                double n = 0;

                for (int y = 0; y < t.rasterized_triangle->getDs()->GetRasterYSize(); y++)
                {
                    for (int x = 0; x < t.rasterized_triangle->getDs()->GetRasterXSize(); x++)
                    {
                        double value = t.rasterized_triangle->getpXpY(x, y);
                        if (!std::isnan(value))
                        {
                            double z = -(a*x+b*y-d)/c; //plane eqn solved for z. allows us to predict z values via x,y coords
                            double diff = z - value;

                            rmse += diff * diff;
                            n++;
                        }
                    }
                }

                //bail, somehow we have no raster cells under our triangle.
                if (n == 0.)
                    return 0;

                rmse /= n;

                rmse = sqrt(rmse);

                return rmse;

            }
            double categoryraster_isok(const typename CDT::Face_handle &fh, const raster& r) const
            {
                triangle t;
                vertex _v0;
                vertex _v1;
                vertex _v2;

                _v0[0] = fh->vertex(0)->point().x();
                _v0[1] = fh->vertex(0)->point().y();

                _v1[0] = fh->vertex(1)->point().x();
                _v1[1] = fh->vertex(1)->point().y();

                _v2[0] = fh->vertex(2)->point().x();
                _v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(_v0, _v1, _v2, r);


                if(t.is_nan)
                    return true; //bail

                if( t.v0[2] != t.v1[2] ||
                    t.v0[2] != t.v2[2] ||
                    t.v1[2] != t.v2[2])
                    return false;

                double n=0;//total tri
                std::map<int,int> lc; //holds the count of each category (e.g., landcover type)
                for (int y = 0; y < t.rasterized_triangle->getDs()->GetRasterYSize(); y++)
                {
                    for (int x = 0; x < t.rasterized_triangle->getDs()->GetRasterXSize(); x++)
                    {
                        double value = t.rasterized_triangle->getpXpY(x, y);
                        if (!std::isnan(value))
                        {
                            int key = (int)value;
                            if(lc.find(key) != lc.end())
                                lc[key]++;
                            else
                                lc[key] = 1;
                            n++;
                        }
                    }
                }

                if (n==0)
                    return true; //bail

                bool isok=false;

                std::vector<double> percents;

                for(auto& itr : lc)
                {
                    double pCover =  itr.second / n;
                    percents.push_back(pCover);

                }
                return *std::max_element(percents.begin(),percents.end()); // the most dominate landcover

            }
            double mean_tolerance(const typename CDT::Face_handle &fh, const raster &r) const
            {
                triangle t;
                vertex v0;
                vertex v1;
                vertex v2;

                v0[0] = fh->vertex(0)->point().x();
                v0[1] = fh->vertex(0)->point().y();

                v1[0] = fh->vertex(1)->point().x();
                v1[1] = fh->vertex(1)->point().y();

                v2[0] = fh->vertex(2)->point().x();
                v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(v0, v1, v2, r);
                if(t.is_nan)
                    return 0; //bail

                // Initialize triangle mean elevation (m)
                double triangle_z_mean = 0;
                // Initialize count of triangle vertices found not-nan
                double tri_count = 0;

                // Check if elevations of vertex is not-nan
                if (!std::isnan(t.v0[2]))
                {
                    triangle_z_mean += t.v0[2];
                    ++tri_count;
                }

                if (!std::isnan(t.v1[2]))
                {
                    triangle_z_mean += t.v1[2];
                    ++tri_count;
                }

                if (!std::isnan(t.v2[2]))
                {
                    triangle_z_mean += t.v2[2];
                    ++tri_count;
                }


                //planar mean elevation
                triangle_z_mean /= tri_count;

                // If no non-nan vertices were found, return zero
                if (tri_count == 0)
                    return 0;

                // Initialize sum and count of grid cell elevations
                double sum = 0;
                double count = 0;

                for (int i = 0; i < t.rasterized_triangle->getDs()->GetRasterYSize(); i++)
                {
                    for (int j = 0; j < t.rasterized_triangle->getDs()->GetRasterXSize(); j++)
                    {
                        double value = t.rasterized_triangle->getpXpY(j, i);
                        if (!std::isnan(value))
                        {
                            sum += value;
                            ++count;
                        }
                    }
                }

                // If none were found return zero
                if (count == 0)
                    return 0;

                // Take mean of all grid cell elevations found
                double mean = sum / count;

                // Take difference of means
                double diff = fabs(triangle_z_mean - mean);

                return diff;

            }
            double max_diff(const typename CDT::Face_handle &fh, const raster & r) const
            {
                triangle t;
                vertex _v0;
                vertex _v1;
                vertex _v2;

                _v0[0] = fh->vertex(0)->point().x();
                _v0[1] = fh->vertex(0)->point().y();

                _v1[0] = fh->vertex(1)->point().x();
                _v1[1] = fh->vertex(1)->point().y();

                _v2[0] = fh->vertex(2)->point().x();
                _v2[1] = fh->vertex(2)->point().y();

                t.make_rasterized(_v0, _v1, _v2, r);
                if(t.is_nan)
                    return 0; //bail

                auto pxpy = t.rasterized_triangle->xy_to_pxpy(t.v0[0],t.v0[1]);
                t.v0[0] = pxpy.first;
                t.v0[1] = pxpy.second;

                pxpy = t.rasterized_triangle->xy_to_pxpy(t.v1[0],t.v1[1]);
                t.v1[0] = pxpy.first;
                t.v1[1] = pxpy.second;


                pxpy = t.rasterized_triangle->xy_to_pxpy(t.v2[0],t.v2[1]);
                t.v2[0] = pxpy.first;
                t.v2[1] = pxpy.second;


                //Xs or Y s are collinear. This seems to happen if we've set a min triangle area less than the pixel resolution of our dem
                //when this happens, all the points end up collinear in 1 axis, and anything that depends on this
                if(  (t.v0[0] == 0 && t.v1[0] == 0 && t.v2[0] == 0) || (t.v0[1] == 0 && t.v1[1] == 0 && t.v2[1] == 0) )
                {
                    return 0; //bail out nothing we can do
                }

                //create the vectors veco0 and
                double u1,u2,u3;
                double v1,v2,v3;
                double o1,o2,o3; //origin of tri

                o1 = t.v0[0];
                o2 = t.v0[1];
                o3 = t.v0[2];

                //if we have only nan z values, bail.
                if(std::isnan(t.v0[2]) || std::isnan(t.v1[2]) || std::isnan(t.v2[2]))
                {
                    t.is_nan=true;
                    return 0;
                }

                //following http://www.had2know.com/academics/equation-plane-through-3-points.html
                //create the two vectors
                u1 = t.v1[0] - o1;
                u2 = t.v1[1] - o2;
                u3 = t.v1[2] - o3;

                v1 = t.v2[0] - o1;
                v2 = t.v2[1] - o2;
                v3 = t.v2[2] - o3;

                //calculate the normal vector via cross product
                double a, b, c;
                a = u2 * v3 - v2 * u3;
                b = v1 * u3 - u1 * v3;
                c = u1 * v2 - v1 * u2;


                //solve for d
                double d =a*o1 + b*o2 + c*o3;

                double max_diff = -1;
                double n = 0;


                for (int y = 0; y < t.rasterized_triangle->getDs()->GetRasterYSize(); y++)
                {
                    for (int x = 0; x < t.rasterized_triangle->getDs()->GetRasterXSize(); x++)
                    {
                        double value = t.rasterized_triangle->getpXpY(x, y);
                        if (!std::isnan(value))
                        {
                            double z = -(a*x+b*y-d)/c; //plane eqn solved for z. allows us to predict z values via x,y coords
                            double diff = fabs(z - value);

                            if(diff > max_diff )
                                max_diff = diff;

                            n++;
                        }
                    }
                }

                //bail, somehow we have no raster cells under our triangle.
                if (n == 0.)
                    return 0;

                return max_diff;

            }

            CGAL::Mesh_2::Face_badness operator()(const Quality q) const
            {
                if (q.area() > max_area)
                    return CGAL::Mesh_2::IMPERATIVELY_BAD; //IMPERATIVELY_BAD

                if (q.sine() < this->B)
                    return CGAL::Mesh_2::BAD;

                if (q.area() <= min_area )
                    return CGAL::Mesh_2::NOT_BAD;

                for(auto& itr : q._tolerance)
                {
                    //first == current tol
                    //second == max tol
                    if(itr.first >= itr.second) //do we violate the max tol
                        return CGAL::Mesh_2::BAD;
                }

                for(auto& itr : q._category_tol)
                {
                    //first == current tol
                    //second == max tol
                    if(itr.first <= itr.second)  // is our max land cover below the required threshold to keep this triangle?
                        return CGAL::Mesh_2::BAD;
                }

                return CGAL::Mesh_2::NOT_BAD;
            }

            CGAL::Mesh_2::Face_badness operator()(const typename CDT::Face_handle &fh,
                                                  Quality &q) const
            {
                typedef typename CDT::Geom_traits Geom_traits;
                typedef typename Geom_traits::Compute_area_2 Compute_area_2;
                typedef typename Geom_traits::Compute_squared_distance_2
                        Compute_squared_distance_2;

                Geom_traits traits;

                Compute_area_2 area_2 = traits.compute_area_2_object();
                Compute_squared_distance_2 squared_distance = traits.compute_squared_distance_2_object();

                const Point_2 &pa = fh->vertex(0)->point();
                const Point_2 &pb = fh->vertex(1)->point();
                const Point_2 &pc = fh->vertex(2)->point();

                double a = CGAL::to_double(squared_distance(pb, pc));
                double b = CGAL::to_double(squared_distance(pc, pa));
                double c = CGAL::to_double(squared_distance(pa, pb));

                //input rasters are in lat/long (ie., geographic) so we must compute the area differently.
                double area = 0;
                if(is_geographic)
                {
                    const char* srs_wkt = std::get<0>(r.at(0))->getDs()->GetProjectionRef();
                    OGRSpatialReference srs;
#if GDAL_VERSION_MAJOR < 3
                    char* srs_wkt_nonconst = const_cast<char*>(srs_wkt);
                    srs.importFromWkt(&srs_wkt_nonconst);
#else
                    srs.importFromWkt(&srs_wkt);
#endif

                    vertex v0;
                    vertex v1;
                    vertex v2;

                    v0[0] = fh->vertex(0)->point().x();
                    v0[1] = fh->vertex(0)->point().y();

                    v1[0] = fh->vertex(1)->point().x();
                    v1[1] = fh->vertex(1)->point().y();

                    v2[0] = fh->vertex(2)->point().x();
                    v2[1] = fh->vertex(2)->point().y();

                    OGRLinearRing ring;
                    ring.addPoint(v0[0],v0[1]);
                    ring.addPoint(v1[0],v1[1]);
                    ring.addPoint(v2[0],v2[1]);
                    ring.addPoint(v0[0],v0[1]); //close it

                    OGRPolygon poly;
                    poly.addRing(&ring);

                    OGRSpatialReference  srs_out;
                    srs_out.importFromProj4("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ");

                    auto poCT = OGRCreateCoordinateTransformation( &srs,
                                                                   &srs_out );

                    poly.transform(poCT);

                    area = poly.get_Area();

                }
                else
                {
                   area = CGAL::to_double(area_2(pa, pb, pc));
                }


                // area = 4 * area^2(triangle)
                double area2 = 2.0 * CGAL::to_double(area_2(pa, pb, pc));
                area2 = area2 * area2;

                if (a < b) if (a < c)
                    q._sine = area2 / (b * c);
                else
                    q._sine = area2 / (a * b);
                else if (b < c)
                    q._sine = area2 / (a * c);
                else
                    q._sine = area2 / (a * b);

                if (max_area != 0)
                {
                    q._area = area;
                }

                //if our triangle already fails either: area or angle, bail now.
                auto current_badness = operator()(q);
                if (current_badness != CGAL::Mesh_2::NOT_BAD)
                    return current_badness;

                if(use_weights)
                {
                    //total weighted score of this triangle
                    double alpha = 0;

                    //do numeric rasters
                    // <0> = raster
                    // <1> = tol
                    // <2> = weight
                    for(auto& itr : r)
                    {
                        //if tolerance == -1, skip the tolerance tests because we want to make a uniform mesh.
                        if( std::get<1>(itr) != -1)
                        {
                            auto t = error_fn(fh,*( std::get<0>(itr)));

                            // Don't push back, we don't want to be checking these, only want to check area and shape, we will manually do the tolerance
                            // q._tolerance.push_back(std::make_pair(t, std::get<1>(itr)));

                            // even with the weights there are a few things which we must fullfill
                            if (q.area() > max_area)
                                return CGAL::Mesh_2::IMPERATIVELY_BAD;

                            //bad angles
                            if (q.sine() < this->B)
                                return CGAL::Mesh_2::BAD;

                            //we've gone to small, we need out
                            if (q.area() <= min_area )
                                return CGAL::Mesh_2::NOT_BAD;


                            if( t <= std::get<1>(itr))
                                current_badness = CGAL::Mesh_2::NOT_BAD;
                            else
                                current_badness = CGAL::Mesh_2::BAD;


                            double weight =  std::get<2>(itr);
                            weight *= (current_badness == CGAL::Mesh_2::NOT_BAD ? 1 : 0);
                            alpha += weight;
                        }
                        else
                        {
                            //give it a good score
                            alpha += std::get<2>(itr);
                        }

                    }

                    for(auto& itr : category_rasters)
                    {
                        auto t = categoryraster_isok(fh,*(std::get<0>(itr)));
//                      q._category_tol.push_back(std::make_pair(t,std::get<1>(itr)));

                        // even with the weights there are a few things which we must fullfill
                        if (q.area() > max_area)
                            return CGAL::Mesh_2::IMPERATIVELY_BAD;

                        //bad angles
                        if (q.sine() < this->B)
                            return CGAL::Mesh_2::BAD;

                        //we've gone to small, we need out
                        if (q.area() <= min_area )
                            return CGAL::Mesh_2::NOT_BAD;

                        if( t >= std::get<1>(itr))
                            current_badness = CGAL::Mesh_2::NOT_BAD;
                        else
                            current_badness = CGAL::Mesh_2::BAD;

                        double weight =  std::get<2>(itr);
                        weight *= (current_badness == CGAL::Mesh_2::NOT_BAD ? 1 : 0);
                        alpha += weight;
                    }

                    if (alpha >= weight_threshold)
                        return CGAL::Mesh_2::NOT_BAD;
                    else
                        return CGAL::Mesh_2::BAD;

                } else //this rigorously ensures each tolerance is met
                {
                    //do numeric rasters
                    // <0> = raster
                    // <1> = tol
                    // <2> = weight
                    for(auto& itr : r)
                    {
                        //if tolerance == -1, skip the tolerance tests because we want to make a uniform mesh.
                        if( std::get<1>(itr) != -1)
                        {
                            auto t = error_fn(fh,*( std::get<0>(itr)));
                            q._tolerance.push_back(std::make_pair(t, std::get<1>(itr)));

                            current_badness = operator()(q);
                            if (current_badness != CGAL::Mesh_2::NOT_BAD)
                                return current_badness;
                        }
                        else
                        {
                            //fake it passing by giving it a 0 rmse
                            q._tolerance.push_back(std::make_pair(0,std::get<1>(itr)));
                        }

                    }

                    for(auto& itr : category_rasters)
                    {
                        auto t = categoryraster_isok(fh,*(std::get<0>(itr)));
                        q._category_tol.push_back(std::make_pair(t,std::get<1>(itr)));

                        current_badness = operator()(q);
                        if (current_badness != CGAL::Mesh_2::NOT_BAD)
                            return current_badness;
                    }
                    return operator()(q);
                }
            }
        };

        Is_bad is_bad_object() const
        {
            return Is_bad(this->bound(), max_area, min_area, r, category_rasters, error_metric, is_geographic, use_weights, weight_threshold, prj_trans, this->traits);
        }
    };
