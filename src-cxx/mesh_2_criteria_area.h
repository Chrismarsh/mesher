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
#include <limits>
#include <algorithm>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <vector>
#include <cmath>
#include "triangle.h"
#include <chrono>

// Per-process instrumentation enabled only for patch operations. No cache is
// shared between raster inputs, criteria configurations, or backend invocations.
struct MesherRasterStats { bool enabled = false; size_t evaluations = 0; double seconds = 0; };
static MesherRasterStats mesher_raster_stats;
struct MesherRasterTimer {
    std::chrono::steady_clock::time_point start;
    // Count one raster evaluation on entry. Scoping the timer to the metric includes early
    // returns for empty or NoData windows.
    MesherRasterTimer() {
        if (mesher_raster_stats.enabled) { start = std::chrono::steady_clock::now(); ++mesher_raster_stats.evaluations; }
    }
    // Add elapsed time on scope exit so every return path contributes to the same patch
    // statistics.
    ~MesherRasterTimer() {
        if (mesher_raster_stats.enabled)
            mesher_raster_stats.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
};

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
        const bool debug;

        OGRCoordinateTransformation* prj_trans;
    public:

        // Keep geometry and raster bounds together so refinement and the final face audit use
        // the same criteria. Geographic input needs an equal-area transform for area
        // comparisons.
        mesh_2_criterion_area(const double aspect_bound = 0.125,
                              const double max_area = 0,
                              const double min_area = 1,
                              const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& rasters = std::vector< std::tuple< boost::shared_ptr<raster>,double, double>> (),
                              const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& category_rasters = std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>(),
                              const std::string& error_metric = std::string(),
                              const bool is_geographic = false,
                              const bool use_weights = false,
                              const double weight_threshold = 0,
                              const bool debug = false,
                              const Geom_traits &traits = Geom_traits())
        : r(rasters),category_rasters(category_rasters),error_metric(error_metric),is_geographic(is_geographic),use_weights(use_weights),weight_threshold(weight_threshold),debug(debug)

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

        // The audit reads the same area limit used to classify faces in the refinement queue.
        inline
        double maxarea() const
        { return max_area; }

        // Update the area limit used when constructing subsequent face predicates.
        inline
        void set_area_bound(const double ab)
        { max_area = ab; }

        // CGAL expresses the shape bound as the squared sine of the smallest angle.
        inline
        double bound() const
        { return B; }

        // Accept the shape threshold in the squared-sine units stored in Quality.
        inline
        void set_bound(const double bound)
        { B = bound; }


        class Quality
        {
        public:

            // Start with an unevaluated face. The face predicate fills these fields for later
            // queue comparisons and cached quality checks.
            Quality()
            {
                _sine = 0;
                _area = 0;
                _size = 0;
                _weighted_raster_bad = false;
            };


            // Return face area in the units used for the absolute area bounds. The queue uses a
            // separate normalized size for ordering.
            const double &area() const
            { return _area; }

            // Expose the cached squared sine for direct comparison with CGAL's shape threshold.
            const double &sine() const
            { return _sine; }


            // Convert the squared sine to degrees for failure reports. Clamping keeps a
            // roundoff excursion outside the asin domain from hiding the angle.
            double min_angle_degrees() const
            {
                const double sine = std::sqrt(std::max(0.0, std::min(1.0, _sine)));
                return std::asin(sine) * 180.0 / std::acos(-1.0);
            }

        // q1<q2 means q1 is prioritised over q2
        // ( q1 == *this, q2 == q )
            bool operator<(const Quality &q) const
            {
                // CGAL stores qualities in an ordered map. Every comparison
                // must be a strict weak ordering, even when area and angle
                // disagree. The old independent comparisons could make both
                // a < b and b < a true, violating the queue's map contract.
                // Prioritize oversized faces, then the worst minimum angle,
                // following CGAL's standard size/shape refinement ordering.
                if (_size > 1 || q._size > 1)
                    return _size > q._size;
                return _sine < q._sine;
            }

            double _size; // area normalized by its maximum, for queue priority
            double _area;
            double _sine;
            // CGAL re-evaluates cached Quality when removing a face from the
            // queue. Preserve weighted rejection as well as geometry there.
            bool _weighted_raster_bad;
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
            const bool debug;
            OGRCoordinateTransformation* prj_trans;
            const Geom_traits &traits;

            boost::function<double(const typename CDT::Face_handle&, const raster&)> error_fn;
        public:

            typedef typename CDT::Point Point_2;

            // Bind the configured numeric error metric once when CGAL asks for a predicate. Its
            // per-face calls reuse the rasters and projection held by the criteria.
            Is_bad(const double aspect_bound,
                   const double area_bound,
                   const double min_area,
                   const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& r,
                   const std::vector< std::tuple< boost::shared_ptr<raster>,double, double>>& category_rasters,
                   const std::string& error_metric,
                   const bool is_geographic=false,
                   const bool use_weights=false,
                   const double weight_threshold=0,
                   const bool debug=false,
                   OGRCoordinateTransformation* prj_trans=nullptr,
                   const Geom_traits &traits = Geom_traits() )
                    : B(aspect_bound), max_area(area_bound), min_area(min_area),
                      r(r), category_rasters(category_rasters),error_metric(error_metric),is_geographic(is_geographic),
                      use_weights(use_weights),weight_threshold(weight_threshold),debug(debug),traits(traits)
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

            // Estimate how well the triangle's vertex plane represents the underlying raster
            // using RMS residuals over its valid samples.
            double rmse_tolerance(const typename CDT::Face_handle &fh, const raster &r) const
            {
                MesherRasterTimer raster_timer;
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
                if(!t.rasterized_triangle || !t.rasterized_triangle->getDs())
                {
                    // std::cout << "max_diff: invalid rasterized triangle, bailing" << std::endl;
                    return 0;
                }
                if(!t.rasterized_triangle || !t.rasterized_triangle->getDs())
                {
                    std::cout << "rmse_tolerance: invalid rasterized triangle, bailing" << std::endl;
                    return 0;
                }

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
            // Different classes at the vertices score zero. Otherwise the dominant valid class
            // fraction measures how homogeneous the triangle is.
            double categoryraster_isok(const typename CDT::Face_handle &fh, const raster& r) const
            {
                MesherRasterTimer raster_timer;
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
                if(!t.rasterized_triangle || !t.rasterized_triangle->getDs())
                {
                    std::cout << "categoryraster_isok: invalid rasterized triangle, bailing" << std::endl;
                    return true;
                }

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
            // Compare the mean of available vertex samples with the mean of valid cells in the
            // triangle. This metric measures bias across the footprint.
            double mean_tolerance(const typename CDT::Face_handle &fh, const raster &r) const
            {
                MesherRasterTimer raster_timer;
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
                if(!t.rasterized_triangle || !t.rasterized_triangle->getDs())
                {
                    std::cout << "mean_tolerance: invalid rasterized triangle, bailing" << std::endl;
                    return 0;
                }

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
            // Use the largest absolute residual against the vertex plane so a single poorly
            // represented cell can trigger refinement.
            double max_diff(const typename CDT::Face_handle &fh, const raster & r) const
            {
                MesherRasterTimer raster_timer;
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

            // CGAL can revisit cached Quality without a face handle. Apply the same area, shape
            // and raster rules during that second decision.
            CGAL::Mesh_2::Face_badness operator()(const Quality& q) const
            {
                if (q.area() > max_area)
                {
                    if (debug)
                        std::cout << "Is_bad area reject: area=" << q.area()
                                  << " max_area=" << max_area << std::endl;
                    return CGAL::Mesh_2::IMPERATIVELY_BAD; //IMPERATIVELY_BAD
                }

                if (q.sine() < this->B)
                {
                    if (debug)
                        std::cout << "Is_bad angle reject: squared_sine=" << q.sine()
                                  << " min_angle_deg=" << q.min_angle_degrees()
                                  << " bound_squared_sine=" << this->B << std::endl;
                    return CGAL::Mesh_2::BAD;
                }

                if (q.area() <= min_area )
                    return CGAL::Mesh_2::NOT_BAD;

                // Raster resolution is mandatory, like maximum area. CGAL
                // may drop ordinary BAD (shape) faces near protected acute
                // constraint clusters; doing so for raster error leaves an
                // apparently finished mesh which still violates tolerance.
                if (q._weighted_raster_bad)
                    return CGAL::Mesh_2::IMPERATIVELY_BAD;

                for(auto& itr : q._tolerance)
                {
                    if (itr.second == -1) continue; // explicitly disabled raster
                    //first == current tol
                    //second == max tol
                    if(itr.first >= itr.second) //do we violate the max tol
                    {
                        if (debug)
                            std::cout << "Is_bad tolerance reject: tol=" << itr.first
                                      << " max_tol=" << itr.second << std::endl;
                        return CGAL::Mesh_2::IMPERATIVELY_BAD;
                    }
                }

                for(auto& itr : q._category_tol)
                {
                    //first == current tol
                    //second == max tol
                    if(itr.first <= itr.second)  // is our max land cover below the required threshold to keep this triangle?
                    {
                        if (debug)
                            std::cout << "Is_bad category reject: frac=" << itr.first
                                      << " min_frac=" << itr.second << std::endl;
                        return CGAL::Mesh_2::IMPERATIVELY_BAD;
                    }
                }

                return CGAL::Mesh_2::NOT_BAD;
            }

            // Evaluate geometry first. A face below the raster refinement floor still has to
            // satisfy shape, and an oversized face can be rejected before any raster reads.
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
                if (debug)
                {
                    auto precision = std::cout.precision(std::numeric_limits<double>::max_digits10);
                    std::cout << "Is_bad triangle: ("
                              << CGAL::to_double(pa.x()) << "," << CGAL::to_double(pa.y()) << ") ("
                              << CGAL::to_double(pb.x()) << "," << CGAL::to_double(pb.y()) << ") ("
                              << CGAL::to_double(pc.x()) << "," << CGAL::to_double(pc.y()) << ")"
                              << std::endl;
                    std::cout.precision(precision);
                }

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
                    q._size = area / max_area;
                }

                // Geometry comes first. Once a well-shaped face reaches the
                // minimum area, accept it before reading any raster: further
                // numeric, categorical or weighted criteria no longer apply.
                auto current_badness = operator()(q);
                if (current_badness != CGAL::Mesh_2::NOT_BAD || q.area() <= min_area)
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

                    q._weighted_raster_bad = !(alpha >= weight_threshold);
                    return operator()(q);

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

        // Give CGAL a face predicate with the current bounds and the same borrowed raster
        // inputs used for auditing the finished mesh.
        Is_bad is_bad_object() const
        {
            return Is_bad(this->bound(), max_area, min_area, r, category_rasters, error_metric, is_geographic, use_weights, weight_threshold, debug, prj_trans, this->traits);
        }
    };
