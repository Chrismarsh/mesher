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
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <tuple>
#include <csignal>


#include "mesh.h"
#include "version.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ogrsf_frmts.h"

#include "raster.h"

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>

namespace pt = boost::property_tree;
namespace po = boost::program_options;



pt::ptree read_json(const std::string& path)
{
 // the imbue appears to fix some oddness in how the json parser parses the str.
    // http://stackoverflow.com/q/35275314/410074
    // putting it here fixed a segfault very deep in the json parse

    std::ifstream in(path);
    in.imbue(std::locale());

    std::stringstream json_file;
    if (in.is_open())
    {
        std::string line;
        while ( getline (in,line) )
        {
            json_file << line << '\n';
        }
        in.close();
    }
    else
    {
        throw std::invalid_argument("Unable to open " + path);

    }

    pt::ptree config;
    try
    {
        pt::read_json( json_file, config);
    }
    catch (pt::json_parser_error &e)
    {
      throw std::invalid_argument( "Error reading file: " + path + " on line: " + std::to_string(e.line()) + " with error: " + e.message());
    }

    return config;
}

int main(int argc, char *argv[])
{
    GDALAllRegister();

//    raise(SIGSTOP);

    std::string version = "mesher " GIT_BRANCH "/" GIT_COMMIT_HASH ;
    std::string poly_file;
    std::string interior_plgs_file;
    double max_area = 0;
    double min_area = 1;
    bool is_geographic = false;
    size_t lloyd_itr = 0;
    std::string error_metric = "rmse"; //default of RMSE
    double weight_threshold=0;// if weights are used, threshold of weighted sum that must be surpassed for triangle to be accepeted as good

    bool use_weights = false; // use the weighted generic method. If any weights are passed in (weights.size() > 0) then this will be set true and the weight methods will be used.

    // for both raster and category raster, the 3rd param are weights to apply a convex combination to the constraints.
    // If weights is emtpy, then mesher will rigorously fullfill each constraint

    //holds all rasters we perform tolerance checking on
    std::vector< std::tuple< boost::shared_ptr<raster>,double, double> > rasters;

    //these are category based rasters, e.g., landcover, so use a fractional % to figure out if we should split on the tri.
    std::vector< std::tuple< boost::shared_ptr<raster> ,double, double> > category_rasters;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "This message")
            ("version,v", "Version number")
            ("is-geographic,g", po::value<bool>(&is_geographic),"Set to true if the input data are in geographic (lat/long) format.")
            ("poly-file,p", po::value<std::string>(&poly_file),
             "PLGS file to use to bound triangulation. Same format as Triangle.")
             ("interior-plgs-file,i", po::value<std::string>(&interior_plgs_file),
             "Interior PLGS file to use to bound triangulation, e.g., rivers")

            ("raster,r", po::value<std::vector<std::string>>(), "If tolerance checking is used,"
                    "this provides a list of the rasters to provide the tolerance checking against. "
                    "Order needs to match the order the tolerances are given in.")
            ("tolerance,t", po::value<std::vector<double>>(), "Tolerances, same units as the method"
                    " Must be give in the same order as the rasters.")

            ("category-raster,R", po::value<std::vector<std::string>>(), "Optional landcover raster to conform mesh to.")
            ("category-frac,T",  po::value<std::vector<double>>(), "Fractional percent of continuous landcover required to not-split a triangle.")
            ("weight,w",  po::value<std::vector<double>>(), "Fractional weights (must all sum to 1). Must be in exactly same order as rasters and category rasters given.")
            ("area,a", po::value<double>(&max_area), "Maximum area a triangle can be. Square unit.")
            ("min-area,m", po::value<double>(&min_area), "Minimum area a triangle can be. Square unit.")
            ("lloyd,l", po::value<size_t>(&lloyd_itr), "Number of Llyod iterations.")
            ("error-metric,M", po::value<std::string>(&error_metric), "Error metric. One of: rmse, mean_tol, max_tol."
                                                         "mean_tol compares the mean triangle vertex value to the mean raster value. "
                                                        "max_tol mimics the ArcGIS TIN tolerance, and is the maximum difference between the triangle and any single raster cell.")
            ("weight-threshold,h", po::value<double>(&weight_threshold),"If weights are used, threshold of weighted sum that must be surpassed for triangle to be accepted as good.");
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        exit(1);
    }
    if (vm.count("version"))
    {
        std::cout << version << std::endl;
        exit(1);
    }

    if ( (vm.count("raster") && ! vm.count("tolerance")) ||
            (!vm.count("raster") &&  vm.count("tolerance")) )
    {
        std::cout << "Both raster and tolerance must be specified!" << std::endl;
        exit(1);
    }

    if ( (vm.count("category-raster") && ! vm.count("category-frac")) ||
         (!vm.count("category-raster") &&  vm.count("category-frac")) )
    {
        std::cout << "Both category rasters and fractions must be specified!" << std::endl;
        exit(1);
    }

    if(max_area <= 0)
    {
        std::cout << "Area must be positive" << std::endl;
        exit(1);
    }

    if(min_area <= 0 )
    {
        std::cout << "Min area must be greater than zero" << std::endl;
        exit(1);
    }

    std::vector<double> weights;
    if(vm.count("weight"))
    {
        use_weights = true;
        weights = vm["weight"].as<std::vector<double>>();
    }

    size_t mw = 0; //counter for getting the offset into the weights array

    //Ensure we have the same number of tolerances as we do raster files
    if( vm.count("raster") && vm.count("tolerance"))
    {
        auto files = vm["raster"].as<std::vector<std::string>>();
        auto tols = vm["tolerance"].as<std::vector<double>>();

        if(files.size() != tols.size())
        {
            std::cout << "Mismatched lengths in rasters and tolerances. Must be equal."<<std::endl;
            exit(1);
        }

//        if(use_weights && weights.size() != tols.size())
//        {
//            std::cout << "Mismatched lengths in tolerances and weights. Must be equal."<<std::endl;
//            exit(1);
//        }

        for(int i = 0 ; i< tols.size();++i)
        {
            auto r = boost::make_shared<raster>();
            r->open(files.at(i));

            double w=1;
            if (weights.size() > 0) w = weights.at(mw++) ; // if there are no weights, just insert a dummy var

            if(w > 1)
            {
                std::cout << "Weights must be <=1!"<<std::endl;
                exit(1);
            }

            rasters.push_back(std::make_tuple( r,tols.at(i), w ));
        }
    }

    //ensure we have the same number of category rasters as category fractions
    if( vm.count("category-raster") && vm.count("category-frac"))
    {
        auto files = vm["category-raster"].as<std::vector<std::string>>();
        auto tols = vm["category-frac"].as<std::vector<double>>();

        if(files.size() != tols.size())
        {
            std::cout << "Mismatched lengths in category-raster and category-frac. Must be equale."<<std::endl;
            exit(1);
        }
//        if(use_weights && weights.size() != tols.size())
//        {
//            std::cout << "Mismatched lengths in category-fra and weights. Must be equal."<<std::endl;
//            exit(1);
//        }

        for(int i = 0 ; i< tols.size();++i)
        {
            auto r = boost::make_shared<raster>();
            r->open(files.at(i));

            if(tols.at(i) > 1)
            {
                std::cout << "Fractional percentage required for category-frac!"<<std::endl;
                exit(1);
            }

            double w=1;
            if (weights.size() > 0) w = weights.at(mw++) ; // if there are no weights, just insert a dummy var
            
            if(w > 1)
            {
                std::cout << "Weights must be <=1!"<<std::endl;
                exit(1);
            }

            category_rasters.push_back(std::make_tuple( r,tols.at(i),w ));
        }
    }

    CDT cdt; //main mesh data structure

    //read exterior PLGS
    std::ifstream infile(poly_file);
    if(!infile)
    {
        std::cout << "Failed to open poly file" << std::endl;
        exit(1);
    }

    boost::filesystem::path p(poly_file);
    auto path = p.parent_path();

    std::string line;
    //header
    std::getline(infile, line);
    std::istringstream iss(line);
    int row,col;
    iss >> row >> col;
    std::cout << std::setprecision(20);
    std::vector<Vertex_handle> vertex;
    if (col != 2)
    {
        std::cout << "Wrong number of columns" <<std::endl;
        exit(1);
    }

    //Build the PLGS that constrains the triangulation
    for(int i = 0 ; i < row -1 ; i++)
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        double id,x,y;
        if (!(iss >> id >> x >> y))
        {
            std::cout << "Error on line " << i << std::endl;
            exit(1);
        }

        Vertex_handle vh = cdt.insert(Point(x,y));
        vh->info()=i;
        vertex.push_back(vh);
    }

    //skip last line of vertexes. this is a duplicate, so don't need it (Triangle does though and these files are Triangle compatible)
    std::getline(infile, line);

    //empty line
    std::getline(infile, line);

    //header
    std::getline(infile, line);
    std::istringstream header2(line);
    double special;
    header2 >> row >> special;

    //ignoring the last 2 items as they are for Triangle, but we aren't using it here
    for(int i=0;i<row-2; i++)
    {
        std::getline(infile, line);

        int rowid,i0,i1;
        std::istringstream iss(line);
        if(! (iss >> rowid >> i0 >> i1))
        {
            std::cout << "Error in section 2 on line " << rowid << std::endl;
            exit(1);
        }
        //1 indexed,fix
        --i0;
        --i1;
        Vertex_handle v0,v1;
        v0 = vertex.at(i0);
        v1 = vertex.at(i1);

        cdt.insert_constraint(v0,v1); //set PLGS constraint
    }

    //and lastly, close the loop
    Vertex_handle v0,v1;
    v0 = vertex.at(row-2);
    v1 = vertex.at(0);
    cdt.insert_constraint(v0,v1);

    //do interior plgs
    auto interior_plgs = read_json(interior_plgs_file);
    std::cout << "Reading interior PLGS file " << interior_plgs_file << std::endl;

    //iterate over all the features
    for(auto& feat : interior_plgs.get_child("features") )
    {
        //this will hold the vertexes we make from all the read in points. Then they will be joined in a line to make the constriant
        std::vector<Vertex_handle> interior_vertexes;

        for(auto& coords : feat.second.get_child("geometry.coordinates"))
        {
            std::vector<double> v;

            for(auto& c: coords.second)
            {
                v.push_back( boost::lexical_cast<double>(c.second.data()));
            }
            Vertex_handle vh = cdt.insert(Point(v[0],v[1]));
            interior_vertexes.push_back(vh);
//            std::cout <<  v[0]<<","<<v[1] << std::endl;

        }

        //go up to end - 1 as we will not make a cycle, just a linear feature
        for(size_t i = 0; i < interior_vertexes.size() - 1; i++)
        {
//            std::cout << i << std::endl;
            Vertex_handle v0,v1;
            v0 = interior_vertexes.at(i);
            v1 = interior_vertexes.at(i+1);


            try
            {
                cdt.insert_constraint(v0,v1);
            }catch(std::runtime_error& e)
            {
                std::cout << "[ Warning ]: The following point lies outside of the outer PLGS. It will be ignored." << std::endl;
                std::cout << "\tv0 x=" << v0->point().x() << ",y=" <<v0->point().y() << std::endl;
                std::cout << "\tv1 x=" << v1->point().x() << ",y=" <<v1->point().y() << std::endl;
            }
        }

    }


    std::cout << "Number of input PLGS vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Meshing the triangulation..." << std::endl;
    CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125 /*internal angle*/,max_area,min_area,rasters,category_rasters,error_metric,is_geographic,use_weights,weight_threshold));

    //run lloyd optimizations if required.
    //if run, 100 is a good pick
    //http://doc.cgal.org/latest/Mesh_2/#title9
    if(lloyd_itr > 0)
    {
        std::cout << "Running " << lloyd_itr << " Lloyd iterations...";
        CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = lloyd_itr);
        std::cout << " done." << std::endl;
    }

    auto nodefilepath = path; //eg PLGSwolf_lidar1.1.node
    nodefilepath /= p.filename().replace_extension(".1.node");
    std::ofstream nodefile( nodefilepath.string());

    //header
    nodefile << cdt.number_of_vertices() << " 2 0 1" << std::endl;
    nodefile << std::setprecision(15);

    size_t mesh_vertex_i=1;

    for(auto itr = cdt.finite_vertices_begin(); itr != cdt.finite_vertices_end(); ++itr)
    {

        itr->info() = mesh_vertex_i;

        nodefile << mesh_vertex_i << "   " << itr->point().x() << "   " << itr->point().y() << "   0" << std::endl;
        ++mesh_vertex_i;
    }


    nodefile.close();


    int elem_i=1; //1 indexing

    for(auto itr = cdt.finite_faces_begin(); itr != cdt.finite_faces_end(); itr++ )
    {
        if(itr->is_in_domain())
        {
            itr->id = elem_i;
            ++elem_i;
        }

    }

    //i has total number of faces that are in domain
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
    std::cout << "Number of triangles: " << elem_i << std::endl;

    //niehgbour file
    auto neighfilepath = path; //eg PLGSwolf_lidar1.1.neigh
    neighfilepath /= p.filename().replace_extension(".1.neigh");
    std::ofstream neighfile(neighfilepath.string());
    neighfile << elem_i-1 << " 3" << std::endl; //elem_i is 1 indexing, convert to 0

    //element file, defines the triangles
    auto elefilepath = path; //eg "PLGSwolf_lidar1.1.ele"
    elefilepath /= p.filename().replace_extension(".1.ele");
    std::ofstream elemfile(elefilepath.string());
    elemfile << elem_i-1 << " 3 0" << std::endl;


    int i=1;

    for(auto itr = cdt.finite_faces_begin(); itr != cdt.finite_faces_end(); itr++ )
    {
        if(itr->is_in_domain())
        {
            size_t v0 = itr->vertex(0)->info();
            size_t v1 = itr->vertex(1)->info();
            size_t v2 = itr->vertex(2)->info();

            elemfile << i << "    " << v0 << "    " << v1 << "    "<< v2 << std::endl;

            auto n0 = itr->neighbor(0);
            auto n1 = itr->neighbor(1);
            auto n2 = itr->neighbor(2);

            neighfile << i <<  "  " << n0->id << "  " << n1->id  <<"  "<< n2->id << std::endl;
          ++i;
        }

    }

    elemfile.close();
    neighfile.close();
    return 0;
}

