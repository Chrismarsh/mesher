#!/usr/bin/env python

from osgeo import ogr, osr
import json
import os
import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--mesh", required=True,
                        help="Mesh file for input.")

    parser.add_argument("-s", "--shapefile", required=True,
                        help="Shp file for output.")

    args = vars(parser.parse_args())

    with open(args["mesh"]) as f:
        mesh = json.load(f)

    write_shp(args["shapefile"], mesh, {}, {})

def write_shp(fname, mesh, parameter_files, initial_conditions):
    # Create the shape file to hold the triangulation

    driver = ogr.GetDriverByName("ESRI Shapefile")

    try:
        os.remove(fname)  # remove if existing
    except OSError:
        pass

    output_usm = driver.CreateDataSource(fname)

    srs_out = osr.SpatialReference()

    srs_out.ImportFromProj4(mesh['mesh']['proj4'])

    layer = output_usm.CreateLayer('mesh', srs_out, ogr.wkbPolygon)

    for key, value in parameter_files.items():
        layer.CreateField(ogr.FieldDefn(key, ogr.OFTReal))

    for key, value in initial_conditions.items():
        layer.CreateField(ogr.FieldDefn(key, ogr.OFTReal))

    for elem in range(mesh['mesh']['nelem']):
        v0 = mesh['mesh']['elem'][elem][0]
        v1 = mesh['mesh']['elem'][elem][1]
        v2 = mesh['mesh']['elem'][elem][2]

        # we need this to do the area calculation
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(mesh['mesh']['vertex'][v0][0], mesh['mesh']['vertex'][v0][1])
        ring.AddPoint(mesh['mesh']['vertex'][v1][0], mesh['mesh']['vertex'][v1][1])
        ring.AddPoint(mesh['mesh']['vertex'][v2][0], mesh['mesh']['vertex'][v2][1])
        ring.AddPoint(mesh['mesh']['vertex'][v0][0],
                      mesh['mesh']['vertex'][v0][1])  # add again to complete the ring.

        # need this for the area calculation
        tpoly = ogr.Geometry(ogr.wkbPolygon)
        tpoly.AddGeometry(ring)

        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(tpoly)

        # for this section,
        # key[0:10] -> if the name is longer, it'll have been truncated when we made the field
        # pdb.set_trace()
        for key, data in parameter_files.items():
            try:
                output = data[elem]
                # pdb.set_trace()
                feature.SetField(key[0:10], float(output))
            except TypeError as e:

                print(e)
                print('---')
                print(key)
                print(data)
                print('---')
                print(elem)
                print(output)
                print('---')
                raise(e)

        for key, data in initial_conditions.items():
            output = data[elem]
            feature.SetField(key[0:10], output)

        layer.CreateFeature(feature)
    output_usm.FlushCache()
    output_usm = None  # close file


if __name__ == "__main__":
    main()