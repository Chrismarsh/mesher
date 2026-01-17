from osgeo import ogr


def extract_polygons(geom):
    if geom is None:
        return []

    gtype = geom.GetGeometryType()
    if gtype in (ogr.wkbPolygon, ogr.wkbPolygon25D):
        return [geom.Clone()]
    if gtype in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D,
                 ogr.wkbGeometryCollection, ogr.wkbGeometryCollection25D):
        polys = []
        for i in range(geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            polys.extend(extract_polygons(g))
        return polys
    return []


def normalize_polygon(geom):
    if geom is None:
        return None
    geom = geom.MakeValid()
    geom = geom.Buffer(0)
    polys = extract_polygons(geom)
    if not polys:
        return None
    return polys[0]


def bbox_to_polygon_geom(bbox):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(bbox[0], bbox[1])
    ring.AddPoint(bbox[2], bbox[1])
    ring.AddPoint(bbox[2], bbox[3])
    ring.AddPoint(bbox[0], bbox[3])
    ring.AddPoint(bbox[0], bbox[1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def load_polygon_geom(shp_path):
    ds = ogr.Open(shp_path)
    if ds is None:
        raise RuntimeError(f'Unable to open polygon shapefile: {shp_path}')
    layer = ds.GetLayer(0)
    geom_union = None
    for feat in layer:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        geom = geom.Clone()
        geom = geom.MakeValid()
        if geom_union is None:
            geom_union = geom
        else:
            geom_union = geom_union.Union(geom)
    ds = None
    geom_union = normalize_polygon(geom_union)
    if geom_union is None:
        raise RuntimeError(f'Polygon shapefile produced no valid polygons: {shp_path}')
    return geom_union


def iter_linestrings(geom):
    if geom is None:
        return []
    gtype = geom.GetGeometryType()
    if gtype in (ogr.wkbPolygon, ogr.wkbPolygon25D):
        return iter_linestrings(geom.Boundary())
    if gtype in (ogr.wkbLineString, ogr.wkbLineString25D):
        return [geom.Clone()]
    if gtype in (ogr.wkbMultiLineString, ogr.wkbMultiLineString25D,
                 ogr.wkbGeometryCollection, ogr.wkbGeometryCollection25D,
                 ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D):
        lines = []
        for i in range(geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            lines.extend(iter_linestrings(g))
        return lines
    return []


def linestring_features_from_geom(geom):
    features = []
    for line in iter_linestrings(geom):
        coords = []
        for i in range(line.GetPointCount()):
            x, y, _ = line.GetPoint(i)
            coords.append([x, y])
        if len(coords) < 2:
            continue
        features.append({
            "type": "Feature",
            "properties": {},
            "geometry": {"type": "LineString", "coordinates": coords}
        })
    return features
