import json
import vtk

# base file name
basename = 'DEM_cuenca_Laura'

#open topology file
with open(basename+'.mesh') as f:
    mesh = json.load(f)

#open parameter file
with open(basename+'.param') as f:
    param = json.load(f)

#init the vtk datastructure
vtu = vtk.vtkUnstructuredGrid()
vtuwriter = vtk.vtkXMLUnstructuredGridWriter()
vtuwriter.SetFileName(basename+'.vtu')
vtuwriter.SetInputData(vtu)

vtu_points = vtk.vtkPoints()
vtu_triangles = vtk.vtkCellArray()

vtu_points.SetNumberOfPoints(len(mesh['mesh']['vertex']))

# elevation isn't stored in the param or mesh file, so compute from the z like mesher does
vtu_cells  = {'Elevation': vtk.vtkFloatArray()}
vtu_cells['Elevation'].SetName('elevation')

# loop over the mesh elements and reconstruct the triangle
for e in mesh['mesh']['elem']:
    v0 = e[0]
    v1 = e[1]
    v2 = e[2]

    vtu_points.SetPoint(v0, mesh['mesh']['vertex'][v0][0], mesh['mesh']['vertex'][v0][1], mesh['mesh']['vertex'][v0][2])
    vtu_points.SetPoint(v1, mesh['mesh']['vertex'][v1][0], mesh['mesh']['vertex'][v1][1], mesh['mesh']['vertex'][v1][2])
    vtu_points.SetPoint(v2, mesh['mesh']['vertex'][v2][0], mesh['mesh']['vertex'][v2][1], mesh['mesh']['vertex'][v2][2])

    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, v0)
    triangle.GetPointIds().SetId(1, v1)
    triangle.GetPointIds().SetId(2, v2)

    vtu_triangles.InsertNextCell(triangle)

    vtu_cells['Elevation'].InsertNextTuple1((mesh['mesh']['vertex'][v0][2] +
                                                             mesh['mesh']['vertex'][v1][2] +
                                                             mesh['mesh']['vertex'][v2][2]) / 3.)

# insert all the parameters from the param file
for key, data in param.items():
    k =  key

    #sanity check to ensure we don't add a 0 length
    if len(param[k]) == 0 :
        continue
    vtu_cells[k]=vtk.vtkFloatArray()
    vtu_cells[k].SetName(k)

    for e in param[k]:
        vtu_cells[k].InsertNextTuple1(e)

# write the proj4 wkt output
proj4 = vtk.vtkStringArray()
proj4.SetNumberOfComponents(1);
proj4.SetName("proj4");
proj4.InsertNextValue(mesh['mesh']['proj4']);

vtu.GetFieldData().AddArray(proj4);

vtu.SetPoints(vtu_points)
vtu.SetCells(vtk.VTK_TRIANGLE,vtu_triangles)
for p in vtu_cells.values():
    vtu.GetCellData().AddArray(p)
vtuwriter.Write()


