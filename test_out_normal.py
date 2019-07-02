
from electric_charge_reconstructed.mesh_data import mesh_data

air_mesh, air_boundary=mesh_data.read()

from electric_charge_reconstructed.outward_normal.outward_normal import outward_normal
def surface_normal(mesh, boundary, marker):
    triangles, points, _, table_l2g, _ = markered_surface.markered_surface(mesh, boundary, marker)
    npoints = len(points)
    ntriangles = len(triangles)
    triangles.resize(ntriangles*3)
    points.resize(npoints*3)
    result = outward_normal(points, triangles)
    return result,points

from electric_charge_reconstructed.markered_surface import markered_surface
result,points =surface_normal(air_mesh, air_boundary, 3)
for i in range(402):
    sum=0
    sum += (points[i*3]/0.002-result[i*3])*(points[i*3]/0.002-result[i*3])
    sum += (points[i*3+1]/0.002-result[i*3+1])*(points[i*3+1]/0.002-result[i*3+1])
    sum += ((points[i*3+2]-0.005)/0.002-result[i*3+2])*((points[i*3+2]-0.005)/0.002-result[i*3+2])
    print(sum)
