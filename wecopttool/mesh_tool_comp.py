import os
import pyvista as pv
import gmsh
import pygmsh
import capytaine as cpt
import numpy as np
import meshio
import time

wb_r1 = 0.88
wb_r2 = 0.35
h1 = 0.17
h2 = 0.37

freeboard = 0.01

t1=1.5
t2=0.355
t3=7.25
ah_r1=1.085
ah_r2=0.405
ah_r3=0.355
ofst=0.1

def get_total_cone_height(r1, r2, h):
    """Calculates the total theoretical cone length needed to achieve a
    flat-topped cone (aka a conical frustum) of the given `r2` at the given `h`.
    Parameters
    ----------
    r1
        Radius of cone at its base
    r2
        Radius of cone desired at `h`
    h
        Desired height of cone
    """
    s = np.sqrt((r1-r2)**2 + h**2)
    theta = np.arccos(h/s)
    return h + r2 / np.tan(theta)


def to_meshio(pv_mesh):
    """Converts a `pyvista.PolyData` mesh object to a `meshio._Mesh` object."""
    timestr = time.strftime('%y%m%d-%H%M%S')
    filename = f'mesh_{timestr}.vtk'
    pv.save_meshio(filename, pv_mesh)
    meshio_mesh = meshio.read(filename)
    os.remove(filename)
    return meshio_mesh

def pygmsh_wb_mesh(mesh_size_factor=0.1):
    """Generate surface mesh of hull.
    Parameters
    ----------
    mesh_size_factor
        Control for the mesh size. Smaller values give a finer mesh.
    """
    with pygmsh.occ.Geometry() as geom:
        gmsh.option.setNumber('Mesh.MeshSizeFactor', mesh_size_factor)
        cyl = geom.add_cylinder([0, 0, 0],
                                [0, 0, -h1],
                                wb_r1)
        cone = geom.add_cone([0, 0, -h1],
                                [0, 0, -h2],
                                wb_r1, wb_r2)
        geom.translate(cyl, [0, 0, freeboard])
        geom.translate(cone, [0, 0, freeboard])
        geom.boolean_union([cyl, cone])
        mesh = geom.generate_mesh()

    return mesh


def pyvista_wb_mesh(resolution=250):
    """Generates surface mesh of the WaveBot hull.

    Parameters
    ----------
    resolution
        Number of points on the circular face of each geometry object.
        Higher number results in a finer mesh.
    """
    h3 = get_total_cone_height(ah_r1, ah_r2, h2)
    cyl = pv.Cylinder(center=[0, 0, freeboard-h1/2.],
                      direction=[0, 0, -1],
                      radius=wb_r1,
                      height=h1,
                      resolution=resolution).triangulate()
    cone = pv.Cone(center=[0, 0, freeboard-h1-h3/2.],
                   direction=[0, 0, -1],
                   radius=wb_r1,
                   height=h3,
                   resolution=resolution).triangulate()
    cylint = pv.Cylinder(center=[0, 0, freeboard-h1-h2/2.],
                         direction=[0, 0, -1],
                         radius=wb_r1,
                         height=h2,
                         resolution=resolution).triangulate()
    mesh = cyl.boolean_union(cone.boolean_intersection(cylint))
    wavebot = to_meshio(mesh)

    return wavebot

t1 = time.perf_counter()
mesh2 = pyvista_wb_mesh()
t2 = time.perf_counter()
mesh1 = pygmsh_wb_mesh()
t3 = time.perf_counter()
print(t2-t1)
print(t3-t2)
fb1 = cpt.FloatingBody.from_meshio(mesh1)
fb2 = cpt.FloatingBody.from_meshio(mesh2)
fb1.show()
fb2.show()

def pygmsh_ah_mesh(mesh_size_factor:float=0.25):
        """
        """
        with pygmsh.occ.Geometry() as geom:
            gmsh.option.setNumber('Mesh.MeshSizeFactor', mesh_size_factor)
            cyl1 = geom.add_cylinder([0,0,0],
                                     [0,0,-t1],
                                     ah_r1)
            cone = geom.add_cone([0,0,-t1],
                                 [0,0,-t2],
                                 ah_r1,
                                 ah_r2)
            cylout = geom.add_cylinder([0,0,-1*(t1+t2)],
                                     [0,0,-t3],
                                     ah_r2)
            cylin = geom.add_cylinder([0,0,-1*(t1+t2)],
                                     [0,0,-t3],
                                     ah_r3)
            cyl2 = geom.boolean_difference(cylout,cylin)[0]
            wecGeom = geom.boolean_union(entities=[cyl1, cone, cyl2],
                                         delete_first=True)[0]

            geom.translate(wecGeom, [0, 0, ofst])
            mesh = geom.generate_mesh()

        return mesh

def pyvista_ah_mesh(resolution=100):
    """Generates surface mesh of the AquaHarmonics hull.

    Parameters
    ----------
    resolution
        Number of points on the circular face of each geometry object.
        Higher number results in a finer mesh.

    Returns
    -------
    mesh
        Surface mesh of hull.
    """

    ah_h = get_total_cone_height(ah_r1, ah_r2, t2)
    cyl1 = pv.Cylinder(center=[0,0,-t1/2.],
                       direction=[0,0,-1],
                       height=t1,
                       radius=ah_r1,
                       resolution=resolution).triangulate()
    cone = pv.Cone(center=[0,0,-(t1+ah_h/2.)],
                   direction=[0,0,-1],
                   height=ah_h,
                   radius=ah_r1,
                   resolution=resolution).triangulate()
    cylout = pv.Cylinder(center=[0,0,-(t1+t2+t3/2.)],
                         direction=[0,0,-1],
                         height=t3,
                         radius=ah_r2,
                         resolution=resolution).triangulate()
    cylin = pv.Cylinder(center=[0,0,-(t1+t2+t3/2.)],
                        direction=[0, 0, -1],
                        height=t3+.1,
                        radius=ah_r3,
                        resolution=resolution).triangulate()
    mesh = cyl1.boolean_union(cone
              ).boolean_union(cylout
              ).boolean_difference(cylin
              ).triangulate()
    mesh.translate([0, 0, ofst], inplace=True)
    aquaharmonics_wec = to_meshio(mesh)

    return aquaharmonics_wec

t1 = time.perf_counter()
mesh4 = pyvista_ah_mesh()
t2 = time.perf_counter()
mesh3 = pygmsh_ah_mesh()
t3 = time.perf_counter()
print(t2-t1)
print(t3-t2)
fb3 = cpt.FloatingBody.from_meshio(mesh3)
fb4 = cpt.FloatingBody.from_meshio(mesh4)
fb3.show()
fb4.show()
