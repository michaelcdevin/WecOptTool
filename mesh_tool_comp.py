import pygalmesh
import gmsh
import pygmsh
import capytaine as cpt
import numpy as np

r1 = 0.88
r2 = 0.35
h1 = 0.17
h2 = 0.37
freeboard = 0.01

def pygmsh_mesh(mesh_size_factor=0.1):
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
                                r1)
        cone = geom.add_cone([0, 0, -h1],
                                [0, 0, -h2],
                                r1, r2)
        geom.translate(cyl, [0, 0, freeboard])
        geom.translate(cone, [0, 0, freeboard])
        geom.boolean_union([cyl, cone])
        mesh = geom.generate_mesh()

    return mesh

def pygalmesh_mesh(max_cell_circumradius= 0.1):
    """Generate surface mesh of hull.

    Parameters
    ----------
    max_cell_circumradius
        Upper bound of the circumradius allowed in each mesh element.
        Lower values result in a finer mesh.
    """

    cyl = pygalmesh.Cylinder(0, h1, r1, 0.1)
    cone = pygalmesh.Cone(r1, -h1, 0.1)
    cyl = pygalmesh.Translate(cyl, [0, 0, freeboard])
    cone = pygalmesh.Translate(cone, [0, 0, freeboard])
    combined_shape = pygalmesh.Union([cyl, cone])
    print('before generate_mesh')
    mesh = pygalmesh.generate_mesh(
        combined_shape,
        max_cell_circumradius=max_cell_circumradius,
        max_edge_size_at_feature_edges=0.15,
        min_facet_angle=25,
        max_radius_surface_delaunay_ball=0.15,
        max_circumradius_edge_ratio=2.0,
        )
    print('finish generating mesh')

    # mesh = pygalmesh.generate_mesh(
    #     pygalmesh.Ball([0.0, 0.0, 0.0], 1.0),
    #     min_facet_angle=30.0,
    #     max_radius_surface_delaunay_ball=0.1,
    #     max_facet_distance=0.025,
    #     max_circumradius_edge_ratio=2.0,
    #     max_cell_circumradius=lambda x: abs(np.sqrt(np.dot(x, x)) - 0.5) / 5 + 0.025,
    # )

    return mesh

mesh1 = pygmsh_mesh()
mesh2 = pygalmesh_mesh()
fb1 = cpt.FloatingBody.from_meshio(mesh1)
fb2 = cpt.FloatingBody.from_meshio(mesh2)
fb1.show()
fb2.show()