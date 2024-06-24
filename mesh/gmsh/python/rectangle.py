import gmsh
import scipy.io
import argparse

def makeMesh(args):
    """
    Generates a rectangular mesh using Gmsh and saves it to a .mat file.

    Parameters:
    args (Namespace): A namespace object containing the mesh parameters.
        - width (float): The width of the rectangle.
        - height (float): The height of the rectangle.
        - nElementsX (int): The number of elements along the x-axis.
        - nElementsZ (int): The number of elements along the z-axis.
        - filename (str): The name of the output .mat file.

    Saves:
    A .mat file containing the mesh data with the following variables:
        - nodes: A list of node coordinates.
        - elements: A list of element connectivity.
        - boundaries: A list of boundary nodes for each edge.
    """
    gmsh.initialize()
    gmsh.model.add("Rectangle")
    
    gmsh.model.geo.addPoint(0,                    0, 0, 1.0, 1)
    gmsh.model.geo.addPoint(0,           args.width, 0, 1.0, 2)
    gmsh.model.geo.addPoint(args.height, args.width, 0, 1.0, 3)  
    gmsh.model.geo.addPoint(args.height,          0, 0, 1.0, 4)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setRecombine(2, 1)

    gmsh.model.geo.mesh.setTransfiniteCurve(1, args.nElementsX, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, args.nElementsX, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, args.nElementsZ, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, args.nElementsZ, "Progression", 1)

    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    
    gmsh.model.mesh.generate(2)
    
    boundaries = []
    for line in [1, 2, 3, 4]:
        nodes, _, _ = gmsh.model.mesh.getNodes(1, line, True)
        boundaries.append((line, nodes))

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()

    nodes = [list(node_coords[i*3:(i+1)*3]) for i in range(len(node_tags))]

    elementIDs, nodeIds = gmsh.model.mesh.getElementsByType(10)
        
    nodeIds = nodeIds.reshape((len(elementIDs), 9))

    elements = [list(nodeids) for (eid, nodeids) in zip(elementIDs, nodeIds)]

    scipy.io.savemat(f'{args.filename}', {
            'nodes': nodes,
            'elements': elements,
            'boundaries': boundaries,
    })

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a rectangular mesh using Gmsh.')
    parser.add_argument('--width', type=float, required=True, default=100, help='Width of the rectangle')
    parser.add_argument('--height', type=float, required=True, default=20, help='Height of the rectangle')
    parser.add_argument('--nElementsX', type=int, required=True, default=20, help='Number of elements along the x-axis')
    parser.add_argument('--nElementsZ', type=int, required=True, default=4, help='Number of elements along the z-axis')
    parser.add_argument('--filename', type=str, required=True, help='Output filename for the .mat file')

    args = parser.parse_args()
    makeMesh(args)