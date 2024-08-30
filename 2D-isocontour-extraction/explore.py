import numpy as np
import vtk

def read_vti_file(filename):
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def write_vtp_file(polydata, filename):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()

def extract_isocontour(grid, isovalue):
    dimensions = grid.GetDimensions()
    spacing = grid.GetSpacing()
    origin = grid.GetOrigin()
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    point_indices = {}
    visited_points = set()

    # Iterate over each cell
    for j in range(dimensions[1] - 1):
        for i in range(dimensions[0] - 1):
            cell_values = [
                grid.GetPointData().GetScalars().GetTuple1((i + di) + (j + dj) * dimensions[0]) 
                for di, dj in [(0, 0), (1, 0), (1, 1), (0, 1)]
            ]
            edge_intersections = []

            # Check for intersections on each edge
            for k in range(4):
                l = (k + 1) % 4
                if (cell_values[k] < isovalue and cell_values[l] >= isovalue) or (cell_values[k] >= isovalue and cell_values[l] < isovalue):
                    t = (isovalue - cell_values[k]) / (cell_values[l] - cell_values[k])
                    x = origin[0] + (i + k // 2 + t) * spacing[0]
                    y = origin[1] + (j + k % 2) * spacing[1]
                    edge_intersections.append((x, y))

            # Connect the edge intersections to form a contour
            if len(edge_intersections) >= 2:
                for idx in range(len(edge_intersections) - 1):
                    point1 = edge_intersections[idx]
                    point2 = edge_intersections[idx + 1]
                    key1 = (round(point1[0], 6), round(point1[1], 6))
                    key2 = (round(point2[0], 6), round(point2[1], 6))

                    # Insert or reuse points
                    if key1 not in visited_points:
                        vertex_idx1 = points.InsertNextPoint(point1[0], point1[1], 0)
                        visited_points.add(key1)
                    else:
                        vertex_idx1 = point_indices[key1]
                    
                    if key2 not in visited_points:
                        vertex_idx2 = points.InsertNextPoint(point2[0], point2[1], 0)
                        visited_points.add(key2)
                    else:
                        vertex_idx2 = point_indices[key2]

                    # Create line segment
                    line = vtk.vtkLine()
                    line.GetPointIds().SetId(0, vertex_idx1)
                    line.GetPointIds().SetId(1, vertex_idx2)
                    polys.InsertNextCell(line)

                    # Update point_indices
                    point_indices[key1] = vertex_idx1
                    point_indices[key2] = vertex_idx2

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(polys)

    return polydata



def main():
    # Read VTI file
    filename = "isabel_2D.vti"
    grid = read_vti_file(filename)

    # User input for isovalue
    isovalue = float(input("Enter isovalue: "))

    # Extract isocontour
    isocontour = extract_isocontour(grid, isovalue)

    # Write VTP file
    output_filename = "output_file.vtp"
    write_vtp_file(isocontour, output_filename)

    print("Isocontour extracted and saved as", output_filename)

if __name__ == "__main__":
    main()
