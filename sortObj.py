import numpy as np
import sys


def parse_obj(obj_filename):
    """Parse the OBJ file and return the vertices, normals, textures, and faces."""
    vertices = []
    normals = []
    textures = []
    faces = []
    other_lines = []  # This will store lines like `mtllib`, `usemtl`, comments, etc.
    
    with open(obj_filename, 'r') as file:
        for line in file:
            if line.startswith('v '):  # Vertex line
                parts = line.strip().split()
                vertex = [float(parts[1]), float(parts[2]), float(parts[3])]
                vertices.append(vertex)
            elif line.startswith('vn '):  # Vertex normal line
                parts = line.strip().split()
                normal = [float(parts[1]), float(parts[2]), float(parts[3])]
                normals.append(normal)
            elif line.startswith('vt '):  # Texture coordinate line
                parts = line.strip().split()
                texture = [float(parts[1]), float(parts[2])]
                textures.append(texture)
            elif line.startswith('f '):  # Face line
                parts = line.strip().split()
                face = [p.split('/') for p in parts[1:]]
                faces.append(face)
            else:
                other_lines.append(line)  # Other data (comments, materials, etc.)

    return np.array(vertices), np.array(normals), np.array(textures), faces, other_lines

def calculate_face_centroid(vertices, face):
    """Calculate the centroid of a triangular face."""
    vertex_indices = [int(f[0]) - 1 for f in face]  # Use the vertex positions only
    return np.mean(vertices[vertex_indices], axis=0)

def float_to_morton(value, max_bits=10):
    """Convert a float in the range [0, 1] to a Morton code part."""
    morton_code = 0
    value_scaled = int(value * ((1 << max_bits) - 1))  # Scale to the number of bits
    
    for i in range(max_bits):
        morton_code |= ((value_scaled >> i) & 1) << (3 * i)  # 3D interleaving for Morton code
    return morton_code

def calculate_morton_code(centroid, bbox_min, bbox_max, max_bits=10):
    """Convert the centroid position to a Morton code."""
    # Normalize the centroid to [0, 1] range
    bbox_d = (bbox_max - bbox_min)
    
    # Put 1 on any 0
    bbox_d[bbox_d == 0] = 1 
    normalized_centroid = (centroid - bbox_min) / bbox_d 

    # Convert each coordinate to Morton code part
    x_morton = float_to_morton(normalized_centroid[0], max_bits)
    y_morton = float_to_morton(normalized_centroid[1], max_bits)
    z_morton = float_to_morton(normalized_centroid[2], max_bits)
    
    # Interleave the Morton codes for x, y, and z
    morton_code = 0
    for i in range(max_bits):
        morton_code |= ((x_morton >> (3 * i)) & 1) << (3 * i)
        morton_code |= ((y_morton >> (3 * i)) & 1) << (3 * i + 1)
        morton_code |= ((z_morton >> (3 * i)) & 1) << (3 * i + 2)
    
    #print(f"Centroid: {normalized_centroid}, Morton code: {x_morton}, {y_morton}, {z_morton}, total: {morton_code}")

    return morton_code

def reorder_faces_by_morton_code(vertices, faces, max_bits=10):
    """Reorder faces by Morton code proximity."""
    # Calculate bounding box of all vertices
    bbox_min = np.min(vertices, axis=0)
    bbox_max = np.max(vertices, axis=0)
    
    # Compute Morton codes for centroids of all faces
    morton_codes_and_faces = []
    for face in faces:
        centroid = calculate_face_centroid(vertices, face)
        morton_code = calculate_morton_code(centroid, bbox_min, bbox_max, max_bits)
        morton_codes_and_faces.append((morton_code, face))
    
    # Sort faces by the Morton code
    morton_codes_and_faces.sort(key=lambda x: x[0])

    #print([morton_code for morton_code, _ in morton_codes_and_faces])
    
    # Return the reordered faces
    return [face for _, face in morton_codes_and_faces]

def write_obj(output_filename, vertices, normals, textures, faces, other_lines):
    """Write the vertices, normals, textures, and reordered faces into a new OBJ file."""
    with open(output_filename, 'w') as file:
        # Write other OBJ data like comments, material, etc.
        for line in other_lines:
            file.write(line)
        
        # Write vertices
        for vertex in vertices:
            file.write(f'v {vertex[0]} {vertex[1]} {vertex[2]}\n')
        
        # Write texture coordinates
        for texture in textures:
            file.write(f'vt {texture[0]} {texture[1]}\n')
        
        # Write vertex normals
        for normal in normals:
            file.write(f'vn {normal[0]} {normal[1]} {normal[2]}\n')
        
        # Write faces
        for face in faces:
            face_str = ' '.join(['/'.join(map(str, [f[i] if f[i] != '' else '' for i in range(3)])) for f in face])
            file.write(f'f {face_str}\n')

def reorder_obj_faces_by_morton_code(input_obj, output_obj, max_bits=10):
    # Parse the OBJ file
    vertices, normals, textures, faces, other_lines = parse_obj(input_obj)
    
    # Reorder faces by Morton code
    reordered_faces = reorder_faces_by_morton_code(vertices, faces, max_bits)
    
    # Write the reordered faces into a new OBJ file
    write_obj(output_obj, vertices, normals, textures, reordered_faces, other_lines)

# Example usage

if len(sys.argv) != 3:
    print("Usage: python sortObj.py input.obj output.obj")
    sys.exit(1)

input_obj_file = sys.argv[1]
output_obj_file = sys.argv[2]
reorder_obj_faces_by_morton_code(input_obj_file, output_obj_file)

print(f"Faces reordered by Morton codes and written to {output_obj_file}")
