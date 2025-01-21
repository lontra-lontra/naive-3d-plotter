#include <stdio.h>
#include <stdlib.h>
#define img_max_x 1080
#define img_max_y 1920

#pragma pack(1) // Ensures there is no padding between struct members
typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;

typedef struct {
	float x;
	float y;
	float z;
} coordinates_3d;

typedef struct {
	float x;
	float y;
} coordinates_2d;


typedef struct {
	int x;
	int y;
} img_coord;

img_coord conv(coordinates_2d c_2d){
    img_coord img;
    img.x = img_max_x - (int)c_2d.y;
    img.y = (int)c_2d.x;
    if (img.x >= img_max_x || img.y >= img_max_y) {
        return (img_coord){-1, -1};
    }
    if (img.x < 0 || img.y < 0) {
        return (img_coord){-1, -1};
    }

    return img;
}









coordinates_2d project(coordinates_3d point3d) {
    coordinates_2d point2d;
    
    // Cabinet projection formulas
    point2d.x = point3d.x + point3d.z/1.41;
    point2d.y = point3d.y + point3d.z/1.41;
    
    return point2d;
}


void write_in_direction(coordinates_3d c3d, coordinates_3d direction,float distance, RGB **matrix) {
    coordinates_2d c2d;
    img_coord img;
    float resolution = 1000;
    for (int i = 0; i < distance*resolution; i++) 
    {              
        c3d.x = c3d.x + direction.x/resolution;
        c3d.y = c3d.y + direction.y/resolution;
        c3d.z = c3d.z + direction.z/resolution;
        c2d = project(c3d);
        img = conv(c2d);
        if (img.x != -1) {
            matrix[img.x][img.y].r = 255;
            matrix[img.x][img.y].g = 255;
            matrix[img.x][img.y].b = 255;
        }
    }
}



void write_cube(RGB **matrix, coordinates_3d reference_vertice, float size){ // reference vertice is the vertice of the cube that is closest to the origin

coordinates_3d vertices[8] = {
    {0, 0, 0}, {size, 0, 0}, {size, size, 0}, {0, size, 0},
    {0, 0, size}, {size, 0, size}, {size, size, size}, {0, size, size}
    };

    int edges[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Bottom face
    {4, 5}, {5, 6}, {6, 7}, {7, 4}, // Top face
    {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Vertical edges
    };

    for (int i = 0; i < 12; i++) {
        coordinates_3d start = vertices[edges[i][0]];
        start.x += reference_vertice.x*size;
        start.y += reference_vertice.y*size;
        start.z += reference_vertice.z*size; 
        coordinates_3d end = vertices[edges[i][1]];
        end.x += reference_vertice.x*size;
        end.y += reference_vertice.y*size;
        end.z += reference_vertice.z*size;   
        coordinates_3d direction = {
            end.x - start.x,
            end.y - start.y,
            end.z - start.z
            };
        write_in_direction(start, direction, 1.0, matrix);
    }
}


void write_cube_with_actual_cooridinates(RGB **matrix, coordinates_3d reference_vertice, float size){ // reference vertice is the vertice of the cube that is closest to the origin

coordinates_3d vertices[8] = {
    {0, 0, 0}, {size, 0, 0}, {size, size, 0}, {0, size, 0},
    {0, 0, size}, {size, 0, size}, {size, size, size}, {0, size, size}
    };

    int edges[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Bottom face
    {4, 5}, {5, 6}, {6, 7}, {7, 4}, // Top face
    {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Vertical edges
    };

    for (int i = 0; i < 12; i++) {
        coordinates_3d start = vertices[edges[i][0]];
        start.x += reference_vertice.x;
        start.y += reference_vertice.y;
        start.z += reference_vertice.z; 
        coordinates_3d end = vertices[edges[i][1]];
        end.x += reference_vertice.x;
        end.y += reference_vertice.y;
        end.z += reference_vertice.z;   
        coordinates_3d direction = {
            end.x - start.x,
            end.y - start.y,
            end.z - start.z
            };
        write_in_direction(start, direction, 1.0, matrix);
    }
}





void write_bmp(const char *filename, RGB **matrix, int width, int height) {
    FILE *f;
    unsigned char bmpfileheader[14] = {'B', 'M'};
    unsigned char bmpinfoheader[40] = {0};
    int filesize = 54 + 3 * width * height;  // File size

    bmpfileheader[ 2] = (unsigned char)(filesize);
    bmpfileheader[ 3] = (unsigned char)(filesize >> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize >> 16);
    bmpfileheader[ 5] = (unsigned char)(filesize >> 24);

    bmpfileheader[10] = 54;  // Data offset

    bmpinfoheader[ 0] = 40;  // Header size
    bmpinfoheader[ 4] = (unsigned char)(width);
    bmpinfoheader[ 5] = (unsigned char)(width >> 8);
    bmpinfoheader[ 6] = (unsigned char)(width >> 16);
    bmpinfoheader[ 7] = (unsigned char)(width >> 24);
    bmpinfoheader[ 8] = (unsigned char)(height);
    bmpinfoheader[ 9] = (unsigned char)(height >> 8);
    bmpinfoheader[10] = (unsigned char)(height >> 16);
    bmpinfoheader[11] = (unsigned char)(height >> 24);
    bmpinfoheader[12] = 1;   // Number of color planes
    bmpinfoheader[14] = 24;  // Bits per pixel

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);

    // Writing the RGB matrix to the BMP file (bottom-up)
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j < width; j++) {
            RGB pixel = matrix[i][j];
            fwrite(&pixel, 3, 1, f);
        }
    }
    fclose(f);
}




void write_face(RGB **matrix, coordinates_3d vertices[4]) {
    int edges[4][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };
    // Paint the interior black
    coordinates_3d min = vertices[0];
    coordinates_3d max = vertices[0];
    for (int i = 1; i < 4; i++) {
        if (vertices[i].x < min.x) min.x = vertices[i].x;
        if (vertices[i].y < min.y) min.y = vertices[i].y;
        if (vertices[i].z < min.z) min.z = vertices[i].z;
        if (vertices[i].x > max.x) max.x = vertices[i].x;
        if (vertices[i].y > max.y) max.y = vertices[i].y;
        if (vertices[i].z > max.z) max.z = vertices[i].z;
    }

    for (float x = min.x; x <= max.x; x += 0.1) {
        for (float y = min.y; y <= max.y; y += 0.1) {
            for (float z = min.z; z <= max.z; z += 0.1) {
                coordinates_3d point = {x, y, z};
                coordinates_2d projected = project(point);
                img_coord img = conv(projected);
                if (img.x != -1) {
                    matrix[img.x][img.y].r = 0;
                    matrix[img.x][img.y].g = 0;
                    matrix[img.x][img.y].b = 0;
                }
            }
        }
    }




    for (int i = 0; i < 4; i++) {
        coordinates_3d start = vertices[edges[i][0]];
        coordinates_3d end = vertices[edges[i][1]];
        coordinates_3d direction = {
            end.x - start.x,
            end.y - start.y,
            end.z - start.z
        };
        write_in_direction(start, direction, 1.0, matrix);
    }

}

void cube_with_faces(RGB **matrix, coordinates_3d reference_vertice, float size){
    coordinates_3d vertices[4] = {
        {0, 0, 0}, {size, 0, 0}, {size, size, 0}, {0, size, 0}
    };
    for (int i = 0; i < 4; i++) {
        vertices[i].x += reference_vertice.x*size;
        vertices[i].y += reference_vertice.y*size;
        vertices[i].z += reference_vertice.z*size;
    }
    printf("vertices[0]: %f, %f, %f\n", vertices[0].x, vertices[0].y, vertices[0].z);
    write_face(matrix, vertices);     
    vertices[0] = (coordinates_3d){0, size, 0};
    vertices[1] = (coordinates_3d){size, size, 0};
    vertices[2] = (coordinates_3d){size, size, size};
    vertices[3] = (coordinates_3d){0, size, size};
    for (int i = 0; i < 4; i++) {
        vertices[i].x += reference_vertice.x*size;
        vertices[i].y += reference_vertice.y*size;
        vertices[i].z += reference_vertice.z*size;
    }

    write_face(matrix, vertices); 
    vertices[0] = (coordinates_3d){size, size, 0};
    vertices[1] = (coordinates_3d){size, size, size};
    vertices[2] = (coordinates_3d){size, 0, size};
    vertices[3] = (coordinates_3d){size, 0, 0};
    for (int i = 0; i < 4; i++) {
        vertices[i].x += reference_vertice.x*size;
        vertices[i].y += reference_vertice.y*size;
        vertices[i].z += reference_vertice.z*size;
    }
    write_face(matrix, vertices); 


}



int main() {
    int width = img_max_y;
    int height = img_max_x;

    // Dynamically allocate memory for the image matrix
    RGB **matrix = malloc(height * sizeof(RGB *));
    for (int i = 0; i < height; i++) {
        matrix[i] = malloc(width * sizeof(RGB));
    }

    // Fill the matrix with gradient colors
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            matrix[i][j].r = (unsigned char)(0);   // Red gradient
            matrix[i][j].g = (unsigned char)(0);  // Green gradient
            matrix[i][j].b = (unsigned char)(0); // Blue gradient
        }
    }

    coordinates_3d face_vertices[4] = {
        {0, 0, 0}, {0, 0, 100}, {0, 50, 100}, {0, 50, 0}
    };

    // Create a 10x10 matrix
    int small_y = 160;
    int small_x = 160;

    // Define the dimensions for the small matrix
    int small_depth = 200;

    // Dynamically allocate memory for the small matrix
    RGB ***small_matrix = malloc(small_x * sizeof(RGB **));
    for (int i = 0; i < small_x; i++) {
        small_matrix[i] = malloc(small_y * sizeof(RGB *));
        for (int j = 0; j < small_y; j++) {
            small_matrix[i][j] = malloc(small_depth * sizeof(RGB));
        }
    }

    // Fill the small matrix with black color
    for (int i = 0; i < small_x; i++) {
        for (int j = 0; j < small_y; j++) {
            for (int k = 0; k < small_depth; k++) {
                small_matrix[i][j][k].r = 0;
                small_matrix[i][j][k].g = 0;
                small_matrix[i][j][k].b = 0;
            }
        }
    }
    for (int x = 0; x < 25; x++) 
    {
        for (int y = 0; y < 2; y++) 
        {
            for (int z = 0; z < 200; z++) 
            {
                //small_matrix[25+y][25+0][25+x].r = 255;
                small_matrix[100+x][y][z].r = 255;
                //small_matrix[25+x][25+y][100-16].r = 255;
            }
        }
    }
for (int sum = -small_depth; sum < small_x + small_y; sum++) {
    for (int i = 0; i < small_x; i++) {
        for (int j = 0; j < small_y; j++) {
            int k = i + j - sum;  // Calculate k directly from sum, i, and j

            // Check if k is within the valid range
            if (k >= 0 && k < small_depth) {
                if (small_matrix[i][j][k].r != 0) {
                    cube_with_faces(matrix, (coordinates_3d){i-1, j-1, k-1}, 12);
                    printf("Cube drawn at (%d, %d, %d)\n", i, j, k);
                }
            }
        }
    }
    printf("Sum: %d\n", sum);
}

    // Draw a cube with a reference vertex at (0, 0, 0) and size 10
    write_bmp("output_image.bmp", matrix, width, height);

    // Free the dynamically allocated memory
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
    }
    free(matrix);

    printf("Image written to 'output_image.bmp'\n");
        // Free the dynamically allocated memory for the small matrix
    for (int i = 0; i < small_x; i++) {
        for (int j = 0; j < small_y; j++) {
            free(small_matrix[i][j]);
        }
        free(small_matrix[i]);
    }
    free(small_matrix);
    return 0;
}

void save_small_matrix_to_txt(const char *filename, RGB ***small_matrix, int width, int height, int depth) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < depth; k++) {
                fprintf(file, "(%d, %d, %d): R=%d, G=%d, B=%d\n", 
                        i, j, k, 
                        small_matrix[i][j][k].r, 
                        small_matrix[i][j][k].g, 
                        small_matrix[i][j][k].b);
            }
        }
    }

    fclose(file);
}
