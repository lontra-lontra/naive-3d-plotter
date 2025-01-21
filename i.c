#include <stdio.h>
#include <stdlib.h>
#include <time.h>




#define img_max_x 1080
#define img_max_y 1920

#pragma pack(1) // Ensures there is no padding between struct members

unsigned int rand_seed = 1;

void srand(unsigned int seed) {
    rand_seed = seed;
}

int rand(void) {
    rand_seed = rand_seed * 1103515245 + 12345;
    return (rand_seed / 65536) % 32768;
}
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


typedef struct {
    coordinates_3d coordinates_3d;
    RGB color;
    int line_type;
    int size;
} cube;

void save_cube_list_to_txt(const char *filename, cube *cubes, int count) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    for (int i = 0; i < count; i++) {
        fprintf(file, "cube %d: Coordinates=(%.2f, %.2f, %.2f), Color=(%d, %d, %d), LineType=%d, Size=%d\n",
                i,
                cubes[i].coordinates_3d.x, cubes[i].coordinates_3d.y, cubes[i].coordinates_3d.z,
                cubes[i].color.r, cubes[i].color.g, cubes[i].color.b,
                cubes[i].line_type, cubes[i].size);
    }

    fclose(file);
}



void retrieve_cube_list_from_txt(const char *filename, cube **cubes, int *count) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        *count = 0;
        return;
    }

    // Count the number of cubes in the file
    char ch;
    *count = 0;
    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            (*count)++;
        }
    }
    rewind(file);

    // Allocate memory for the cubes
    *cubes = malloc(*count * sizeof(cube));
    if (*cubes == NULL) {
        perror("Failed to allocate memory");
        fclose(file);
        *count = 0;
        return;
    }

    // Read the cubes from the file
    for (int i = 0; i < *count; i++) {
        fscanf(file, "cube %*d: Coordinates=(%f, %f, %f), Color=(%hhu, %hhu, %hhu), LineType=%d, Size=%d\n",
               &(*cubes)[i].coordinates_3d.x, &(*cubes)[i].coordinates_3d.y, &(*cubes)[i].coordinates_3d.z,
               &(*cubes)[i].color.r, &(*cubes)[i].color.g, &(*cubes)[i].color.b,
               &(*cubes)[i].line_type, &(*cubes)[i].size);
    }

    fclose(file);
}



void retrieve_cube_list_from_txt__(const char *filename, cube *cubes, int count) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    for (int i = 0; i < count; i++) {
        fscanf(file, "cube %*d: Coordinates=(%f, %f, %f), Color=(%hhu, %hhu, %hhu), LineType=%d, Size=%d\n",
               &cubes[i].coordinates_3d.x, &cubes[i].coordinates_3d.y, &cubes[i].coordinates_3d.z,
               &cubes[i].color.r, &cubes[i].color.g, &cubes[i].color.b,
               &cubes[i].line_type, &cubes[i].size);
    }

    fclose(file);
}


void save_small_matrix_to_txt(const char *filename, RGB ***small_matrix, int width, int height, int depth) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
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

void retrieve_small_matrix_from_txt(const char *filename, RGB ***small_matrix, int height, int width, int depth) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    int i, j, k, r, g, b;
    while (fscanf(file, "(%d, %d, %d): R=%d, G=%d, B=%d\n", &i, &j, &k, &r, &g, &b) == 6) {
        if (i >= 0 && i < height && j >= 0 && j < width && k >= 0 && k < depth) {
            small_matrix[i][j][k].r = (unsigned char)r;
            small_matrix[i][j][k].g = (unsigned char)g;
            small_matrix[i][j][k].b = (unsigned char)b;
        }
    }

    fclose(file);
}





coordinates_2d project(coordinates_3d point3d) {
    coordinates_2d point2d;
    
    // Cabinet projection formulas
    point2d.x = point3d.x + point3d.z/1.41;
    point2d.y = point3d.y + point3d.z/1.41;
    
    return point2d;
}


int write_in_direction(coordinates_3d c3d, coordinates_3d direction, float distance, RGB **matrix, int dotted) {
    coordinates_2d c2d;
    img_coord img;
    int wrote = 0;
    int resolution = 999;
    for (int i = 0; i < distance * resolution; i++) {
        c3d.x = c3d.x + direction.x / resolution;
        c3d.y = c3d.y + direction.y / resolution;
        c3d.z = c3d.z + direction.z / resolution;
        c2d = project(c3d);
        img = conv(c2d);
        if (img.x != -1) {
            if (dotted == 1) {
                if (i % (resolution/3) < 111 ) { // Adjust the modulus value to change the dot pattern
                    matrix[img.x][img.y].r = 255;
                    matrix[img.x][img.y].g = 255;
                    matrix[img.x][img.y].b = 255;
                }
                else{
                    //printf("didnt write\n");
                }
            } else {
                matrix[img.x][img.y].r = 255;
                matrix[img.x][img.y].g = 255;
                matrix[img.x][img.y].b = 255;
            }
            wrote = 1;
        }
    }
    return wrote;
}



int write_cube(RGB **matrix, coordinates_3d reference_vertice, float size, int line_type) { 
    int wrote = 0;
    // reference vertice is the vertice of the cube that is closest to the origin
    reference_vertice.x = reference_vertice.x * size;
    reference_vertice.y = reference_vertice.y * size;
    reference_vertice.z = reference_vertice.z * size;

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
        wrote += write_in_direction(start, direction, 1.0, matrix, line_type)*(1-wrote);
    }
    return wrote;
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
        write_in_direction(start, direction, 1.0, matrix,1);
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




int write_face(RGB **matrix, coordinates_3d vertices[4]) {
    int wrote = 0;
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
        wrote += write_in_direction(start, direction, 1.0, matrix,0)*(1-wrote);
    }
    return wrote;
}

int cube_with_faces(RGB **matrix, coordinates_3d reference_vertice, float size){
    int wrote = 0;
    coordinates_3d vertices[4] = {
        {0, 0, 0}, {size, 0, 0}, {size, size, 0}, {0, size, 0}
    };
    for (int i = 0; i < 4; i++) {
        vertices[i].x += reference_vertice.x*size;
        vertices[i].y += reference_vertice.y*size;
        vertices[i].z += reference_vertice.z*size;
    }
    //printf("vertices[0]: %f, %f, %f\n", vertices[0].x, vertices[0].y, vertices[0].z);
    wrote += write_face(matrix, vertices)*(1-wrote);

    vertices[0] = (coordinates_3d){0, size, 0};
    vertices[1] = (coordinates_3d){size, size, 0};
    vertices[2] = (coordinates_3d){size, size, size};
    vertices[3] = (coordinates_3d){0, size, size};
    for (int i = 0; i < 4; i++) {
        vertices[i].x += reference_vertice.x*size;
        vertices[i].y += reference_vertice.y*size;
        vertices[i].z += reference_vertice.z*size;
    }

    wrote += write_face(matrix, vertices)*(1-wrote);
    vertices[0] = (coordinates_3d){size, size, 0};
    vertices[1] = (coordinates_3d){size, size, size};
    vertices[2] = (coordinates_3d){size, 0, size};
    vertices[3] = (coordinates_3d){size, 0, 0};
    for (int i = 0; i < 4; i++) {
        vertices[i].x += reference_vertice.x*size;
        vertices[i].y += reference_vertice.y*size;
        vertices[i].z += reference_vertice.z*size;
    }
    wrote += write_face(matrix, vertices)*(1-wrote);


}



void add_cube(const char *filename, float x, float y, float z) {
    // Open the file in append mode
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    // Create a new cube with the given coordinates and default values
    cube new_cube;
    new_cube.coordinates_3d.x = x;
    new_cube.coordinates_3d.y = y;
    new_cube.coordinates_3d.z = z;
    new_cube.color.r = 255;
    new_cube.color.g = 255;
    new_cube.color.b = 255;
    new_cube.line_type = 1;
    new_cube.size = 12;

    // Write the new cube to the file
    fprintf(file, "cube: Coordinates=(%.2f, %.2f, %.2f), Color=(%d, %d, %d), LineType=%d, Size=%d\n",
            new_cube.coordinates_3d.x, new_cube.coordinates_3d.y, new_cube.coordinates_3d.z,
            new_cube.color.r, new_cube.color.g, new_cube.color.b,
            new_cube.line_type, new_cube.size);

    // Close the file
    fclose(file);
}


void write_image_from_cube_list_file(const char *cube_list_filename, const char *output_image_filename) {
    int width = img_max_y;
    int height = img_max_x;

    // Dynamically allocate memory for the image matrix
    RGB **matrix = malloc(height * sizeof(RGB *));
    for (int i = 0; i < height; i++) {
        matrix[i] = malloc(width * sizeof(RGB));
    }

    // Initialize the matrix with black color
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            matrix[i][j].r = 0;
            matrix[i][j].g = 0;
            matrix[i][j].b = 0;
        }
    }

    // Retrieve the cube list from the text file
    int cube_count = 0;
    cube *cubes = NULL;
    FILE *file = fopen(cube_list_filename, "r");
    if (file != NULL) {
        // Count the number of cubes in the file
        char ch;
        while ((ch = fgetc(file)) != EOF) {
            if (ch == '\n') {
                cube_count++;
            }
        }
        rewind(file);

        // Allocate memory for the cubes
        //cubes = malloc(cube_count * sizeof(cube));

        // Read the cubes from the file
        retrieve_cube_list_from_txt(cube_list_filename, &cubes, &cube_count);

        fclose(file);
    } else {
        perror("Failed to open cube list file");
        return;
    }
    // Sort the cubes by x + y - z
    for (int i = 0; i < cube_count - 1; i++) {
        for (int j = 0; j < cube_count - i - 1; j++) {
            float value1 = cubes[j].coordinates_3d.x + cubes[j].coordinates_3d.y - cubes[j].coordinates_3d.z;
            float value2 = cubes[j + 1].coordinates_3d.x + cubes[j + 1].coordinates_3d.y - cubes[j + 1].coordinates_3d.z;
            if (value1 > value2) {
                cube temp = cubes[j];
                cubes[j] = cubes[j + 1];
                cubes[j + 1] = temp;
            }
        }
    }
    
    // Draw the cubes on the matrix
    for (int i = 0; i < cube_count; i++) {
        if (cubes[i].line_type == 0)
        {
            if (cube_with_faces(matrix, cubes[i].coordinates_3d, cubes[i].size) != 1) {
                //printf("Cube not drawn at (%.2f, %.2f, %.2f)\n", cubes[i].coordinates_3d.x, cubes[i].coordinates_3d.y, cubes[i].coordinates_3d.z);
            } else {
                //printf("Cube drawn at (%.2f, %.2f, %.2f)\n", cubes[i].coordinates_3d.x, cubes[i].coordinates_3d.y, cubes[i].coordinates_3d.z);
            }
        }
        else
        {
            if (write_cube(matrix, cubes[i].coordinates_3d, cubes[i].size,1) != 1) {
                //printf("Cube not drawn at (%.2f, %.2f, %.2f)\n", cubes[i].coordinates_3d.x, cubes[i].coordinates_3d.y, cubes[i].coordinates_3d.z);
            } else {
                //printf("Cube drawn at (%.2f, %.2f, %.2f)\n", cubes[i].coordinates_3d.x, cubes[i].coordinates_3d.y, cubes[i].coordinates_3d.z);
            }
        }
    }

    // Write the matrix to a BMP file
    write_bmp(output_image_filename, matrix, width, height);

    // Free the dynamically allocated memory
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(cubes);
}

void remove_duplicate_cubes_and_sort(const char *filename) {
    int cube_count = 0;
    cube *cubes = NULL;
    FILE *file = fopen(filename, "r");
    if (file != NULL) {
        // Count the number of cubes in the file
        char ch;
        while ((ch = fgetc(file)) != EOF) {
            if (ch == '\n') {
                cube_count++;
            }
        }
        rewind(file);

        // Allocate memory for the cubes
        //cubes = malloc(cube_count * sizeof(cube));

        // Read the cubes from the file
        retrieve_cube_list_from_txt(filename, &cubes, &cube_count);

        fclose(file);
    } else {
        perror("Failed to open file");
        return;
    }

    // Remove duplicates
    int new_count = 0;
    for (int i = 0; i < cube_count; i++) {
        int duplicate = 0;
        for (int j = 0; j < new_count; j++) {
            if (cubes[i].coordinates_3d.x == cubes[j].coordinates_3d.x &&
                cubes[i].coordinates_3d.y == cubes[j].coordinates_3d.y &&
                cubes[i].coordinates_3d.z == cubes[j].coordinates_3d.z) {
                duplicate = 1;
                break;
            }
        }
        if (!duplicate) {
            cubes[new_count++] = cubes[i];
        }
    }

    // Sort the cubes by x + y - z
    for (int i = 0; i < new_count - 1; i++) {
        for (int j = 0; j < new_count - i - 1; j++) {
            float value1 = cubes[j].coordinates_3d.x + cubes[j].coordinates_3d.y - cubes[j].coordinates_3d.z;
            float value2 = cubes[j + 1].coordinates_3d.x + cubes[j + 1].coordinates_3d.y - cubes[j + 1].coordinates_3d.z;
            if (value1 > value2) {
                cube temp = cubes[j];
                cubes[j] = cubes[j + 1];
                cubes[j + 1] = temp;
            }
        }
    }

    // Save the updated cube list to the text file
    save_cube_list_to_txt(filename, cubes, new_count);

    // Free the dynamically allocated memory
    free(cubes);
}

void add_cube_and_save_img(const char *filename, float x, float y, float z, int t) {

    remove_duplicate_cubes_and_sort("cubes.txt");
    // Parse the coordinates from the command line arguments

    // Open the cube list file
    const char *cube_list_filename = "cubes.txt";
    const char *output_image_filename = "output_image.bmp";

    // Retrieve the cube list from the text file
    int cube_count = 0;
    cube *cubes = NULL;
    FILE *file = fopen(cube_list_filename, "r");
    if (file != NULL) {
        // Count the number of cubes in the file
        char ch;
        while ((ch = fgetc(file)) != EOF) {
            if (ch == '\n') {
                cube_count++;
            }
        }
        rewind(file);

        // Allocate memory for the cubes
        //cubes = malloc(cube_count * sizeof(cube));

        // Read the cubes from the file
        //retrieve_cube_list_from_txt(cube_list_filename, cubes, cube_count);
        retrieve_cube_list_from_txt(cube_list_filename, &cubes, &cube_count);
        

        fclose(file);
    } else {
        perror("Failed to open cube list file");
    }

    // Check if the cube already exists
    int cube_exists = 0;
    for (int i = 0; i < cube_count; i++) {
        if (cubes[i].coordinates_3d.x == x && cubes[i].coordinates_3d.y == y && cubes[i].coordinates_3d.z == z) {
            cube_exists = 1;
            break;
        }
    }

    // If the cube does not exist, add it to the list
    if (!cube_exists) {
        cubes = realloc(cubes, (cube_count + 1) * sizeof(cube));
        cubes[cube_count].coordinates_3d.x = x;
        cubes[cube_count].coordinates_3d.y = y;
        cubes[cube_count].coordinates_3d.z = z;
        cubes[cube_count].color.r = 255;
        cubes[cube_count].color.g = 255;
        cubes[cube_count].color.b = 255;
        cubes[cube_count].line_type = t;
        cubes[cube_count].size = 12;
        cube_count++;

        // Save the updated cube list to the text file
        save_cube_list_to_txt(cube_list_filename, cubes, cube_count);
    }

    // Create the image from the cube list
    write_image_from_cube_list_file(cube_list_filename, output_image_filename);

    // Free the dynamically allocated memory
    free(cubes);
}





int main(int argc, char *argv[])
{
    float x, y, z;
    int t;
    if (argc == 4) {
        x = atof(argv[1]);
        y = atof(argv[2]);
        z = atof(argv[3]);
        t = 's';
    }
    else 
    {
        if (argc == 5){
            x = atof(argv[1]);
            y = atof(argv[2]);
            z = atof(argv[3]);
            t = atof(argv[4]);
        }
        else{
            fprintf(stderr, "Usage: %s <x> <y> <z>\n", argv[0]);
            return 1;
        }
    }
    add_cube_and_save_img("cubes.txt", x, y, z, t);
    int steps = 20; 
    float step_size = 1;

        coordinates_3d current_position = {0, 0, 0};
        srand(time(NULL));

        // Keep track of visited positions
        int visited_count = 0;
        coordinates_3d *visited_positions = malloc(steps * sizeof(coordinates_3d));
        visited_positions[visited_count++] = current_position;

        for (int i = 0; i < steps; i++) {
            printf("Step %d: Moving to (%.2f, %.2f, %.2f)\n", i + 1, current_position.x, current_position.y, current_position.z);
            float dx = 0, dy = 0, dz = 0;
            int direction = rand() % 3;
            if (direction == 0) {
            dx = step_size; // Always move in positive x direction
            } else if (direction == 1) {
            dy = step_size; // Always move in positive y direction
            } else {
            dz = (rand() % 2 ? 1 : -1) * step_size; // Move in either positive or negative z direction
            }

            for (int repeat = 0; repeat < 2; repeat++) {
                coordinates_3d new_position = {
                    current_position.x + dx,
                    current_position.y + dy,
                    current_position.z + dz
                };

                // Check if the new position has already been visited
                int already_visited = 0;
                for (int j = 0; j < visited_count; j++) {
                    if (visited_positions[j].x == new_position.x &&
                        visited_positions[j].y == new_position.y &&
                        visited_positions[j].z == new_position.z) {
                        already_visited = 1;
                        break;
                    }
                }

                if (!already_visited) {
                    current_position = new_position;
                    visited_positions[visited_count++] = current_position;
                    add_cube("cubes.txt", current_position.x, current_position.y, current_position.z);
                }
            }
        }
        // Create the image from the cube list
        printf("foi");
        add_cube_and_save_img("cubes.txt", x, y, z, t);
        free(visited_positions);
    return 0;
}
