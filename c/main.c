#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

struct index_to_species
{
  const char *index;
  const char *label;
};
typedef struct index_to_species index_to_species;
struct Atoms
{
  index_to_species *index_map;
  int nspecies;
  char *vectors;
  int counter;
};
typedef struct Atoms Atoms;

size_t getIndexOf(const char a[], size_t size, char value)
{
  size_t index = 2;
  while (index < size && a[index] != value)
    ++index;
  return (index == size ? -1 : index);
}

int main(int argc, char *argv[])
{
  if (sizeof(*argv) / sizeof(char *) < 1 || argc < 1)
  {
    printf("Error: no arguments supplied, exiting...\n");
    return EXIT_FAILURE;
  }
  // CONQUEST coords_next.dat file has the format
  // 3x3 array of numbers: lattice vectors
  // total number of atoms in unit cell
  // vectors of atom coords, then species index and T/F T/F T/F to determine whether atom can move in x,y,z respectively
  // the species index is used alongside command line arguments, or  optionally a file, to find the elements
  // T/F flags are removed
  // There are less than 256 elements so int is fine.
  // Expect program arguments to be supplied as ./C2V coord_file 1 <element1> 2 <element2>..., i.e. number of arguments is even
  // argv[0] = program name
  // argv[1] = path to CONQUEST coord file
  char **savedArgv = argv;
  int savedArgc = argc;
  if ((savedArgc - 2) % 2 != 0)
  {
    printf("Error: number of arguments supplied must be even.");
    return EXIT_FAILURE;
  }
  // Check every even index, that is argv[2], argv[4] etc. is a positive integer
  int bit = (1 << (sizeof(int) * 8 - 1));
  for (int i = 2; i < savedArgc; i += 2)
  {
    printf("%d", i);
    int msb = atoi(savedArgv[i]) & bit;
    if (msb)
    {
      printf("Error: species index at argument %d must be a positive integer.", i);
      return EXIT_FAILURE;
    }
  }
  // scanf("%s", coord_file);
  FILE *coord_file_pointer = fopen(savedArgv[1], "rb");
  // file_exists(input_file_pointer, input_string);
  if (coord_file_pointer == NULL)
  {
    printf("Error: unable to open file %s.", savedArgv[1]);
    return EXIT_FAILURE;
  }

  // char coord_file[100];
  char coord_output_path[1024];
  printf("Please enter absolute path to destination of output file.\n");
  int user_coord_path = scanf("%s", coord_output_path);
  printf("user input path %s\n", coord_output_path);
  FILE *output_file_pointer = fopen(coord_output_path, "wb");
  // file_exists(input_file_pointer, input_string);
  if (output_file_pointer == NULL)
  {
    printf("Error: unable to create output file %s.", coord_output_path);
    return EXIT_FAILURE;
  }
  fprintf(output_file_pointer, "\n");

  char *buffer = (char *)malloc(sizeof(char) * 1024);
  long numbytes;
  // fseek(coord_file_pointer, 0L, SEEK_END);
  // /*buffer to hold the input text */
  // buffer = (char *)calloc(numbytes, sizeof(char));

  // /* memory error */
  if (buffer == NULL)
    return 1;
  int natoms;
  int yes_natoms = fscanf(coord_file_pointer, "%d", &natoms);
  // int *nspecies = (int *)calloc((savedArgc - 2) / 2, sizeof(int));
  // if (nspecies == NULL)
  // {
  //   return EXIT_FAILURE;
  // }
  int total_species = sizeof((savedArgc - 2) / 2) / sizeof(int);
  printf("total species %d", total_species);
  // double sorted_atom_vectors[natoms][3];
  char *x, *y, *z;
  char species_index;
  // a struct thats keeps track of the passed in index to the passed in species label
  index_to_species(*index_map) = (index_to_species *)calloc(total_species, sizeof(index_to_species));
  for (int i = 0; i <= total_species; i++)
  {
    int j = (i == total_species ? i - 1 : i);
    index_map[j].index = (char *)malloc(sizeof(char) * 16);
    index_map[j].label = (char *)malloc(sizeof(char) * 16);
    printf("argv 2*i+1: %s\n", savedArgv[2 * i + 3]);
    printf("argv %s\n", savedArgv[2 * i + 2]);
    index_map[j].index = savedArgv[2 * i + 2];
    index_map[j].label = savedArgv[2 * i + 3];
  }
  // the T T T flags don't exist in other file formats, so we just keep them here but don't use them outside of temp storage
  char move_x, move_y, move_z;
  // /* copy all the text into the buffer */
  // fread(buffer, sizeof(char), numbytes, coord_file_pointer);
  // numbytes = ftell(coord_file_pointer);
  // argv[3], argv[5] etc are strings for elements
  // this is mean to make a string of the form "A Mn O" esntc. This is used to identify species
  // User may pass in something like 1 A 3 X 2 B 4 C where A,X,B,C are elements and the integers
  // are their indices in CONQUEST_input
  // So we require the map of species index to element label, and then we write them to VASP file

  char *ele_string = (char *)calloc(32, sizeof(char));
  // size_t arg_strlen = sizeof("Bi") / sizeof(char);
  // char *final_strng = (char *)malloc(sizeof(char) * (2 * arg_strlen) * (argc));
  if (ele_string == NULL)
  {
    printf("Error allocating string buffers.");
    return EXIT_FAILURE;
  }
  // strncpy(ele_string, savedArgv[3], sizeof(savedArgv[3]) / sizeof(char));
  // strncat(ele_string, " ", sizeof(" "));
  for (int i = 2; i <= argc - 2; i += 2)
  {
    // snprintf(final_strng, arg_strlen, "%s ", savedArgv[i + 1]);
    strcat(ele_string, savedArgv[i + 1]);
    strcat(ele_string, " ");
  }
  printf("ele_strng %s\n", ele_string);
  // return 0;
  fprintf(output_file_pointer, "%s\n", ele_string);
  // printf("ele_strng %s\n", ele_string);

  fprintf(output_file_pointer, "1.0\n");
  // first 3 lines of CONQUEST file always the lattice coords
  printf("??????");
  char *coord_buff = (char *)calloc(1024, sizeof(char) * sizeof(double));
  if (coord_buff == NULL)
  {
    printf("Error allocating buffer to store coordinates.");
    return EXIT_FAILURE;
  }
  while (fgets(coord_buff, sizeof(coord_buff), coord_file_pointer) != NULL)
  {
    fprintf(output_file_pointer, "  %s\n", coord_buff);
  }
  // for (int i = 0; i < 3; i++)
  // {
  // char *xd = fgets(coord_buff, sizeof(coord_buff), coord_file_pointer);
  // printf("xd %s\n", coord_buff);
  // if (fscanf(coord_file_pointer, "%s %s %s", x, y, z) != 3)
  // if (fgets(coord_buff, sizeof(coord_buff), coord_file_pointer))
  // {
  //   fprintf(output_file_pointer, "  %s\n", coord_buff);
  // }
  // }
  //   free(buffer);
  // fclose(output_file_pointer);
  // fclose(coord_file_pointer);
  // printf("Error writing CONQUEST lattice vectors to destination.\n");
  // return EXIT_FAILURE;

  Atoms All_Atoms_Vectors[total_species];
  for (int i = 0; i < total_species; i++)
  {
    All_Atoms_Vectors[i].index_map = &index_map[i];
    // All_Atoms_Vectors[i].nspecies = nspecies[i];
    All_Atoms_Vectors[i].nspecies = 0;
  }
  int atom_index = 1;
  char tempstr[10];
  char tempcoord[512];
  // size_t what = fread(tempcoord, sizeof(tempcoord), sizeof(tempcoord) * sizeof(char), coord_file_pointer);
  // int what = getline(&tempcoord, sizeof(&tempcoord), coord_file_pointer);
  char *what = fgets(tempcoord, sizeof(tempcoord), coord_file_pointer);
  // printf("what %d", what);
  // for (int n = 0; n < 1; ++n)
  // {
  printf("%s", what);
  // }
  return 0;
  while (EOF != fscanf(coord_file_pointer, "%s %s %s %c %c %c %c", x, y, z, &species_index, &move_x, &move_y, &move_z))
  {
    // Store the totals of each species, necessary for VASP
    // size_t index_in_argv = getIndexOf(*argv, sizeof(int), species_index);
    for (int i = 0; i < total_species; i++)
    {
      if (species_index == snprintf(tempstr, sizeof(tempstr), "%s", (All_Atoms_Vectors[i].index_map)->index))
      {
        All_Atoms_Vectors[i].nspecies++;
      }
    }
    // nspecies[index_in_argv - 2]++;
  }
  for (int i = 0; i < total_species; i++)
  {
    All_Atoms_Vectors[i].vectors = (char *)malloc(All_Atoms_Vectors[i].nspecies * 3 * sizeof(char) * 256); // initialise memory to store vector: [x1, y1, z1, x2, y2, z2, ...]
  }
  // check number of atoms in each species = natoms
  int total_num_species = 0;
  for (int i = 0; i < total_species; i++)
  {
    total_num_species += All_Atoms_Vectors[i].nspecies;
  }
  if (total_num_species != natoms)
  {
    printf("Error: total number of atoms summed from species not equal to initial total number of atoms.\n");
    return EXIT_FAILURE;
  }
  // Reset file pointer back to beginning
  fseek(coord_file_pointer, 0, SEEK_SET);
  // Skip first 4 lines, which are the lattice vectors and the total number of atoms
  // We care only about extracting the atom coordinates now
  char *temp = fgets(buffer, sizeof(buffer), coord_file_pointer);
  temp = fgets(buffer, sizeof(buffer), coord_file_pointer);
  temp = fgets(buffer, sizeof(buffer), coord_file_pointer);
  temp = fgets(buffer, sizeof(buffer), coord_file_pointer);
  while (EOF != fscanf(coord_file_pointer, "%s %s %s %c %c %c %c", x, y, z, &species_index, &move_x, &move_y, &move_z))
  {
    for (int i = 0; i < total_species; i++)
    {
      if (strncmp(All_Atoms_Vectors[i].index_map->index, &species_index, sizeof(char) * sizeof(species_index) * sizeof(All_Atoms_Vectors[i].index_map->index)))
      {
        All_Atoms_Vectors[i].vectors[All_Atoms_Vectors[i].counter] = *x;
        All_Atoms_Vectors[i].vectors[All_Atoms_Vectors[i].counter + 1] = *y;
        All_Atoms_Vectors[i].vectors[All_Atoms_Vectors[i].counter + 2] = *z;
        All_Atoms_Vectors[i].counter += 3;
      }
      continue;
    }
    continue;
    // fprintf(output_file_pointer, "  %lf %lf %lf\n", &x, &y, &z);
  }

  fprintf(output_file_pointer, "%s", ele_string);
  // Write number of each species here
  fprintf(output_file_pointer, "\n");
  fprintf(output_file_pointer, "Direct");
  // We now have all the structs of each atom and their atom coordinates
  // We wrote the species to output in order they were passed in to the program
  for (int i = 0; i < total_species; i++)
  {
    for (int j = 0; j < sizeof(All_Atoms_Vectors[i].vectors) / sizeof(char); j += 3)
    {
      fprintf(output_file_pointer, "  %s, %s, %s\n", &All_Atoms_Vectors[i].vectors[j], &All_Atoms_Vectors[i].vectors[j + 1], &All_Atoms_Vectors[i].vectors[j + 2]);
    }
  }
  // if (fgets(buffer, sizeof(buffer), coord_file_pointer) == NULL)
  // {
  //   // printf("The input file did not have the correct number of values and/or type of values. The order is (one value on one line): long N_x, N_y, N_a; double length_x, length_y, t_f, lambda, t_D");

  //   fclose(coord_file_pointer);
  //   free(buffer);
  //   fclose(output_file_pointer);
  //   return EXIT_FAILURE;
  // }
  free(buffer);
  free(ele_string);
  free(coord_buff);
  fclose(output_file_pointer);
  fclose(coord_file_pointer);
  return 0;
}