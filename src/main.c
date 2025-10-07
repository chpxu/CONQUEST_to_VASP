#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

struct index_to_species
{
  int index;
  const char *label;
};

struct Atoms
{
  struct index_to_species *index_map;
  int nspecies;
  double *vectors;
  int counter;
};
size_t getIndexOf(const char a[], size_t size, int value)
{
  size_t index = 2;
  while (index < size && a[index] != value)
    ++index;
  return (index == size ? -1 : index);
}
int main(int argc, char *argv[])
{
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
  // if ((argc - 1) % 2 != 0)
  // {
  //   printf("Error: number of arguments supplied must be even.");
  //   return EXIT_FAILURE;
  // }
  // Check every even index, that is argv[2], argv[4] etc. is a positive integer
  int bit = (1 << (sizeof(int) * 8 - 1));
  for (int i = 1; i < floor(argc / 2); i = 2 * i)
  {
    int msb = atoi(argv[i]) & bit;
    if (msb)
    {
      printf("Error: species index at argument %d must be a positive integer.", i);
      return EXIT_FAILURE;
    }
  }
  // scanf("%s", coord_file);
  FILE *coord_file_pointer = fopen(argv[1], "r");
  // file_exists(input_file_pointer, input_string);
  if (coord_file_pointer == NULL)
  {
    fclose(coord_file_pointer);
    printf("Error: unable to open file %s.", argv[1]);
    return EXIT_FAILURE;
  }

  // char coord_file[100];
  char coord_output_path[1024];
  printf("Please enter absolute path to destination of output file.\n");
  int user_coord_path = scanf("%s", coord_output_path);

  FILE *output_file_pointer = fopen(coord_output_path, "rw");
  // file_exists(input_file_pointer, input_string);
  if (output_file_pointer == NULL)
  {
    fclose(output_file_pointer);
    printf("Error: unable to create output file %s.", coord_output_path);
    return EXIT_FAILURE;
  }
  char *buffer[1024];
  long numbytes;
  // fseek(coord_file_pointer, 0L, SEEK_END);
  // /*buffer to hold the input text */
  // buffer = (char *)calloc(numbytes, sizeof(char));

  // /* memory error */
  if (buffer == NULL)
    return 1;
  int natoms;
  fscanf(coord_file_pointer, "%d", &natoms);
  int *nspecies = (int *)calloc((argc - 2) / 2, sizeof(int));
  int total_species = sizeof(nspecies) / sizeof(int);
  double sorted_atom_vectors[natoms][3];
  double x, y, z;
  int species_index;
  // a struct thats keeps track of the passed in index to the passed in species label
  struct index_to_species(*index_map)[total_species];
  for (int i = 1; i <= total_species; i++)
  {
    index_map[i]->index = (int)argv[2 * i];
    index_map[i]->label = argv[2 * i + 1];
  }
  // the T T T flags don't exist in other file formats, so we just keep them here but don't use them outside of temp storage
  char move_x, move_y, move_z;
  // /* copy all the text into the buffer */
  // fread(buffer, sizeof(char), numbytes, coord_file_pointer);
  // numbytes = ftell(coord_file_pointer);
  // argv[3], argv[5] etc are strings for elements
  // this is mean to make a string of the form "A Mn O" etc. This is used to identify species
  // User may pass in something like 1 A 3 X 2 B 4 C where A,X,B,C are elements and the integers
  // are their indices in CONQUEST_input
  // So we require the map of species index to element label, and then we write them to VASP file

  char ele_string[6 * argc * sizeof(char)];
  strncpy(ele_string, argv[3], sizeof(argv[3]));
  strncat(ele_string, " ", sizeof(" "));
  for (int i = 2; i < floor(argc / 2); i = 2 * i + 1)
  {
    size_t strlen = sizeof(argv[i]);
    char final_strng[10] = "";
    snprintf(final_strng, strlen, "%s ", argv[i]);
    strncat(ele_string, final_strng, strlen);
  }

  fprintf(output_file_pointer, "%s", ele_string);
  fprintf(output_file_pointer, "\n");
  fprintf(output_file_pointer, "1.0\n");
  // first 3 lines of CONQUEST file always the lattice coords
  for (int i = 0; i < 3; i++)
  {
    double x, y, z;
    if (fscanf(coord_file_pointer, "%lf %lf %lf", &x, &y, &z) != 3)
    {
      free(buffer);
      fclose(output_file_pointer);
      fclose(coord_file_pointer);
      printf("Error writing CONQUEST lattice vectors to destination.\n");
      return EXIT_FAILURE;
    }
    else
    {
      fprintf(output_file_pointer, "  %lf %lf %lf\n", x, y, z);
    }
  }

  int atom_index = 1;
  while (EOF != fscanf(coord_file_pointer, "%lf %lf %lf %d %c %c %c", &x, &y, &z, &species_index, &move_x, &move_y, &move_z))
  {
    // Store the totals of each species, necessary for VASP
    size_t index_in_argv = getIndexOf(*argv, sizeof(int), species_index);
    nspecies[index_in_argv - 2]++;
  }
  struct Atoms All_Atoms_Vectors[total_species];
  for (int i = 0; i < total_species; i++)
  {
    All_Atoms_Vectors[i].index_map = index_map[i];
    All_Atoms_Vectors[i].nspecies = nspecies[i];
    All_Atoms_Vectors[i].nspecies = 0;
    All_Atoms_Vectors[i].vectors = (double *)calloc(All_Atoms_Vectors[i].nspecies * 3, sizeof(double)); // initialise memory to store vector: [x1, y1, z1, x2, y2, z2, ...]
  }
  // Reset file pointer back to beginning
  fseek(coord_file_pointer, 0, SEEK_SET);
  // Skip first 4 lines, which are the lattice vectors and the total number of atoms
  // We care only about extracting the atom coordinates now
  fgets(*buffer, sizeof(buffer), coord_file_pointer);
  fgets(*buffer, sizeof(buffer), coord_file_pointer);
  fgets(*buffer, sizeof(buffer), coord_file_pointer);
  fgets(*buffer, sizeof(buffer), coord_file_pointer);

  while (EOF != fscanf(coord_file_pointer, "%lf %lf %lf %d %c %c %c", &x, &y, &z, &species_index, &move_x, &move_y, &move_z))
  {
    for (int i = 0; i < total_species; i++)
    {
      if (All_Atoms_Vectors[i].index_map->index == species_index)
      {
        All_Atoms_Vectors[i].vectors[All_Atoms_Vectors[i].counter] = x;
        All_Atoms_Vectors[i].vectors[All_Atoms_Vectors[i].counter + 1] = y;
        All_Atoms_Vectors[i].vectors[All_Atoms_Vectors[i].counter + 2] = z;
        All_Atoms_Vectors[i].counter += 3;
      }
      continue;
    }
    continue;
    // fprintf(output_file_pointer, "  %lf %lf %lf\n", &x, &y, &z);
  }

  fprintf(output_file_pointer, "\n");
  // Write number of each species here
  fprintf(output_file_pointer, "\n");
  fprintf(output_file_pointer, "Direct");
  // We now have all the structs of each atom and their atom coordinates
  // We wrote the species to output in order they were passed in to the program
  for (int i = 0; i < total_species; i++)
  {
    for (int j = 0; j < sizeof(All_Atoms_Vectors[i].vectors) / sizeof(double); j += 3)
    {
      fprintf(output_file_pointer, "  %lf, %lf, %lf\n", All_Atoms_Vectors[i].vectors[j], All_Atoms_Vectors[i].vectors[j + 1], All_Atoms_Vectors[i].vectors[j + 2]);
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
  fclose(output_file_pointer);
  fclose(coord_file_pointer);
  return 0;
}