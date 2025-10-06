#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
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
  scanf("%s", coord_output_path);

  FILE *output_file_pointer = fopen(coord_output_path, "rw");
  // file_exists(input_file_pointer, input_string);
  if (output_file_pointer == NULL)
  {
    fclose(output_file_pointer);
    printf("Error: unable to create output file %s.", output_file_pointer);
    return EXIT_FAILURE;
  }
  char *buffer[1024];
  long numbytes;
  // fseek(coord_file_pointer, 0L, SEEK_END);
  // /*buffer to hold the input text */
  // buffer = (char *)calloc(numbytes, sizeof(char));

  // /* memory error */
  // if (buffer == NULL)
  //   return 1;

  // /* copy all the text into the buffer */
  // fread(buffer, sizeof(char), numbytes, coord_file_pointer);
  // numbytes = ftell(coord_file_pointer);
  // argv[3], argv[5] etc are strings for elements
  // this is mean to make a string of the form "A Mn O" etc. This is used to identify species
  char ele_string[6 * argc];
  strncpy(ele_string, argv[3], sizeof(argv[3]));
  strncat(ele_string, " ", sizeof(" "));
  for (int i = 2; i < floor(argc / 2); i = 2 * i + 1)
  {
    size_t strlen = sizeof(argv[i]);
    char final_strng[10] = "";
    snprintf(final_strng, strlen, "%s ", argv[i]);
    strncat(ele_string, final_strng, strlen);
  }
  fprintf(output_file_pointer, ele_string);
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
  int natoms;
  fscanf(coord_file_pointer, "%d", &natoms);
  int nspecies[(argc - 2) / 2] = {0};
  double x, y, z;
  int species_index;
  char move_x, move_y, move_z;
  while (EOF != fscanf(coord_file_pointer, "%lf %lf %lf %d %c %c %c", &x, &y, &z, &species_index, &move_x, &move_y, &move_z))
  {
    fprintf(output_file_pointer, "  %lf %lf %lf\n", &x, &y, &z);
  }
  fprintf(output_file_pointer, ele_string);
  fprintf(output_file_pointer, "\n");
  // Write number of each species here
  fprintf(output_file_pointer, "\n");
  fprintf(output_file_pointer, "Direct");

  if (fgets(buffer, sizeof(buffer), coord_file_pointer) == NULL)
  {
    // printf("The input file did not have the correct number of values and/or type of values. The order is (one value on one line): long N_x, N_y, N_a; double length_x, length_y, t_f, lambda, t_D");

    fclose(coord_file_pointer);
    free(buffer);
    return EXIT_FAILURE;
  }
  free(buffer);
  fclose(output_file_pointer);
  fclose(coord_file_pointer);
  return 0;
}