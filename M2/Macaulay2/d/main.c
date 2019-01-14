#include <stdio.h>

extern int Macaulay2_main(int argc, char **argv);
extern char timestamp[];

int main(int argc, char **argv)
{
  printf("*** RZ's test code! timestamp is: %s\n\n", timestamp);
  return Macaulay2_main(argc, argv);
}
