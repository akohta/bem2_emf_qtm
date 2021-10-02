#include "bem2_emf_qtm.h"

int main(int argc,char **argv)
{
  DOMD md;
  
  read_data(argc,argv,&md);
  print_data(&md);
  //print_data_MKSA(&md);
  initialize_domd(&md);
  output_node_particles(argv[5],&md);

  solve_bieq(&md);
  dat_write(argv[5],&md);

  mfree_domd(&md);
  return 0;
}
