#include <graphlab.hpp>

int main(int argc, char** argv)
{
  //Graphlab initialization
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  //Main body
  dc.cout() << "Hello World!\n";
  std::cout << "Output per core! (From every core)" << std::endl;

  //Main body ends
  graphlab::mpi_tools::finalize();
  return 0;
}
