#include <iostream>
#include "../include/lz77.h"

int main(int argc, char *argv[])
{

  LZ77::Encoder* lz77_encoder = new LZ77::Encoder();
  char* file_name = argv[1];
  std::string out_file = "out.lz77";

  lz77_encoder->FillBuffer(file_name);
  lz77_encoder->Encode();

  return 0;  
}