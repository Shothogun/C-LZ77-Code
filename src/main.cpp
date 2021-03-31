#include <iostream>
#include "../include/lz77.h"
#include "../include/huffman.h"

int main(int argc, char *argv[])
{

  LZ77::Encoder* lz77_encoder = new LZ77::Encoder();
  LZ77::Decoder* lz77_decoder = new LZ77::Decoder();

  char* file_name = argv[1];
  std::string out_file = "out.lz77";

  lz77_encoder->FillBuffer(file_name);
  lz77_encoder->Encode();
  lz77_encoder->CompressToFile(out_file);

  lz77_decoder->DecompressFromFile(out_file);
  lz77_decoder->Decode("offset");
  lz77_decoder->Decode("length");

  return 0;  
}