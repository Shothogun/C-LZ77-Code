#include <iostream>
#include "../include/lz77.h"
#include "../include/huffman.h"

int main(int argc, char *argv[])
{

  LZ77::Encoder *lz77_encoder = new LZ77::Encoder();
  LZ77::Decoder *lz77_decoder = new LZ77::Decoder();

  if(argc < 2)
  {
    throw std::invalid_argument("Less than 2 arguments");
  }

  char *file_name = argv[1];
  char *out_file = argv[2];
  std::string compressed_file = out_file;
  compressed_file += ".lz77";

  std::string decompressed_file = out_file;
  decompressed_file += ".decompressed";

  lz77_encoder->FillBuffer(file_name);
  lz77_encoder->Encode();
  lz77_encoder->CompressToFile(compressed_file);

  lz77_decoder->DecompressFromFile(compressed_file);
  lz77_decoder->Decode("offset");
  lz77_decoder->Decode("length");
  lz77_decoder->DecompressLZ77Code();
  lz77_decoder->DecompressToFile(decompressed_file);

  return 0;
}