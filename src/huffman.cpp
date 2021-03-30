#include "../include/huffman.h"
#include "../include/bitstream.h"
#define DEBUG 0
#define DECODE_DEBUG 0

int Huffman::Encoder::HowManyCharacters()
{
  return this->character_counter_;
}

int Huffman::Encoder::CharactersQuantity()
{
  return this->character_counter_;
}

void Huffman::Encoder::CountCharacters(int n_characters)
{
  this->character_counter_ += n_characters;
}

void Huffman::Encoder::FillBuffer(std::string file_path)
{

  // Read file as binary data
  std::ifstream f(file_path, std::ios::binary | std::ios::in);
  if (!f)
  {
    std::cout << "File not found\n";
    exit(0);
  };
  char c;

  // Gets each byte until no more
  // is left
  while (f.get(c))
  {
    // Gets each bit from the byte
    for (int i = 7; i >= 0; i--)
    {
      this->file_content_.push_back(((c >> i) & 1));
    }

    // Each byte is a character
    this->CountCharacters(1);
  }
}

void Huffman::Encoder::FillBuffer(std::vector<int> buffer)
{
  // Gets each single string from buffer
  for (auto const &content : buffer)
  {
    // Gets each bit from the byte
    for (int i = 7; i >= 0; i--)
    {
      this->file_content_.push_back(((content >> i) & 1));
    }

    // Each byte is a character
    this->CountCharacters(1);
  }

  // Returns the string sequence
  // representing the bits from the file
  std::vector<bool> file_buffer = this->GetBuffer();

  // A byte bitstream as a string
  std::string byte_bitstream = "";

  // Counts symbols
  for (uint base = 0; base < file_buffer.size(); base += 8)
  {
    byte_bitstream.clear();

    // Reads a byte sequence, represents a character
    for (uint i = 0; i < 8; i++)
    {
      byte_bitstream = byte_bitstream + std::to_string(file_buffer[base + i]);
    }

    this->CountSymbol(byte_bitstream);
  }

  this->ComputeProbabilityTable();
}

std::vector<bool> Huffman::Encoder::GetBuffer()
{
  return this->file_content_;
}

void Huffman::Encoder::CountSymbol(std::string character)
{
  this->symbol_table_[character]++;
}

std::map<std::string, double> Huffman::Encoder::GetSymbolTable()
{
  return this->symbol_table_;
}

std::map<std::string, std::string> Huffman::Encoder::GetSymbolEncode()
{
  return this->symbol_encode_;
}

void Huffman::Encoder::ComputeProbabilityTable()
{
  for (auto &x : this->GetSymbolTable())
  {
    this->symbol_table_[x.first] /= this->HowManyCharacters();
  }
}

void Huffman::Encoder::ComputeHuffmanCode()
{
  // Vector of pair probability and symbol
  std::vector<std::pair<double, std::string>> symbol_vector;

  // Vector that stores the symbols pairs chosen
  // during the combination, following its choice order
  std::vector<std::pair<std::string, std::string>> combined_pairs;

  // Stores the father node of the tree generated,
  // determines its child's code appending 0 or 1
  std::vector<std::string> father;

  // Maps the symbol to its code
  std::map<std::string, std::string> code_map;

  // Counts the ith alpha generated symbol from combination
  // of primitive symbol
  int alpha_count = 0;

  // Fill symbol vector
  for (auto const &x : this->GetSymbolTable())
  {
    symbol_vector.push_back(std::make_pair(x.second, x.first));
  }

  // Combines symbols
  while (symbol_vector.size() > 2)
  {
    // Sort in ascending order
    std::sort(symbol_vector.begin(), symbol_vector.end());

    // Get the 2 least probability symbols
    double first_symbol_probability = symbol_vector.front().first;
    std::string first_symbol = symbol_vector.front().second;
    symbol_vector.erase(symbol_vector.begin());

    double second_symbol_probability = symbol_vector.front().first;
    std::string second_symbol = symbol_vector.front().second;
    symbol_vector.erase(symbol_vector.begin());

    combined_pairs.push_back(std::make_pair(first_symbol, second_symbol));

    double new_probability = first_symbol_probability;
    new_probability += second_symbol_probability;

    symbol_vector.push_back(
        std::make_pair(new_probability,
                       "alpha_" + std::to_string(alpha_count)));

    father.push_back("alpha_" + std::to_string(alpha_count));

    alpha_count++;
  }

  // The remaining symbols will initiate
  // the code pattern
  auto it = symbol_vector.begin();
  auto next = symbol_vector.begin() + 1;
  code_map[it->second] = "1";
  code_map[next->second] = "0";

  std::reverse(father.begin(), father.end());
  std::reverse(combined_pairs.begin(), combined_pairs.end());
  auto current_father = father.begin();

  // Produces code to each symbol
  for (auto p : combined_pairs)
  {
    std::string code = code_map[*current_father];
    code_map[p.first] = code + "1";
    code_map[p.second] = code + "0";
    current_father++;
  }

  for (auto p : code_map)
  {
    if (p.first[0] == '0' || p.first[0] == '1')
    {
      this->symbol_encode_[p.first] = p.second;
    }
  }

  if (DEBUG)
  {
    std::cout << "-----------------------------\n"
              << "------ Combined Pair --------\n"
              << "-----------------------------\n";

    for (auto it : combined_pairs)
    {
      std::cout << it.first
                << " : "
                << it.second
                << "\n";
    }

    std::cout << "\n";

    std::cout << "-----------------------------\n"
              << "---------- Father -----------\n"
              << "-----------------------------\n";

    for (auto it : father)
    {
      std::cout << it
                << "\n";
    }

    std::cout << "\n";

    std::cout << "-----------------------------\n"
              << "---------- Code Map ---------\n"
              << "-----------------------------\n";

    for (auto const &x : code_map)
    {
      std::cout << x.first
                << ':'
                << x.second
                << "\n";
    }

    std::cout << "\n";
  }
}

void Huffman::Encoder::Encode()
{
  std::string byte_bitstream = "";
  std::string encoded_symbol = "";
  std::string bit = "";

  for (uint base = 0; base < this->file_content_.size(); base += 8)
  {
    byte_bitstream.clear();
    encoded_symbol.clear();

    // Constructs the byte symbol
    for (uint i = 0; i < 8; i++)
    {
      byte_bitstream += std::to_string(this->file_content_[base + i]);
    }

    encoded_symbol = this->symbol_encode_[byte_bitstream];

    // Gets the encoded symbol, bit by bit,
    // and concatenates to the encoded bitstream
    for (uint i = 0; i < encoded_symbol.size(); i++)
    {
      bit = encoded_symbol[i];
      this->encoded_data_.push_back(stoi(bit));
    }
  }

  if (DEBUG)
  {

    if (this->encoded_data_.size() > 1000)
    {
      std::cout << "WARNING: A TOO LONG BITSTREAM TO PRINT";
    }

    else
    {
      std::cout << "-----------------------------\n"
                << "---- Original bitstream -----\n"
                << "-----------------------------\n";

      for (auto x : this->file_content_)
      {
        std::cout << to_string(x);
      }
      std::cout << "\n\n";

      std::cout << "-----------------------------\n"
                << "---- Encoded bitstream ------\n"
                << "-----------------------------\n";

      for (auto x : this->encoded_data_)
      {
        std::cout << to_string(x);
      }
      std::cout << "\n\n";
    }

    std::cout << "-----------------------------\n"
              << "------- Symbol Encode -------\n"
              << "-----------------------------\n";

    for (auto const &x : this->symbol_encode_)
    {
      std::cout << x.first
                << ':'
                << x.second
                << std::endl;
    }

    std::cout << "\n";

    std::cout << "-----------------------------\n"
              << "- Symbols number: "
              << this->symbol_encode_.size()
              << "\n"
              << "-----------------------------\n\n";
  }

  // Bits per symbol in new encoding
  double average_rate = 0;

  for (auto const &x : this->GetSymbolTable())
  {
    average_rate += (this->symbol_encode_[x.first].size()) * x.second;
  }

  this->average_rate_ = average_rate;

  std::cout
      << "Average size:\t"
      << average_rate
      << " bits/symbol\n";

  std::cout
      << "Difference:\t"
      << (this->average_rate_ - this->entropy_)
      << " bits/symbol\n";

  double compression_rate = 1;
  compression_rate -= (double)this->encoded_data_.size() /
                      (double)this->file_content_.size();

  compression_rate *= 100;

  std::cout
      << "Original file size:\t\t\t"
      << this->file_content_.size() / 8
      << " bytes\n";

  std::cout
      << "Compressed file size without overhead:\t"
      << this->encoded_data_.size() / 8
      << " bytes\n";

  std::cout
      << "Liquid Compression rate: "
      << compression_rate
      << "%\n";
}

std::vector<std::string> Huffman::Encoder::GetEncodedContent()
{
  std::vector<std::string> encoded_content_vector;
  std::string byte_bitstream = "";
  std::string encoded_symbol = "";
  std::string bit = "";

  for (uint base = 0; base < this->file_content_.size(); base += 8)
  {
    byte_bitstream.clear();
    encoded_symbol.clear();

    // Constructs the byte symbol
    for (uint i = 0; i < 8; i++)
    {
      byte_bitstream += std::to_string(this->file_content_[base + i]);
    }

    encoded_symbol = this->symbol_encode_[byte_bitstream];

    encoded_content_vector.push_back(encoded_symbol);
  }

  return encoded_content_vector;
}

void Huffman::Encoder::CompressToFile(std::string file_name)
{
  // 2 Bytes: symbols number
  // tuples: array of tuples [a size code]
  // a    -> 1 Byte:  symbol
  // size -> 1 Byte:  encode size
  // code -> n Bytes: symbol encode

  // Empty Bitstream object
  Bitstream bstream;
  std::vector<bool> bstream_vector;

  // Inserts symbols number as bits
  for (int i = 0; i < 16; i++)
  {
    bstream.writeBit((this->symbol_encode_.size() >> (15 - i)) & 1);
    bstream_vector.push_back((this->symbol_encode_.size() >> (15 - i)) & 1);
  }

  std::string bit;

  // Inserts array of tuples as bits
  for (auto p : this->symbol_encode_)
  {
    // Inserts Symbol
    for (int i = 0; i < 8; i++)
    {
      bit = p.first[i];
      bstream.writeBit(stoi(bit));
      bstream_vector.push_back(stoi(bit));
    }

    // Inserts encode size
    for (uint i = 0; i < 8; i++)
    {
      bstream.writeBit((p.second.size() >> (7 - i)) & 1);
      bstream_vector.push_back((p.second.size() >> (7 - i)) & 1);
    }

    // Inserts symbol encode
    for (uint i = 0; i < p.second.size(); i++)
    {
      bit = p.second[i];
      bstream.writeBit(stoi(bit));
      bstream_vector.push_back(stoi(bit));
    }
  }

  // Inserts the encoded data
  for (uint i = 0; i < this->encoded_data_.size(); i++)
  {
    bstream.writeBit(this->encoded_data_[i]);
    bstream_vector.push_back(this->encoded_data_[i]);
  }

  bstream.flushesToFile(file_name);

  if (DEBUG)
  {
    std::cout << "-----------------------------\n"
              << "------- Encoded data --------\n"
              << "-----------------------------\n";

    for (auto x : bstream_vector)
    {
      std::cout << to_string(x);
    }
    std::cout << "\n\n";
  }

  std::cout
      << "Compressed file size with overhead:\t"
      << bstream_vector.size() / 8
      << " bytes\n";

  double compression_rate = 1;
  compression_rate -= (double)bstream_vector.size() /
                      (double)this->file_content_.size();

  compression_rate *= 100;
  std::cout
      << "Brute Compression rate: "
      << compression_rate
      << "%\n";
}

void Huffman::Decoder::DecompressFromFile(std::string file_name)
{
  // Empty Bitstream object
  Bitstream bstream = Bitstream(file_name);
  bool bit;

  // Inserts file content in decoder buffer
  while (bstream.hasBits())
  {
    bit = bstream.readBit();
    this->encoded_content_buffer_.push_back(bit);
  }

  if (DEBUG)
  {
    std::cout << "-----------------------------\n"
              << "-- Buffer Read from .huff ---\n"
              << "-----------------------------\n";
    for (auto x : this->encoded_content_buffer_)
    {
      std::cout << to_string(x);
    }
    std::cout << "\n\n";
  }
}

void Huffman::Decoder::Decode()
{
  // How many symbols are in the header
  int n_symbols = 0;

  // Current bit set to read
  // in the encoded_content_buffer
  int current_bit = 0;

  // Symbol read from the header
  std::string symbol;

  // Size read from header
  int symbol_size = 0;

  // Encode read from the header
  std::string encode;

  // Encode read from the header
  std::map<std::string, std::string> code_to_symbol;

  // Gets the number of symbols
  // contained in the header
  for (int i = 0; i < 16; i++)
  {
    // Converts the 2 byte binary
    // n_symbols form to a decimal count
    n_symbols |= ((this->encoded_content_buffer_[current_bit]) << (15 - i));
    current_bit++;
  }

  // Parses header to an map
  // of code to original symbol, for each symbol
  while (n_symbols--)
  {
    symbol.clear();
    encode.clear();
    symbol_size = 0;

    // Read symbol
    for (int i = 0; i < 8; i++)
    {
      symbol = symbol +
               std::to_string(this->encoded_content_buffer_[current_bit]);
      current_bit++;
    }

    // Read size
    for (int i = 0; i < 8; i++)
    {
      // Converts the 2 byte binary
      // symbol size form to a decimal count
      symbol_size |= ((this->encoded_content_buffer_[current_bit]) << (7 - i));
      current_bit++;
    }

    // Read code
    for (int i = 0; i < symbol_size; i++)
    {
      encode = encode +
               std::to_string(this->encoded_content_buffer_[current_bit]);
      current_bit++;
    }

    code_to_symbol[encode] = symbol;
  }

  this->current_bit_ = current_bit;
  this->code_to_symbol_ = code_to_symbol;

  if (DECODE_DEBUG)
  {
    std::cout << "Symbols Number: "
              << n_symbols << "\n";

    std::cout << "-----------------------------\n"
              << "------- Code to Symbol ------\n"
              << "-----------------------------\n";

    for (auto const &x : code_to_symbol)
    {
      std::cout << x.first
                << ':'
                << x.second
                << std::endl;
    }
  }
}

void Huffman::Decoder::DecompressHuffmanCode()
{
  // Express the current code expression
  // read from then compressed file
  std::string code = "";

  // Bit read from the compressed file
  std::string bit = "";

  // Single bit read from the file
  std::string read_from_file = "";

  std::map<std::string, std::string>::iterator it = this->code_to_symbol_.end();

  while (this->current_bit_ <
         this->encoded_content_buffer_.size())
  {
    read_from_file = std::to_string(this->encoded_content_buffer_[this->current_bit_]);
    code = code + read_from_file;

    // Theres a corresponding code read from the
    // compressed file?
    if ((it = this->code_to_symbol_.find(code)) != this->code_to_symbol_.end())
    {
      for (auto character : it->second)
      {
        bit = character;
        this->decompressed_content_buffer.push_back(stoi(bit));
      }
      code.clear();
    }

    this->current_bit_++;
  }

  if (DECODE_DEBUG)
  {
    std::cout << "-----------------------------\n"
              << "----- Decompressed content --\n"
              << "-----------------------------\n";
    for (auto bit : this->decompressed_content_buffer)
    {
      std::cout << bit;
    }
    std::cout << "\n";
  }
}

void Huffman::Decoder::DecompressToFile(std::string file_name)
{

  // Empty Bitstream object
  Bitstream bstream;

  for (auto bit : this->decompressed_content_buffer)
  {
    bstream.writeBit(bit);
  }

  //Grava o bitstream no arquivo.
  bstream.flushesToDecompressedFile(file_name);
}
