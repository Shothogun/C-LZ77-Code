#include "../include/lz77.h"
#include "../include/huffman.h"
#include "../include/bitstream.h"
#include "errno.h"

#define FOR 0
#define DEBUG 0
#define DEBUG_DECODE 0
#define TRIPLES_DEBUG 0
#define DEBUG_COMPRESSED_BSTREAM 0

void LZ77::Encoder::CountSymbol(std::string character)
{
  this->symbol_table_[character]++;
}

std::map<std::string, double> LZ77::Encoder::GetSymbolTable()
{
  return this->symbol_table_;
}

void LZ77::Encoder::ComputeProbabilityTable()
{
  for (auto &x : this->GetSymbolTable())
  {
    this->symbol_table_[x.first] /= this->file_content_.size();
  }
}

void LZ77::Encoder::FillBuffer(std::string file_path)
{
  // Read file as binary data
  std::ifstream f(file_path);

  // File's input line read
  std::string line_read;

  // Single character from file content
  std::string character;

  char c;

  // Error, file no found
  if (!f)
  {
    std::cout << "File not found\n";
    exit(0);
  };

  // Initializes file content as empty
  this->file_content_ = "";

  while (f.get(c))
  {
    this->file_content_ += c;
  }
  f.close();

  // Fills Symbol table
  for (auto &c : this->file_content_)
  {
    character = c;
    this->CountSymbol(character);
    character.clear();
  }

  this->ComputeProbabilityTable();

  this->FlushProbabilityTableAsCSV();

#if DEBUG
  std::cout << "-----------------------------\n"
            << "-- Symbol Probabilities -----\n"
            << "-----------------------------\n";

  // As bytes
  for (auto const &x : this->GetSymbolTable())
  {
    std::cout << x.first
              << ':'
              << x.second
              << std::endl;
  }
#endif
}

void LZ77::Encoder::Encode()
{
  int offset, length;

  int const last_one = this->file_content_.size() - 1;

  // Transmited symbol's index
  int next_symbol_index;

  triple_struct triple;

  std::string symbol = "";
  this->look_ahead_buffer_ = "";

#if FOR
  for (this->current_character_index_ = 0;
       this->current_character_index_ <
       4;
       this->current_character_index_++)
#else
  for (this->current_character_index_ = 0;
       this->current_character_index_ <
       this->file_content_.size();
       this->current_character_index_++)
#endif
  {
#if DEBUG
    std::cout << "Current index:" << this->current_character_index_ << std::endl;
#endif

    // If lookahead exceed bitstream end
    if (this->current_character_index_ + this->look_ahead_buffer_size_ - 1 >
        this->file_content_.size() - 1)
    {
      // Updates lookahead
      this->look_ahead_buffer_ =
          this->file_content_.substr(this->current_character_index_,
                                     this->file_content_.size() - this->current_character_index_);
    }

    else
    {
      // Updates lookahead
      this->look_ahead_buffer_ =
          this->file_content_.substr(this->current_character_index_,
                                     this->look_ahead_buffer_size_);
    }

    // std::cout << this->look_ahead_buffer_ << std::endl;

    // Matching process
    std::tie(offset, length) = this->MatchPattern();

    next_symbol_index = this->current_character_index_ + length;

    // Transmit encode
    this->output_encoding_ += std::to_string(offset) + std::to_string(length);

    this->offset_sequence_buffer_.push_back(offset);
    this->length_sequence_buffer_.push_back(length);

    if (next_symbol_index >= this->file_content_.size())
    {
      symbol = this->file_content_[last_one];
    }

    // If there's not match
    else if (length == 0 && offset == 0)
    {
      symbol = this->file_content_[this->current_character_index_];
    }

    // Otherwise, transmits the next sequence following the
    // pattern matching
    else
    {
      symbol =
          this->file_content_[next_symbol_index];
    }

    this->output_encoding_ += symbol;

    // Next symbol index
    this->current_character_index_ = next_symbol_index;

    triple_struct triple = {offset, length, symbol};
    this->triples_vector_.push_back(triple);

#if DEBUG
    std::cout << "<" << offset << "," << length << ",";
    std::cout << symbol << ">\n";
    std::cout << "Search Buffer tree:"
              << "\n";
    for (auto const &a : this->search_buffer_tree_)
    {
      std::cout << a
                << std::endl;
    }
#endif
  }
#if TRIPLES_DEBUG
  std::cout << "Output Triples:"
            << "\n";
  for (auto const &a : this->triples_vector_)
  {
    std::cout << "Offset:"
              << a.offset
              << "\nLength:"
              << a.length
              << "\nSymbol:"
              << a.codeword
              << std::endl;
  }
#endif
}

std::tuple<int, int> LZ77::Encoder::MatchPattern()
{
  // The current match sequence in the lookahead buffer
  // being analyzed
  std::string current_match_sequence = "";

  // Match sequence in the search buffer absolute position.
  // This position is the first character from the sequence
  // position index in the file_content_.
  int match_position;

  // Offset returning value
  int offset = 0;

  // Length returning value
  int length = 0;

  // Search matching
  std::tie(offset, length) = this->SearchMatching();

  // Updates the Binary Search Tree
  this->UpdateSearchBufferTree(length);

  return std::make_tuple(offset, length);
}

void LZ77::Encoder::UpdateSearchBufferTree(int length)
{
  std::string next_to_delete;

  // Next node to be added
  // Current index plus
  std::string add_node;

  // How many nodes exceeding are occupying
  // the search buffer tree
  int n_exceed_nodes;

  for (int i = this->current_character_index_;
       i < this->current_character_index_ + length + 1; i++)
  {
    if (i >= this->file_content_.size())
    {
      break;
    }

    add_node = this->file_content_.substr(
        i,
        this->look_ahead_buffer_size_);

    // Inserts in the list of nodes to be deleted
    this->nodes_to_exclude.push_back(add_node);

    // Inserting the node
    this->search_buffer_tree_.insert(add_node);
    sequence_position_[add_node] = i;
  }

  if (this->search_buffer_tree_.size() > this->search_buffer_size_)
  {
    n_exceed_nodes = this->search_buffer_tree_.size() - this->search_buffer_size_;
    for (int i = 0;
         i < n_exceed_nodes;
         i++)
    {

      next_to_delete = this->nodes_to_exclude.front();
      this->nodes_to_exclude.erase(this->nodes_to_exclude.begin());
      this->search_buffer_tree_.erase(next_to_delete);
    }
  }
}

std::tuple<int, int> LZ77::Encoder::SearchMatching()
{
  // Offset returning value
  int offset = 0;

  // Length returning value
  int length = 0;

  // First character char from the
  // matching string sequence from search buffer
  std::string match_first_character;

  // First character always transmit <0,0,symbol>
  if (this->search_buffer_tree_.size() == 0)
  {
    return std::make_tuple(0, 0);
  }

  // Matching search on tree
  else
  {
    std::string match = this->SearchBestMatch();

    // No match!
    if (match == "")
    {
      return std::make_tuple(0, 0);
    }

    // Match, computes offset and length
    else
    {
      std::tie(offset, length) = this->LargestMatch(match);
    }
  }

  return std::make_tuple(offset, length);
}

std::tuple<int, int> LZ77::Encoder::LargestMatch(std::string match_string)
{
  int match_position = this->sequence_position_[match_string];
  int offset = this->current_character_index_ - match_position;

#if DEBUG
  std::cout << match_string
            << " "
            << this->look_ahead_buffer_
            << std::endl;
#endif

  int length = 0;
  int const last_one = this->file_content_.size() - 1;

  // Compares each character from each string
  for (int i = 0; i < match_string.size(); i++)
  {
    if (this->look_ahead_buffer_[i] == match_string[i] &&
        this->current_character_index_ + length < last_one)
    {
      length++;
    }
    else
    {
      break;
    }
  }

  return std::make_tuple(offset, length);
}

std::string LZ77::Encoder::SearchBestMatch()
{
  std::string current_sequence = "";

  int match_length = 0;
  int i = this->current_character_index_;
  std::set<std::string>::iterator it, match;

  do
  {
    current_sequence += this->file_content_[i];

    it = this->search_buffer_tree_
             .lower_bound(current_sequence);

    // Found a match!
    if (it != this->search_buffer_tree_.end() &&
        (*it)[match_length] == current_sequence[match_length])
    {
      match = it;
      match_length++;
      i++;
    }

    else
    {
      break;
    }

  } while (match_length < this->look_ahead_buffer_size_ &&
           it != this->search_buffer_tree_.end());

  // No match found!
  if (match_length == 0)
  {
    return "";
  }

  else
  {
    return *match;
  }
}

void LZ77::Encoder::FlushProbabilityTableAsCSV()
{
  std::ofstream myfile;
  myfile.open("histogram.csv");
  std::vector<std::pair<double, int>> symbols;
  char single_character;
  int i = 0;

  myfile << "Symbol,Probability\n";

  for (auto const &x : this->GetSymbolTable())
  {
    i++;
    symbols.push_back(std::make_pair(x.second, i));
  }

  std::sort(symbols.begin(),
            symbols.end(),
            std::greater<std::pair<double, int>>());

  i = 0;

  for (auto const &x : symbols)
  {
    i++;
    myfile << i << "," << (x.first * 1000000) << "\n";
  }

  myfile.close();
}

void LZ77::Encoder::CompressToFile(std::string file_path)
{

  // Empty Bitstream object
  Bitstream bstream;
  std::string bit;

  // ***** Debug *******
  std::vector<bool> bstream_vector;

  Huffman::Encoder *huffman_encoder_offset = new Huffman::Encoder();
  Huffman::Encoder *huffman_encoder_length = new Huffman::Encoder();

  // These sequence buffer is the offsets and lengths values
  // written in a single concatenated string(without space).
  // The Huffman code is pass trough these buffer, encoding
  // the offset and lengths to be send.
  huffman_encoder_offset->FillBuffer(this->offset_sequence_buffer_);
  huffman_encoder_length->FillBuffer(this->length_sequence_buffer_);

  huffman_encoder_offset->ComputeHuffmanCode();
  huffman_encoder_length->ComputeHuffmanCode();

  int offset_symbol_number = huffman_encoder_offset->GetSymbolTable().size();

  // Inserts symbols number as bits
  for (int i = 0; i < 16; i++)
  {
    bstream.writeBit((offset_symbol_number >> (15 - i)) & 1);
    bstream_vector.push_back((offset_symbol_number >> (15 - i)) & 1);
  }

  std::map<std::string, std::string> offset_symbol_encode =
      huffman_encoder_offset->GetSymbolEncode();

  // Inserts array of tuples as bits
  for (auto p : offset_symbol_encode)
  {
    // Inserts Symbol
    for (int i = 0; i < 16; i++)
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

  int length_symbol_number = huffman_encoder_length->GetSymbolTable().size();

  // Inserts symbols number as bits
  for (int i = 0; i < 16; i++)
  {
    bstream.writeBit((length_symbol_number >> (15 - i)) & 1);
    bstream_vector.push_back((length_symbol_number >> (15 - i)) & 1);
  }

  std::map<std::string, std::string> length_symbol_encode =
      huffman_encoder_length->GetSymbolEncode();

  // Inserts array of tuples as bits
  for (auto p : length_symbol_encode)
  {
    // Inserts Symbol
    for (int i = 0; i < 8; i++)
    {
      bit = p.first[8 + i];
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

  // Content write
  for (auto const &triple : this->triples_vector_)
  {
    std::string offset = LZ77::IntToBinString(triple.offset, 16);
    std::string offset_code = offset_symbol_encode[offset];
    std::string length = LZ77::IntToBinString(triple.length, 16);
    std::string length_code = length_symbol_encode[length];
    std::string bit = "";

    // Inserts offset size
    for (uint i = 0; i < offset_code.size(); i++)
    {
      bit = offset_code[i];
      bstream.writeBit(stoi(bit));
      bstream_vector.push_back((stoi(bit)));
    }

    // Inserts length
    for (uint i = 0; i < length_code.size(); i++)
    {
      bit = length_code[i];
      bstream.writeBit(stoi(bit));
      bstream_vector.push_back(stoi(bit));
    }

    char char_codeword[2];
    std::strcpy(char_codeword, triple.codeword.c_str()); // or pass &s[0]

    // Inserts codeword
    for (uint i = 0; i < 8; i++)
    {
      bstream.writeBit((char_codeword[0] >> (7 - i)) & 1);
      bstream_vector.push_back((char_codeword[0] >> (7 - i)) & 1);
    }
  }
#if DEBUG_COMPRESSED_BSTREAM

  for (auto x : bstream_vector)
  {
    std::cout << x;
  }

  std::cout << "\n";
#endif

  bstream.flushesToFile(file_path);
}

void LZ77::Decoder::DecompressFromFile(std::string file_path)
{
  // Initialize current decoder reader position
  this->current_bit_ = 0;

  // Empty Bitstream object
  Bitstream bstream = Bitstream(file_path);
  bool bit;

  // Inserts file content in decoder buffer
  while (bstream.hasBits())
  {
    bit = bstream.readBit();
    this->encoded_content_buffer_.push_back(bit);
  }

#if DEBUG
  {
    std::cout << "-----------------------------\n"
              << "-- Buffer Read from .lz77 ---\n"
              << "-----------------------------\n";
    for (auto x : this->encoded_content_buffer_)
    {
      std::cout << to_string(x);
    }
    std::cout << "\n\n";
  }
#endif
}

void LZ77::Decoder::Decode(std::string option)
{
  // How many symbols are in the header
  int n_symbols = 0;

  // Current bit set to read
  // in the encoded_content_buffer
  int current_bit = this->current_bit_;

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

  int n_offset_bits, n_size_bits;

  if (option == "offset")
  {
    n_offset_bits = 16;
    n_size_bits = 8;
  }

  else if (option == "length")
  {
    n_offset_bits = 8;
    n_size_bits = 8;
  }

  else
  {
    throw std::invalid_argument("Not valid option");
  }

  // Parses header to an map
  // of code to original symbol, for each symbol
  while (n_symbols--)
  {
    symbol.clear();
    encode.clear();
    symbol_size = 0;

    // Read symbol
    for (int i = 0; i < n_offset_bits; i++)
    {
      symbol = symbol +
               std::to_string(this->encoded_content_buffer_[current_bit]);
      current_bit++;
    }

    // Read size
    for (int i = 0; i < n_size_bits; i++)
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

  if (option == "offset")
  {
    this->offset_code_to_symbol_ = code_to_symbol;
  }

  else if (option == "length")
  {
    this->length_code_to_symbol_ = code_to_symbol;
  }

  else
  {
    throw std::invalid_argument("Not valid option");
  }

#if DEBUG_DECODE
  for (auto const &x : code_to_symbol)
  {
    std::cout << "Code: " << x.first << " "
              << "Symbol: " << x.second << std::endl;
  }
#endif
}

void LZ77::Decoder::DecompressLZ77Code()
{

  enum DECODING_STATE
  {
    kOffset,
    kLength,
    kSymbol,
    kDecode
  };

  DECODING_STATE state = kOffset;

  // Express the current code expression
  // read from then compressed file
  std::string code = "";

  // Bit read from the compressed file
  std::string bit = "";

  std::map<std::string, std::string>::iterator it;

  int length, offset;

  int current_bit = 0;

  while (this->current_bit_ <
         this->encoded_content_buffer_.size())
  {
    bit = std::to_string(
        this->encoded_content_buffer_[this->current_bit_]);
    code = code + bit;
#if DEBUG_DECOMPRESS_STREAM
    std::cout << "Current bit: "
              << this->current_bit_
              << "\tcode: "
              << code
              << std::endl;
#endif
    if (state == kOffset)
    {
      // Theres a corresponding code read from the
      // compressed file?
      if ((it = this->offset_code_to_symbol_.find(code)) !=
          this->offset_code_to_symbol_.end())
      {
        std::string string_offset = it->second;
#if DEBUG_DECOMPRESS_STREAM
        std::cout << string_offset << std::endl;
#endif
        offset = 0;
        current_bit = 0;

        for (int current_bit = 0; current_bit < 16; current_bit++)
        {
          bit = string_offset[current_bit];
          offset |= (stoi(bit) << (15 - current_bit));
        }
#if DEBUG_DECOMPRESS_STREAM

        std::cout << "Offset:"
                  << offset
                  << std::endl;
#endif
        code.clear();
        state = kLength;
      }
    }

    else if (state == kLength)
    {
      // Theres a corresponding code read from the
      // compressed file?
      if ((it = this->length_code_to_symbol_.find(code)) !=
          this->length_code_to_symbol_.end())
      {
        std::string string_length = it->second;
        length = 0;
        current_bit = 0;

        for (int i = 0; i < 8; i++)
        {

          if (i >= string_length.size())
          {
            length |= (0 << (7 - i));
            continue;
          }

          bit = string_length[current_bit];
          length |= (stoi(bit) << (7 - i));
          current_bit++;
        }
#if DEBUG_DECOMPRESS_STREAM

        std::cout << "Length:"
                  << length
                  << std::endl;
#endif
        code.clear();
        state = kSymbol;
      }
    }

    else
    {

      std::string match = "";
      std::string symbol = "";
      bool bool_bit;
      int last_bit = this->decompressed_content_buffer.size() - 1;
      int replicate_begining = last_bit -
                               (offset * 8) + 1;
#if DEBUG_DECOMPRESS_STREAM

      std::cout << "Replicate_begining:"
                << replicate_begining
                << std::endl;
#endif
      // Matching writting
      if (length > 0)
      {
        // Codeword to write
        for (int i = 0; i < (length * 8); i++)
        {
          bool_bit = this->decompressed_content_buffer[replicate_begining + i];
          this->decompressed_content_buffer.push_back(bool_bit);
          match += std::to_string(bool_bit);
        }
      }

      // Codeword to write
      for (int i = 0; i < 8; i++)
      {
        bit = std::to_string(
            this->encoded_content_buffer_[this->current_bit_ + i]);
        this->decompressed_content_buffer.push_back(stoi(bit));
        symbol += bit;
      }

      this->current_bit_ += 8;

#if DEBUG_DECOMPRESS_STREAM
      std::cout << "Code:"
                << symbol
                << std::endl;
#endif

      symbol.clear();
      code.clear();
      match.clear();
      state = kOffset;
      continue;
    }

    this->current_bit_++;
  }
}

std::string LZ77::IntToBinString(int value, int string_size)
{
  std::string bin = "";
  int offset;

  for (int i = 0; i < string_size; i++)
  {
    offset = string_size - 1 - i;

    bin += std::to_string(((value >> offset) & 1));
  }

  return bin;
}

void LZ77::Decoder::DecompressToFile(std::string file_name)
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
