#include "../include/lz77.h"
#include "../include/huffman.h"
#define DEBUG 0
#define TRIPLES_DEBUG 0

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

  // Error, file no found
  if (!f)
  {
    std::cout << "File not found\n";
    exit(0);
  };

  // Initializes file content as empty
  this->file_content_ = "";

  // Fill file content as a great string
  while (getline(f, line_read))
  {
    this->file_content_ += line_read;
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

#if DEBUG
  for (this->current_character_index_ = 0;
       this->current_character_index_ <
       this->file_content_.size();
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

    std::cout << "<" << offset << "," << length << ",";
    std::cout << symbol << ">\n";
    triple_struct triple = {offset, length, symbol};
    this->triples_vector_.push_back(triple);

#if DEBUG

    std::cout << "Search Buffer tree:"
              << "\n";
    for (auto const &a : this->search_buffer_tree_)
    {
      std::cout << a
                << std::endl;
    }
#endif

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

  // Current character being analysed from the
  // content file
  std::string current_character;
  current_character = this->file_content_[this->current_character_index_];

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
    auto it = std::lower_bound(
        this->search_buffer_tree_.begin(),
        this->search_buffer_tree_.end(),
        current_character);

    // If there's a match, gets a first
    // character from the match
    if (it != this->search_buffer_tree_.end())
    {
      match_first_character = (*it)[0];
    }

    // No match!
    if (it == this->search_buffer_tree_.end())
    {
      return std::make_tuple(0, 0);
    }

    // The found sequence doesn't have
    // the same first character
    else if (current_character != match_first_character)
    {
      return std::make_tuple(0, 0);
    }

    // Matching lenght computation
    else
    {
      std::tie(offset, length) = this->LargestMatch(*it);
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

void LZ77::Encoder::EncodeOffsetLength()
{
  Huffman::Encoder *huffman_encoder_offset = new Huffman::Encoder();
  Huffman::Encoder *huffman_encoder_length = new Huffman::Encoder();

  std::vector<int> offset_vector;
  std::vector<int> length_vector;

  for (auto const &triple : this->triples_vector_)
  {
    offset_vector.push_back(triple.offset);
    length_vector.push_back(triple.length);
  }

  huffman_encoder_offset->FillBuffer(offset_vector);
  huffman_encoder_length->FillBuffer(length_vector);

  huffman_encoder_offset->Encode();
  huffman_encoder_length->Encode();
}