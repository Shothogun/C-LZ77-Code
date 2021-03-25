#include "../include/lz77.h"
#define DEBUG 1

void LZ77::Encoder::FillBuffer(std::string file_path)
{
  // Read file as binary data
  std::ifstream f(file_path);

  // File's input line read
  std::string line_read;

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
}

void LZ77::Encoder::Encode()
{
  int offset, length;

  // Transmited symbol's index
  int next_symbol_index;

  std::string symbol = "";
  this->look_ahead_buffer_ = "";

#if DEBUG
  for (this->current_character_index_ = 0;
       this->current_character_index_ <=
       3;
       this->current_character_index_++)
#else
  for (this->current_character_index_ = 0;
       this->current_character_index_ <=
       this->file_content_.size() - this->look_ahead_buffer_size_;
       this->current_character_index_++)
#endif
  {
#if DEBUG
    std::cout << "Current index:" << this->current_character_index_ << std::endl;
#endif

    // If lookahead don't exceed bitstream end
    if (this->current_character_index_ + this->look_ahead_buffer_size_ <=
        this->file_content_.size())
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

    // If there's not match
    if (length == 0 && offset == 0)
    {
      symbol = this->file_content_[this->current_character_index_];
    }

    // Otherwise, transmits the next sequence following the
    // pattern matching
    else
    {
      symbol =
          this->file_content_[next_symbol_index];

      this->current_character_index_ = next_symbol_index;
    }

    this->output_encoding_ += symbol;

    // Next character
    if (length > 1)
    {
      this->current_character_index_ += (1 + length);
    }

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

  // if (this->current_character_index_ <
  //     this->look_ahead_buffer_size_)
  // {
  //   // Inserts a single character in the search buffer
  //   // tree and sends a single symbol in the output
  //   current_match_sequence = this->file_content_[this->current_character_index_];

  //   auto match_sequence = this->search_buffer_tree_.find(current_match_sequence);

  //   // Found a match!
  //   if (match_sequence != this->search_buffer_tree_.end())
  //   {
  //     match_position = sequence_position_.at(*match_sequence);
  //     // Distance between current character and
  //     // match sequence position
  //     offset = this->current_character_index_ - match_position;
  //     length = 1;
  //   }

  //   this->search_buffer_tree_.insert(current_match_sequence);
  //   sequence_position_[current_match_sequence] = this->current_character_index_;
  // }

  // Search matching
  std::tie(offset, length) = this->SeachMatching();

  // Updates the Binary Search Tree
  this->UpdateSearchBufferTree();

  return std::make_tuple(offset, length);
}

void LZ77::Encoder::UpdateSearchBufferTree()
{
  std::string next_to_delete;

  // Next node to be added
  // Current index plus
  std::string add_node = this->file_content_.substr(
      this->current_character_index_,
      this->look_ahead_buffer_size_);

  // Inserts in the list of nodes to be deleted
  this->nodes_to_exclude.push_back(add_node);

  // Inserting the node
  this->search_buffer_tree_.insert(add_node);
  sequence_position_[add_node] = this->current_character_index_;

  if (this->search_buffer_tree_.size() > this->search_buffer_size_)
  {
    next_to_delete = this->nodes_to_exclude.front();

    this->nodes_to_exclude.erase(this->nodes_to_exclude.begin());
    this->search_buffer_tree_.erase(next_to_delete);
  }
}

std::tuple<int, int> LZ77::Encoder::SeachMatching()
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
      std::tie(offset, length) = this->LargestMatch();
    }
  }

  return std::make_tuple(offset, length);
}

std::tuple<int, int> LZ77::Encoder::LargestMatch()
{
  return std::make_tuple(0, 0);
}