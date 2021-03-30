#ifndef LZ77_H
#define LZ77_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <tuple>
#include <map>

namespace LZ77
{
  struct triple_struct
  {
    int offset;
    int length;
    std::string codeword;
  };

  typedef triple_struct triple_struct;

  //! Encoder class
  /*
    * Coder
    * 
    * Class dedicated to convert text file to create
    * the LZ77 compression
    */
  class Encoder
  {
  private:
    //! Log values to print
    /*
     * Entropy and average bits per symbol rate
    */
    double entropy_;
    double average_rate_;

    //!
    /*
     * Maps a symbol to its probability
    */
    std::vector<triple_struct> triples_vector_;

    //! Symbol Table
    /*
     * Symbol table
     * 
     * Maps a symbol to its probability
    */
    std::map<std::string, double> symbol_table_;

    //! File string stream
    /* Stores the sequence string index position in the file content.
     * This sequence is one the all possible sequences inside the
     * search buffer, as known as the binary search tree node.
     */
    std::map<std::string, int> sequence_position_;

    //! File string stream
    /*
     * Its the content read from the file
     * as a vector of bits
    */
    std::string file_content_;

    //! File content encoded
    /*
     * Sequence of output triple: offset, length and symbol
    */
    std::string output_encoding_;

    //! Offset sequence buffer
    /*
     * Sequence of offsets produced on LZ77 encoding
    */
    std::vector<int> offset_sequence_buffer_;

    //! Length sequence buffer
    /*
     * Sequence of lengths produced on LZ77 encoding
    */
    std::vector<int> length_sequence_buffer_;

    //! Encoded offset codes
    /*
     * Sequence of offsets encoded from Huffman.
     * Each index correspond to the ith triple.
    */
    std::vector<std::string> encoded_offset_codes;

    //! Lengths offset codes
    /*
     * Sequence of lengths encoded from Huffman.
     * Each index correspond to the ith triple.
    */
    std::vector<std::string> encoded_length_codes;

    //! File string stream
    /*
     * Its the list of all nodes to be deleted, preserving
     * the nodes insertion order
    */
    std::vector<std::string> nodes_to_exclude;

    //! Search buffer size
    /*
     *  Search buffer size value used especially
     *  in the encoding process
    */
    const int search_buffer_size_ = 2048;

    //! Look Ahead Buffer
    /*
     *  Look ahead buffer size value used especially
     *  in the encoding process  
    */
    const int look_ahead_buffer_size_ = 255;

    //! Match position
    /*
     *  Indicates the file's index position
     *  where occurs the match 
    */
    int match_position_;

    //! Current charater index
    /*
     *  Its the current index in the file_content_ 
     *  buffer
    */
    int current_character_index_;

    // Look ahead buffer consulted on
    // symbols sequence matching
    std::string look_ahead_buffer_;

    // Search buffer consulted on
    // symbols sequence matching
    std::multiset<std::string> search_buffer_tree_;

  public:
    //! characters counter
    /*
     * Increases n_characters on the character_counter variable.
    */
    void CountCharacters(int n_characters);

    //! How many Characters
    /*
     * Returns the number of characters from read
     * from the file
    */
    int HowManyCharacters();

    //! Characters quantity
    /*
     * Returns the number of characters from read
     * from the file
    */
    int CharactersQuantity();

    //! Fill stream
    /*
     * Fill the Coder buffer with the file's content
    */
    void FillBuffer(std::string file_path);

    //! Count Symbol
    /*
     * Count the read symbol in the symbol_table
    */
    void CountSymbol(std::string character);

    //! Get Symbol Table
    /*
     * Returns Symbol Table
    */
    std::map<std::string, double> GetSymbolTable();

    //! Computes probability table
    /*
     * This function computes the probability table 
     * from of the symbol table after populate it with its
     * counting in the source input file.
    */
    void ComputeProbabilityTable();

    //! Encode function
    /*
     * Gets the input file buffer file_content_ and encodes it
     * following the LZ77 pattern <offset, length, codeword>
    */
    void Encode();

    //! Update Search Buffer Tree
    /*
     * According to the current file content index position
     * and search buffer content updates the binary search tree
     * adding the (search buffer + look ahead buffer ending) node
     * and removing (search buffer + look ahead buffer beginning)
    */
    void UpdateSearchBufferTree(int length);

    //! Match Pattern function
    /*
     * Compares lookahead buffer with search buffer
     * and return the offset and length values
    */
    std::tuple<int, int> MatchPattern();

    //! Search Matching
    /*
     * Seeks on the search buffer matching sequences
     * with current encoding character in the lookahead buffer
    */
    std::tuple<int, int> SearchMatching();

    //! Search Matching
    /*
     * Seeks on the search buffer matching sequences
     * with current encoding character in the lookahead buffer
    */
    std::tuple<int, int> LargestMatch(std::string match_string);

    //! Flush Probability Table As CSV
    /*
     * Writes symbols and its frequency 
     * as csv file to histogram plotation
    */
    void FlushProbabilityTableAsCSV();

    //! Compress to file
    /*
     * Encodes the offset and length values
     * from the triples, and writes the 
     * compressed file .lz77
    */
    void CompressToFile(std::string file_path);

    //! Search Best Match
    /*
     * Search on tree the best sequence match
     * to the current lookahead buffer
    */
    std::string SearchBestMatch();
  };

  //! Decoder class
  /*
    * Decoder
    * 
    * Class dedicated to decompress .lz77 compressed file to 
    * the original file and format
    */
  class Decoder
  {
  private:
    Decoder();

  public:
  };
} // namespace LZ77

#endif