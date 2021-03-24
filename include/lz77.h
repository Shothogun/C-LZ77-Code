#ifndef LZ77_H
#define LZ77_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <tuple>
#include <map>

namespace LZ77
{
  struct output_triple
  {
    int offset;
    int length;
    int codeword;
  };

  typedef output_triple output_triple;

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
    //! File string stream
    /* Stores the sequence string index position in the file content.
     * This sequence is the all possible sequences inside the
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
    const int search_buffer_size_ = 6;

    //! Look Ahead Buffer
    /*
     *  Look ahead buffer size value used especially
     *  in the encoding process  
    */
    const int look_ahead_buffer_size_ = 5;

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
    //! Fill stream
    /*
     * Fill the Coder buffer with the file's content
    */
    void FillBuffer(std::string file_path);

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
    void UpdateSearchBufferTree();

    //! Match Pattern function
    /*
     * Compares lookahead buffer with search buffer
     * and return the offset and length values
    */
    std::tuple<int, int> MatchPattern();
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