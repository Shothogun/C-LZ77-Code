#ifndef HUFFMAN_HPP
#define HUFFMAN_HPP

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <stack>
#include <iterator>
#include <cmath>

namespace Huffman
{
  //! Coder class
  /*
    * Coder
    * 
    * Class dedicated to convert text file to create
    * the Huffman Coding compression
    */
  class Encoder
  {
  private:
    //! Symbol Table
    /*
     * Symbol table
     * 
     * Maps a symbol to its probability
    */
    std::map<std::string, double> symbol_table_;

    //! characters counter
    /*
     * Counts how many characters exists on the source
    */
    int character_counter_;

    //! Largest code size
    /*
     * Indicates the bits size require do represent
     * all codes.
    */
    int codes_size_;

    //! File string stream
    /*
     * Its the content read from the file
     * as a vector of bits
    */
    std::vector<bool> file_content_;

    //! Encoded data
    /*
     * Contains the encoded data bitstream
    */
    std::vector<bool> encoded_data_;

    //! Symbol Encode map
    /*
     * Maps the symbol to its code
    */
    std::map<std::string, std::string> symbol_encode_;

  public:
    //! Log values to print
    /*
     * Entropy and average bits per symbol rate
    */
    double entropy_;
    double average_rate_;

    //! characters counter
    /*
     * Increases n_characters on the character_counter variable.
    */
    void CountCharacters(int n_characters);

    //! Fill stream
    /*
     * Fill the Coder buffer with the file's content
    */
    void FillBuffer(std::string file_path);

    //! Fill stream
    /*
     * Fill the Coder buffer with a buffer
    */
    void FillBuffer(std::vector<int> buffer);

    //! Fill stream
    /*
     * Fill the Coder buffer with a buffer
    */
    void FillBuffer(std::vector<std::string> buffer);

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

    //! Get buffer
    /*
     * Returns buffer content
    */
    std::vector<bool> GetBuffer();

    //! Get code size
    /*
     * Returns codes size
    */
    int GetCodesSize();

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

    //! Get Symbol Encode
    /*
     * Returns Symbol Encode pair
    */
    std::map<std::string, std::string> GetSymbolEncode();

    //! Computes probability table
    /*
     * This function computes the probability table 
     * from of the symbol table after populate it with its
     * counting in the source input file.
    */
    void ComputeProbabilityTable();

    //! Flush Probability Table As CSV
    /*
     * Writes symbols and its frequency 
     * as csv file to histogram plotation
    */
    void FlushProbabilityTableAsCSV(std::string file_name);

    //! Compute Huffman Code
    /*
     * This function computes the Huffman code for a given
     * symbol alphabet and its frequency in a source.
    */
    void ComputeHuffmanCode();

    //! Write Compress File
    /*
     * Get the encoded file and writes
     * it in the compressed file.
     * Pattern:
     * - Header:
     *    2 Bytes: symbols number
     *    tuples: array of tuples [a size code]
     *    a    -> 1 Byte:  symbol
     *    size -> 1 Byte:  encode size
     *    code -> n Bytes: symbol encode
     * - Content: bitstream
    */
    void WriteCompressFile();

    //! Encode function
    /*
     * Gets the input file buffer file_content_ and encodes it
     * to the created Huffman Code buffer encoded_data_
    */
    void Encode();

    //! Get Encoded content
    /*
     * Returns a vector of the enconded content.
     * Each index is the i-th byte from the original file,
     * storing the corresponding code.
    */
    std::vector<std::string> GetEncodedContent();

    //! Compress to File function
    /*
     * Writes the encoded_data_ to a 
     * compressed file, .huff
    */
    void CompressToFile(std::string file_name);
  };

  //! Decoder class
  /*
    * Decoder
    * 
    * Class dedicated to decompress .huff compressed file to 
    * the original file and format
    */
  class Decoder
  {
  private:
    //! Current Bit
    /*
     * Express the current bit 
     * read from the compressed file
    */
    int current_bit_;

    //! Encoded Content Buffer
    /*
     * .huff file's bits content buffer
    */
    std::vector<bool> encoded_content_buffer_;

    //! Decompressed file buffer
    /*
     * Out decompressed file buffer
    */
    std::vector<bool> decompressed_content_buffer;

    //! Code to symbol
    /*
     * In: Code Out: Original Symbol
    */
    std::map<std::string, std::string> code_to_symbol_;

  public:
    //! Decompress to File function
    /*
     * Get the coded file content(encoded_content_buffer_)
     * and translates it to the original decompressed 
     * decompressed_content_buffer
    */
    void DecompressHuffmanCode();

    //! Decompress from File function
    /*
     * Read the data from .huff file
     * to encoded_content_buffer_ vector
    */
    void DecompressFromFile(std::string file_name);

    //! Decompress to File function
    /* 
     * Writes the decompressed_content_buffer 
     * to a output file 
    */
    void DecompressToFile(std::string file_name);

    //! Decompress from File function
    /*
     * Read the data from .huff file
     * to encoded_content_buffer_ vector
    */
    void Decompress(std::string file_name);

    //! Decode
    /*
     * Gets the encoded_content_buffer_ bits and 
     * decodes it to its original content decompressed
    */
    void Decode();
  };
} // namespace Huffman

#endif