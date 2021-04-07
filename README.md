# C-LZ77-Code
A LZ77 compression algorithm implemetation in C++, to the discipline Signal Compression at Universidade de Brasilia

## Execution

To execute this program, just compile the project by:
```
$ make clean && make
```
And execute:
```
$ ./LZ77 ArquivosParaComprimir/<input_filename> <output_filename>
```
# Results

After executing the above command and running the program, two files will be 
generated acording to the input file: output_filename.lz77(compressed file) and 
output_filename.decompressed (the decompressed file from the out.lz77, that should 
be equal to `<input_file>`)
