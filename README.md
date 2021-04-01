# C-LZ77-Code
A LZ77 compression algorithm implemetation in C++, to the discipline Signal Compression at Universidade de Brasilia

## Execution

To execute this program, just compile the project by:
```
$ make clean && make
```
And execute:
```
$ ./LZ77 ArquivosParaComprimir/<input_file>
```
# Results

After executing the above command and running the program, two files will be 
generated acording to the input file: out.lz77(compressed file) and decompressed
(the decompressed file from the out.huff, that should be equal `<input_file>`)
