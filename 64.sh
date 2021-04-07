make clean&&make&& ./LZ77 ArquivosParaComprimir/dom_casmurro.txt dom_casmurro >resultsLB/dom_casmurro_64
./LZ77 ArquivosParaComprimir/fonte.txt fonte >resultsLB/fonte_64
./LZ77 ArquivosParaComprimir/fonte0.txt fonte0 >resultsLB/fonte_0_64
./LZ77 ArquivosParaComprimir/fonte1.txt fonte1 >resultsLB/fonte_1_64
./LZ77 ArquivosParaComprimir/TEncEntropy.txt TEncEntropy >resultsLB/TEncEntropy_64
./LZ77 ArquivosParaComprimir/TEncSearch.txt TEncSearch >resultsLB/TEncSearch_64
mv *.lz77 resultsCompression
rm *.decompressed   