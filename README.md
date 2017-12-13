# PAN
Phylogenetic Asignment Namer, using SILVAdb plus out of the shelf denoising(vsearch), alignment (infernal) and placement (EPA-ng)


 Input files are supposed to stay in the folder where the library is called. The library assume that the system has:

* cutadapt
* vsearch 
* raxml/EPA-ng
* infernal

the path to get those program could be set at the start of defintion of the class in the library as all other details of the different program.
Further the python library should include:
- pandas
- biopython 
        
 This library work as standalone in in three modes:
 
1. Denoise+Align to ref+place on SILVA tree+ Name read. input are file *.fastq pairends with illumina naming convention
    python3 PanPipe.py FULL
2. Denoise only (trim+demultiplicate+cluster+filter singleton+group across sample) 
  
  python3 PanPipe.py DENOISE

3. Give a name based on NCBI taxonomy and a given jplace file based on the SILVA tree distributed with library.
  python3 PanPipe.py NAME jplacefilename
