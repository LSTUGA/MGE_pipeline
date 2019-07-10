# Introduction
This is a pipeline for the detection of mobile genetic elements (MGEs) from complete genomes or genome assemblies. The types of detected MGEs include prophage and plasmid. This pipeline integrates two tools (i.e, ProphET and Phigaro) for the detection of prophage and one tool (i.e., MOB-suite) for the detection of plasmid.

# Dependencies and installation instructions
1. Python 2.7 and Python 3.4
2. NCBI [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+) 2.7.1+
3. [seqtk](https://github.com/lh3/seqtk) 1.3-r106
4. [Prokka](https://github.com/tseemann/prokka) 1.12
5. [MOB-suite](https://github.com/phac-nml/mob-suite) 1.4.9 (dependencies: blast, circlator, mash)
```
# install all the required dependencies
$ pip3 install circlator #(dependencies: BWA, prodigal, SAMtools, MUMmer, Canu/SPAdes)
$ pip3 install mob_suite
$ mob_init
```
6. [Phigaro](https://github.com/lpenguin/phigaro) 0.1.5.0 (dependencies: MetaGeneMark, HMMER, locate)
```
# install all the required dependencies
$ sudo -H pip install phigaro
# please create and modify the configuration file ~/.phigaro/config.yml according to Phigaro's instruction.
```
7. [ProphET](https://github.com/facebook/prophet) 0.5.1 (dependencies: EMBOSS, BEDTools, BLAST, Perl modules)
```
# install all the required dependencies and Perl modules
$ git clone https://github.com/jaumlrc/ProphET.git
$ cd ProphET
# Change download parameters in the file: /user/path/to/ProphET/UTILS.dir/fetch_genomes_based_on_taxid.pl
  my $DOWNLOAD_INCREMENTS = 200;
$ ./INSTALL.pl
# Add ProphET to ~/.bashrc
  export PATH=$PATH:/user/path/to/ProphET
# Add GFF lib to ~/.bashrc
  export PERL5LIB=$PERL5LIB:/home/shaoting/tools/ProphET/UTILS.dir/GFFLib
  export PATH=$PATH:/home/shaoting/tools/ProphET/UTILS.dir/GFFLib
```
# Usage
```
mge_pipeline.py -i <assembly_path> -o <output_path> -t <threads> -c <minimum_coverage> -p <minimum_ident>
-i: path to input assembly
-o: path to output directory
-t: threads
-c: minimum coverage of each blast hit for MGE verification, default 50
-p: minimum identical percentage of each hit for MGE verification, default 90
--db_prophage: path to prophage database
--db_plasmid: path to plasmid database
--check: use '--check' flag to check the required dependencies
```

# Database
**Option 1** \
Default database for prophages: http://phaster.ca/downloads/z_DNA_fragment_DB.gz \
Download and uncompress this database, put it into the folder /user/path/to/mge_pipeline/database, and then rename it as "PHASTER_prophage_database.fna" \
Default database for plasmids: ftp://ftp.ncbi.nih.gov/refseq/release/plasmid \
Download and uncompress this database, put it into the folder /user/path/to/mge_pipeline/database, and then rename it as "Refseq_plasmid_database.fna"

**Option 2** \
Specify your own database for either prophages or plasmids using --db_prophage or --db_plasmid. The database should be in fasta format.

# Output description
Contig_ID: contig ID of the genome assembly \
Start_pos: start position of the detected MGE \
End_pos: end position of the detected MGE \
MGE_length: length of the detected MGE \
MGE_type: type of MGE \
Method: method used for MGE detection \
Database: database used for MGE verification \
Blast_hit: blast hit in the database \
Hit_identity: identical percentage of the blast hit \
Alignment_length: alignment length of the blast hit \
Query_coverage: query coverage of the blast hit \
Hit_length: length of the hit in the database \
