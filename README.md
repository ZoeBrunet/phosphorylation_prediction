# Prediction of active protein phosphorylation sites using machine learning (student projetc)

Cells store and transmit information via post-translationa modification (PTM) of proteins. A prominent example of PTMs is reversible phosphorylation of serine, threonine or tyrosine side chains catalyzed by hundreds of protein kinases in human. Studies using low- and high-throughput techniques have reported hundreds of thousands of phosphorylation sites located in approximately 45% of eukaryotic proteomes. Phosphorylation regulates various cellular processes and its deregulation has been associated to diseases such as cancers or diabetes. As only a small percentage of human substrates and modification sites are linked to responsible protein kinases, novel computer algorithms able to to identify kinase-substrate relationships are required for the better understanding of such deleterious phosphorylation-dependent regulations. The work will focus on the analysis of homolog protein sequences in a broad set of organisms (from bacteria to mammals) to identify phosphorylation sites. Various parameters will be investigated for relevant scoring schemes and machine learning classification (e.g. site-specific conservation, window-based conservation, conservation across subsets of closely related species, and local sequence context). For validation, the results will be contrasted to unpublished experimental datasets of kinase- substrate relationships. We expect results to outperform a preliminary approach.


## Warning 

This project is a student project. This is a beta version . 


## Prerequisites

Make sure you have Python3. You also have to install Biopython available on http://biopython.org/


## How to use it ?

To align fasta file run run_muscle.py
```console
foo@bar:~$ python run_muscle.py file.fasta
```
It will create an align directory in which one you will find the multiple alignment

To get the frequency of a pattern in a sequence you can run freq_pattern.py
```console
foo@bar:~$ python freq_pattern.py pattern file.fasta
```
It will return you a list with the frequency of the pattern you choose in each position of the alignment.

To get the information content you can run information_content.py
```console
foo@bar:~$ python information_content.py pattern file.fasta
```

To create a data set you will need a dump of Phospho.ELM database. The program align_fasta_file.py, will automatically request the orthoDB database to find ortholog of each protein.
It will return a csv file where you can find alignment score for each pull of protein
```console
foo@bar:~$ python align_fasta_file.py path pattern file
```
/!\ For the moment align_fasta_file.py is a very slow program

## Parameters 


| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| pattern       | regular expression | Amino acid sequence you want to detect | |
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |


## Example

```console
foo@bar:~$ python run_muscle.py example.fasta
MUSCLE v3.8.31 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

example 6 seqs, max length 19, avg  length 9
00:00:00    23 MB(-6%)  Iter   1  100.00%  K-mer dist pass 1
00:00:00    23 MB(-6%)  Iter   1  100.00%  K-mer dist pass 2
00:00:00    24 MB(-7%)  Iter   1  100.00%  Align node       
00:00:00    24 MB(-7%)  Iter   1  100.00%  Root alignment
00:00:00    24 MB(-7%)  Iter   2  100.00%  Refine tree   
00:00:00    24 MB(-7%)  Iter   2  100.00%  Root alignment
00:00:00    24 MB(-7%)  Iter   2  100.00%  Root alignment
00:00:00    24 MB(-7%)  Iter   3  100.00%  Refine biparts
path2example_align.fasta
```

```console
foo@bar:~$ python freq_pattern.py T example.fasta
[0, 0, 0, 0, 0.2, 0.6000000000000001, 0.4, 0.2, 0.2]
```

```console
foo@bar:~$ python freq_pattern.py T.A example.fasta
[0, 0, 0, 0.5333333333333333, 0.26666666666666666, 0.8, 0.13333333333333333, 0.26666666666666666, 0]
```

```console
foo@bar:~$ python information_content.py example_align.fasta
[0, 0, 0, 0.5333333333333333, 0.26666666666666666, 0.8, 0.13333333333333333, 0.26666666666666666, 0]
18.097364729664267
```

```console
foo@bar:~$ python create_training_set.py sample.csv 
Import O14543 from the csv file
Import O14543 from the csv file
request the orthodb API for gene id 1/2 = 50.0%
request the orthodb API for gene id 2/2 = 100.0%
run muscle 1/2 = 50.0%
run muscle 2/2 = 100.0%
```

## License

This project is licensed under the GNU General Public License.
