# Prediction of active protein phosphorylation sites using machine learning

Cells store and transmit information via post-translationa modification (PTM) of proteins. A prominent example of PTMs is reversible phosphorylation of serine, threonine or tyrosine side chains catalyzed by hundreds of protein kinases in human. Studies using low- and high-throughput techniques have reported hundreds of thousands of phosphorylation sites located in approximately 45% of eukaryotic proteomes. Phosphorylation regulates various cellular processes and its deregulation has been associated to diseases such as cancers or diabetes. As only a small percentage of human substrates and modification sites are linked to responsible protein kinases, novel computer algorithms able to to identify kinase-substrate relationships are required for the better understanding of such deleterious phosphorylation-dependent regulations. The work will focus on the analysis of homolog protein sequences in a broad set of organisms (from bacteria to mammals) to identify phosphorylation sites. Various parameters will be investigated for relevant scoring schemes and machine learning classification (e.g. site-specific conservation, window-based conservation, conservation across subsets of closely related species, and local sequence context). For validation, the results will be contrasted to unpublished experimental datasets of kinase- substrate relationships. We expect results to outperform a preliminary approach.


## Prerequisites

Make sure you have Python3. You also have to install Biopython available on http://biopython.org/


## How to use it ?

To run the program you have to put your fasta file in the same file than scoring.py.
The pattern you want scoring can be any Python regular expression
```console
foo@bar:~$ python3 scoring.py pattern file.fasta
```
It will return you a list with the frequency of the pattern you choose in each position of the alignment.


## Parameters 

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| pattern       | regular expression | amino acid sequence you want to detect | |
| file          |       fasta file   | sequence of orthologs protein you want to compare     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which on pattern can be find| 15 |


## Example

```console
foo@bar:~$ python3 scoring.py T example.fasta
[0, 0, 0, 0, 0.2, 0.6000000000000001, 0.4, 0.2, 0.2]
```

```console
foo@bar:~$ python3 scoring.py T.A example.fasta
[0, 0, 0, 0.5333333333333333, 0.26666666666666666, 0.8, 0.13333333333333333, 0.26666666666666666, 0]
```

## Run test

Copy paste this command line 

```console
foo@bar:~$ python3 test_scoring.py
```

## License

This project is licensed under the GNU General Public License.
