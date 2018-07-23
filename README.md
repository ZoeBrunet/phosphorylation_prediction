# Prediction of active protein phosphorylation sites using machine learning (student project)

Cells store and transmit information via post-translationa modification (PTM) of proteins. A prominent example of PTMs is reversible phosphorylation of serine, threonine or tyrosine side chains catalyzed by hundreds of protein kinases in human. Studies using low- and high-throughput techniques have reported hundreds of thousands of phosphorylation sites located in approximately 45% of eukaryotic proteomes. Phosphorylation regulates various cellular processes and its deregulation has been associated to diseases such as cancers or diabetes. As only a small percentage of human substrates and modification sites are linked to responsible protein kinases, novel computer algorithms able to to identify kinase-substrate relationships are required for the better understanding of such deleterious phosphorylation-dependent regulations. The work will focus on the analysis of homolog protein sequences in a broad set of organisms (from bacteria to mammals) to identify phosphorylation sites. Various parameters will be investigated for relevant scoring schemes and machine learning classification (e.g. site-specific conservation, window-based conservation, conservation across subsets of closely related species, and local sequence context). For validation, the results will be contrasted to unpublished experimental datasets of kinase- substrate relationships. We expect results to outperform a preliminary approach.


## Warning 

This project is a student project. This is a beta version. 


## Prerequisites

Make sure you have Python3. You can find all the needed packages in requierement.txt, just run :
```console
foo@bar:~$ pip install -r requirements.txt
```

## Create Dataset

### How to use it

To create a data set you will need a dump of Phospho.ELM or of dbtpm database. The programs create_training_set_phospho_ELM.py and create_training_dbtpm.py, will automatically request the orthoDB database to find ortholog of each protein.

It will return a csv files where you can find alignment score for each pull of protein. This file is in data/csv/pattern/file_pattern_phospho_sites.csv
```console
foo@bar:~$ python create_training_set_phospho_ELM.py pattern file
foo@bar:~$ python create_training_dbtpm.py pattern file
```
The dataset will be in data/csv/pattern and its name will be input_pattern_phospho_sites.csv.

### Parameters

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| pattern       | regular expression | Amino acid sequence you want to detect | |
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |
| --nthread          | int (optional)  | Number of thread to execute program | 1 |
| --color    | bool (optional)     | Enable color in console output | False |
| --ortholog    | bool (optional)     | Show the orthologs in the window | False |

### Example

```console
foo@bar:~$ python create_training_set_phospho_ELM.py --pogression --ortholog Y path/2/data/csv/sample.csv 13
Parsing csv
Extracting Y phosphorylation site
Preparing queries
import O14543 from csv file 1/1 = 100.0%

Info :

UniprotID : O14543   GeneID : 9021   TaxID : 9606   Position : 203   Metazoa : True   Phosphorylation : True

sequence            :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9986:003a73         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
246437:003f86       :N    G    H    L    G    S    Y    E    K    V    T    Q    L     
 
10160:002b02        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
1026970:004209      :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9544:00362a         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
29078:00306b        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
885580:0044d1       :N    E    R    L    D    S    Y    E    K    V    T    Q    L     
 
43179:001c9c        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
10141:000de1        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
34839:000c7e        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
51337:001ca6        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
10116:0005e0        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
10090:003af8        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
10036:001185        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
10042:000bff        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9402:002fe2         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
482537:002c70       :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
30611:0045bd        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
379532:003fbc       :N    G    H    L    D    S    Y    Q    K    V    T    Q    L     
 
61622:00417a        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9483:0047d8         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9531:00106d         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9545:000d28         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9555:000076         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9595:00389e         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9597:001dc4         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9598:001a6e         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
9601:002008         :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
54131:003c06        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
60711:004821        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
37293:0023a2        :N    G    H    L    D    S    Y    E    K    V    T    Q    L     
 
freq                :0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0 
shanon entropy      :0.0, 0.2, 0.2, 0.0, 0.2, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0 
 
information content :            6.3            ,15.4,            6.4
 
ACH                 :            -0.7            ,-1.5,            -1.0 

foo@bar:~$ python create_training_set_dbtpm.py --color Y path/2/data/csv/sample.csv 
```

## Get frequence of a pattern

### How to use it

To get the frequency of a pattern in a sequence you can run freq_pattern.py
```console
foo@bar:~$ python freq_pattern.py pattern file.fasta
```
It will return you a list with the frequency of the pattern you choose in each position of the alignment.

### Parameters

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| pattern       | regular expression | Amino acid sequence you want to detect | |
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |


### Example

```console
foo@bar:~$ python freq_pattern.py T path/2/data/align/example_align.fasta
[0, 0, 0, 0, 0.2, 0.6000000000000001, 0.4, 0.2, 0.2]
```

## Information Content 

### How to use it

We use this formula :

![IC](https://latex.codecogs.com/gif.latex?IC%3D%20%5Csum%5Climits_%7Bj%3D1%7D%5E%7Bmax%5C_window%7D%5Csum%5Climits_%7Bi%3D1%7D%5E%7B20%7D%20p_%7Bij%7D%20log_%7B10%7D%28p_%7Bij%7D%20*%2020%29)

Where P<sub>ij</sub> is the frequency of a particular amino acid i in the j-th column

To get the information content you can run information_content.py
```console
foo@bar:~$ python information_content.py pattern file.fasta
```

### Parameters

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |

### Example

```console
foo@bar:~$ python information_content.py path/2/data/align/example_align.fasta
example_align 6 seqs, max length 19, avg  length 9
00:00:00    23 MB(-6%)  Iter   1  100.00%  K-mer dist pass 1
00:00:00    23 MB(-6%)  Iter   1  100.00%  K-mer dist pass 2
00:00:00    24 MB(-7%)  Iter   1  100.00%  Align node       
00:00:00    24 MB(-7%)  Iter   1  100.00%  Root alignment
00:00:00    24 MB(-7%)  Iter   2  100.00%  Refine tree   
00:00:00    24 MB(-7%)  Iter   2  100.00%  Root alignment
00:00:00    24 MB(-7%)  Iter   2  100.00%  Root alignment
00:00:00    24 MB(-7%)  Iter   3  100.00%  Refine biparts
23.54473994433018
```

## Machine learning 

### How to use it

We use the autoML function of H2O.ai. The best models create will be in a model directory in a data directory
```console
foo@bar:~$ python create_models.py file.csv max_model
```

### Parameters

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| file          | absolute path to csv file   | Dataset before machine learning     | |
| -max_models    | int (optional)     | Number of models to test| None |
| -max_time    | int (optional)     |  The maximum runtime in seconds that you want to allot in order to complete the model | 3600 |
| -max_mem_size    | string (optional)     |  the maximum size, in bytes, of the memory allocation pool to H2O. This value must a multiple of 1024 greater than 2MB. Append the letter m or M to indicate megabytes, or g or G to indicate gigabytes.  | 1g |

### Example

```console
foo@bar:~$ python create_models.py path/2/data/csv/sample.csv -max_model 10
foo@bar:~$ python create_models.py path/2/data/csv/sample.csv -max_time 7200
```

## Add information in training set

### How to use it

To enrich your data set you may want add some information. Put your information in another csv file. Be sure that your csv has 'uniprotID' and 'position' field.
```console
foo@bar:~$ python merge_csv.py csv1 csv2
```
The output file csv1+csv2.csv will be in the same directory than csv1.

### Parameters

| Name          |     type           |           description              | 
| ------------- |    -------------   | -------------                    | 
| csv1          | absolute path to csv file   | First csv you want to merge    | 
| csv2          | absolute path to csv file   | Second csv you want to merge    | 

### Example

```console
foo@bar:~$ python merge_csv.py path/2/data/csv/csv1 path/2/data/csv/csv2
```

## Plot pie chart 

### How to use it

To get info on the dump of phospho.ELM or dbtpm you can plot pie chart. Use print_info_phospho_ELM.py or print_info_dbtpm.py
```console
foo@bar:~$ python print_info_phospho_ELM.py file column caption
foo@bar:~$ python print_info_dbtpm.py file column caption
```

### Parameters

| Name          |     type           |           description              |
| ------------- |    -------------   | -------------                    | 
| file          | absolute path to csv file   | Dump of phospho.ELM   | 
| column          | string | Name of the csv column you want to plot  | 
| caption          | string | Title of your figure  | 

### Example

```console
foo@bar:~$ python print_info_phospho_ELM.py path/2/data/csv/sample.csv column caption
foo@bar:~$ python print_info_dbtpm.py path/2/data/csv/sample.csv column caption
```

## Align orthologs

### How to use it

To align fasta file run run_muscle.py
```console
foo@bar:~$ python run_muscle.py file.fasta
```
It will create an align directory in which one you will find the multiple alignment. The name of the output file will be input_align.fasta 

### Parameters 

| Name          |     type           |           description              |
| ------------- |    -------------   | -------------                    | 
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare  |


### Example

```console
foo@bar:~$ python run_muscle.py path/2/data/fastas/example.fasta
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
path/2/data/align/example_align.fasta
```

## Compute Shanon entropy

### How to use it

To compute Shanon entropy we use the formula :

![IC](https://latex.codecogs.com/gif.latex?-%5Csum%5Climits_%7Bi%3D1%7D%5E%7B20%7D%20p_i%20log_2%28p_i%29)

Where p<sub>i</sub> is the frequency of a particular amino acid.

```console
foo@bar:~$ python shanon_entropy.py file.fasta
```
It will create an align directory in which one you will find the multiple alignment. The name of the output file will be inputname_align.fasta 

### Parameters 

| Name          |     type           |           description              |
| ------------- |    -------------   | -------------                    | 
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare  |


### Example

```console
foo@bar:~$ python shanon_entropy.py path/2/data/align/example.fasta
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2516291673878228, 0.0, 1.0, 1.0, 1.4591479170272448, 1.5219280948873621, 0.9709505944546686, 1.9219280948873623, 0.8112781244591328]
```

## License

This project is licensed under the GNU General Public License.
