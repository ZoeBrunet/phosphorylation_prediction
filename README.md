# Prediction of active protein phosphorylation sites using machine learning (student project)

Cells store and transmit information via post-translationa modification (PTM) of proteins. A prominent example of PTMs is reversible phosphorylation of serine, threonine or tyrosine side chains catalyzed by hundreds of protein kinases in human. Studies using low- and high-throughput techniques have reported hundreds of thousands of phosphorylation sites located in approximately 45% of eukaryotic proteomes. Phosphorylation regulates various cellular processes and its deregulation has been associated to diseases such as cancers or diabetes. As only a small percentage of human substrates and modification sites are linked to responsible protein kinases, novel computer algorithms able to to identify kinase-substrate relationships are required for the better understanding of such deleterious phosphorylation-dependent regulations. The work will focus on the analysis of homolog protein sequences in a broad set of organisms (from bacteria to mammals) to identify phosphorylation sites. Various parameters will be investigated for relevant scoring schemes and machine learning classification (e.g. site-specific conservation, window-based conservation, conservation across subsets of closely related species, and local sequence context). For validation, the results will be contrasted to unpublished experimental datasets of kinase- substrate relationships. We expect results to outperform a preliminary approach.


## Warning 

This project is a student project. This is a beta version. 


## Prerequisites

Make sure you have Python3. You can find all the needed packages in requierement.txt, just run :
```console
foo@bar:~$ pip install -r requirements.txt
```

## Dataset

To produce our training set we use dump of dbPTM available on <a href="http://dbptm.mbc.nctu.edu.tw/download.php">http://dbptm.mbc.nctu.edu.tw/download.php</a>
We have several intermediate files to store information.

| Name          |     info           |           needed file              |
| :------------- |    :-------------   | :-------------                    | 
| index      | Sequence and position of positive sites in each protein   | dump of dbPTM  |
| filter_index      | Index without redundant proteins   | index  |
| index_neg      | Sequence and position of negative sites in each protein   | filter_index for each studied residue  |
| final_index      | Sequence and position of negative and positive sites in each protein   | filter_index and index_neg  |
| dataset      | Features for each site   | final_index  |

### Create dataset 

#### How to use it

Use create_dataset.py to create csv dataset needed for model training.

If necessary the program create an data, fastas, metazoa, non_metazoa, sorted_fastas, align and csv directory. The index and the dataset are stored in data/csv/pattern/
```console
foo@bar:~$ python create_dataset.py pattern1,pattern2 file1,file2
```

#### Parameters

|Name          |     type           |           description              | Default value|
| :------------- |    :-------------   | :-------------                    | :------------- |
| pattern       | list(regexpr) | Amino acid sequence you want to detect | |
| file          | list(string)   | dump of dbPTM     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |
| --nthread          | int (optional)  | Number of thread to execute program | 1 |
| --color    | bool (optional)     | Enable color in console output | False |
| --ortholog    | bool (optional)     | Show the orthologs in the window | False |
| --species    | list(int) (optional)   | Species you want to include in your dataset | All |


#### Example

```console
foo@bar:~$ python create_dataset.py --ortholog --nthread 4 H,T path/2/phosphorylation_prediction/data_for_test/csv/data_H.csv,path/2/phosphorylation_prediction/data_for_test/csv/data_T.csv --color
-0.92
```

### Create index

### Final index to dataset

### Run muscle

## Training set and Benchmark

### Split dataset

## Machine learning

### Create models

## Model validation

### Compare tool

## Figures

### Piechart

### Boxplot

### Model info

## Features

You can test each feature we use to enrich our dataset.

### ACH

#### How to use it

Use ACH.py to get Average Cumulative Hydrophobicity of a protein sequence. 

We use hydrophobicity index proposed by Sweet and Eisenberg

The program return a float
```console
foo@bar:~$ python ACH.py sequence
```

#### Parameters

| Name          |     type           |           description              | 
| :------------- |    :-------------   | :------------- | 
| seq       | string | Amino acid sequence | 

#### Example

```console
foo@bar:~$ python ACH.py MTPTTPRLS
-0.92
```

### Freq

#### How to use it

Use freq_pattern.py to get frequency of a pattern.

If necessary the program create an align directory in which one you will find the multiple alignment. The name of the output file will be inputname_align.fasta. It also return a list of float which is the frequency of the choosen patern for each position in the alignment
```console
foo@bar:~$ python freq_pattern.py pattern file
```

#### Parameters

| Name          |     type           |           description              |
| :------------- |    :-------------   | :-------------                    | 
| pattern       | regular expression | Amino acid sequence you want to detect | 
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare     | 


#### Example

```console
foo@bar:~$ python freq_pattern.py S path/2/phosphorylation_prediction/data_for_test/fastas/example.fasta
[0, 0, 0, 0, 0, 0.8333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3333333333333333, 0, 0, 0]
```

## Information Content 

### How to use it

Use information_content.py to get Information Cotent of a fasta file
We use this formula :

![IC](https://latex.codecogs.com/gif.latex?IC%3D%20%5Csum%5Climits_%7Bj%3D1%7D%5E%7Bmax%5C_window%7D%5Csum%5Climits_%7Bi%3D1%7D%5E%7B20%7D%20p_%7Bij%7D%20log_%7B10%7D%28p_%7Bij%7D%20*%2020%29)

Where P<sub>ij</sub> is the frequency of a particular amino acid i in the j-th column

If necessary the program create an align directory in which one you will find the multiple alignment. The name of the output file will be inputname_align.fasta. It also return a float which is the information content for the alignment.
```console
foo@bar:~$ python information_content.py file
```

### Parameters

| Name          |     type           |           description              |
| :------------- |    :-------------   | :-------------                    | 
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare     | 

### Example

```console
foo@bar:~$ python information_content.py path/2/phosphorylation_prediction/data_for_test/fastas/example.fasta
21.6141621720123
```

### Shanon entropy

#### How to use it

Use shanon_entropy.py to get Shanon Entropy of a fasta file
We use this formula :

![IC](https://latex.codecogs.com/gif.latex?-%5Csum%5Climits_%7Bi%3D1%7D%5E%7B20%7D%20p_i%20log_2%28p_i%29)

Where p<sub>i</sub> is the frequency of a particular amino acid.

If necessary the program create an align directory in which one you will find the multiple alignment. The name of the output file will be inputname_align.fasta. It also return a list of float which is the shanon entropy for each position in the alignment.

```console
foo@bar:~$ python shanon_entropy.py file
```

#### Parameters 

| Name          |     type           |           description              |
| :------------- |    :-------------   | :-------------                    | 
| file          | absolute path to fasta file   | Sequence of orthologs protein you want to compare  |


#### Example

```console
foo@bar:~$ python shanon_entropy.py path/2/phosphorylation_prediction/data_for_test/fastas/example.fasta
[0.6500224216483541, 0.6500224216483541, 0.6500224216483541, 0.6500224216483541, 
0.6500224216483541, 0.6500224216483541, 0.6500224216483541, 0.7219280948873623, 
0.7219280948873623, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```

## License

This project is licensed under the GNU General Public License.
