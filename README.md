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

If necessary the program create an data, fastas, metazoa, non_metazoa, sorted_fastas, align and csv directory. The index and the dataset are stored in data/csv/pattern/.
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
foo@bar:~$ python create_dataset.py H path/2/phosphorylation_prediction/data_for_test/csv/data_H.csv --species 559292
```

### Create index

#### How to use it

Use create_index.py to create csv index needed for creation of dataset.

If necessary the program create an data, fastas and csv directory. The index and the dataset are stored in data/csv/pattern/. The program will create index, filtered_index, neg_index and final_index
```console
foo@bar:~$ python create_index.py pattern1,pattern2 file1,file2
```

#### Parameters

|Name          |     type           |           description              | Default value|
| :------------- |    :-------------   | :-------------                    | :------------- |
| pattern       | list(regexpr) | Amino acid sequence you want to detect | |
| file          | list(string)   | dump of dbPTM     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |
| --nthread          | int (optional)  | Number of thread to execute program | 1 |
| --species    | list(int) (optional)   | Species you want to include in your dataset | All |


#### Example

```console
foo@bar:~$ python create_index.py --ortholog --nthread 4 H,T path/2/phosphorylation_prediction/data_for_test/csv/data_H.csv,path/2/phosphorylation_prediction/data_for_test/csv/data_T.csv
foo@bar:~$ python create_index.py H path/2/phosphorylation_prediction/data_for_test/csv/data_H.csv --species 559292
```

### Final index to dataset

#### How to use it

Use final_index_to_dataset.py to create csv dataset from final_index.

If necessary the program create an  metazoa, non_metazoa, sorted_fastas and align directory. The dataset is stored in data/csv/pattern/
```console
foo@bar:~$ python final_index_to_dataset.py pattern file
```

#### Parameters

|Name          |     type           |           description              | Default value|
| :------------- |    :-------------   | :-------------                    | :------------- |
| pattern       | regexpr | Amino acid sequence you want to detect | |
| file          | string   | final_index     | |
| max_window    | int (optional)     | Max size of the amino acid sequence in which the pattern can be find| 15 |
| --nthread          | int (optional)  | Number of thread to execute program | 1 |
| --color    | bool (optional)     | Enable color in console output | False |
| --ortholog    | bool (optional)     | Show the orthologs in the window | False |



#### Example

```console
foo@bar:~$ python final_index_to_dataset.py --ortholog --nthread 4 H path/2/phosphorylation_prediction/data_for_test/csv/H/final_index_H.csv
```

### Run muscle

Use run_muscle.py to align proteines sequence in a fasta file.

If necessary the program create an align directory. The fasta is stored in align
```console
foo@bar:~$ python run_muscle.py file
```

#### Parameters

|Name          |     type           |           description              | 
| :------------- |    :-------------   | :-------------                    | 
| file          | string   | fasta file     | |



#### Example

```console
foo@bar:~$ python run_muscle.py path/2/phosphorylation_prediction/data_for_test/fastas/example.fasta
MUSCLE v3.8.31 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

example 6 seqs, max length 18, avg  length 12
00:00:00    23 MB(-6%)  Iter   1  100.00%  K-mer dist pass 1
00:00:00    23 MB(-6%)  Iter   1  100.00%  K-mer dist pass 2
00:00:00    25 MB(-7%)  Iter   1  100.00%  Align node       
00:00:00    25 MB(-7%)  Iter   1  100.00%  Root alignment
00:00:00    25 MB(-7%)  Iter   2  100.00%  Root alignment
00:00:00    25 MB(-7%)  Iter   3  100.00%  Refine biparts
/path/2/phosphorylation_prediction/data_for_test/align/example_align.fasta
```

## Training set and Benchmark

### Split dataset

#### How to use it

Use split_dataset.py to create training and validation set from dataset. The validation test is create with all the proteines who are not used in musite and PhosphoSVM training set.

The validation and training set are stored in the same directory than the input file. The have the same name with suffix "_benchmark.csv" and "_training.csv" there are also 2 info files to sum up the important information

```console
foo@bar:~$ python split_dataset.py dataset used_proteines 
```

#### Parameters

|Name          |     type           |           description              | Default value|
| :------------- |    :-------------   | :-------------                    | :------------- |
| dataset      | string | path to dataset | |
| used_protein          | string   | path to text file in which path to forbiden protein id are stored     | |
| -convert    | string (optional)     | If somme protein have no uniprotID you can provide a csv in which you hava manually translate id. First column is the UniprotID and second column is the id of the protein in the other file| None |




#### Example

/!\ Before using these command line, please change the absolute path in /pat/2/phosphorylation_prediction/data_for_test/forbiden_id/convert.csv
```console
foo@bar:~$ python split_dataset.py  -convert /path/2/phosphorylation_prediction/data_for_test/forbiden_id/convert.csv /path/2/phosphorylation_prediction/data_for_test/csv/T/phospho_sites_T.csv /path/2/phosphorylation_prediction/data_for_test/forbiden_id/list_used_prot.txt
```

## Machine learning

### Create models

#### How to use it

Use create_models.py to create models.

The models are all stored in the same directory. A file info.txt give some metrics for each model.

```console
foo@bar:~$ python create_models.py file directory_name
```

#### Parameters

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| file          | string   | Dataset before machine learning     | |
| directory_name          | string   | Path to the directory where models will be stored     | |
| -max_models    | int (optional)     | Number of models to test| None |
| -max_time    | int (optional)     |  The maximum runtime in seconds that you want to allot in order to complete the model | 3600 |
| -max_mem_size    | string (optional)     |  the maximum size, in bytes, of the memory allocation pool to H2O. This value must a multiple of 1024 greater than 2MB. Append the letter m or M to indicate megabytes, or g or G to indicate gigabytes.  | 1g |


#### Example

```console
foo@bar:~$ python create_models.py -max_time 60 /path/2/phosphorylation_prediction/data_for_test/csv/T/phospho_sites_T_training.csv /path/2/phosphorylation_prediction/data_for_test/models
```

## Model validation

### Compare tool

#### How to use it

Use compare_tools.py to get performance of models.

If necessary the program create a prediction directory in which you will find all the results

```console
foo@bar:~$ python compare_tool.py benchmark model_directory
```

#### Parameters

| Name          |     type           |           description              | Default value|
| ------------- |    -------------   | -------------                    | :-------------: |
| benchmark          | string   | Validation set     | |
| model_directory          | string   | Path to the directory where models are stored     | |
| --musite    | bool     | Enable this bool to have Musite prediction| False |
| --rfp    | bool     | Enable this bool to have Musite prediction| False |
| -step    | int (optional)     |  Step for threshold | 100 |
| -models   | list(string) (optional)     |  List of models you want to compare | best of each family |


#### Example

```console
foo@bar:~$ python compare_tools.py /path/2/phosphorylation_prediction/data_for_test/csv/T/phospho_sites_T_benchmark.csv /path/2/phosphorylation_prediction/data_for_test/models
foo@bar:~$ python compare_tools.py /path/2/phosphorylation_prediction/data_for_test/csv/T/phospho_sites_T_benchmark.csv /path/2/phosphorylation_prediction/data_for_test/models -models  GBM_grid_0_AutoML_20180917_124949_model_46,XRT_0_AutoML_2
```

### model info

#### How to use it

Use print_models.py to get info and variable importance of models.

The program create an info file in which you find the parameters of the model. Il also plot the variable importance of the model for GBM, GLM, RDF and XRT

```console
foo@bar:~$ python print_models.py model
```

#### Parameters

| Name          |     type           |           description              |
| ------------- |    -------------   | -------------                    | 
| model          | string   | Path to the model     | 


#### Example

```console
foo@bar:~$ python print_models.py /path/2/phosphorylation_prediction/data_for_test/models/GBM_grid_0_AutoML_20180917_124949_model_46,XRT_0_AutoML_2
```

## Figures

### Piechart

#### How to use it

Use pie-chart.py to plot pie chart with index or dataset

```console
foo@bar:~$ python pie-chart.py feature file filename
```

#### Parameters

| Name          |     type           |           description              | 
| ------------- |    -------------   | -------------                    | 
| feature          | string   | Name of the column to plot | 
| files          | string   | Path to the csv  | 
| filename   | string     | Path to the plot | 

#### Example

```console
foo@bar:~$ python pie-chart.py phosphorylation_site /path/2/phosphorylation_prediction/data_for_test/csv/T/phospho_sites_T.csv /path/2/phosphorylation_prediction/data_for_test/figures/pie_chart_example.html 
```

### Boxplot

#### How to use it

Use boxplot.py to plot boxplot with index or dataset

```console
foo@bar:~$ python boxplot.py feature file filename
```

#### Parameters

| Name          |     type           |           description              | 
| ------------- |    -------------   | -------------                    | 
| feature          | string   | Name of the column to plot | 
| files          | list(string)   | Path to csv  | 
| filename   | string     | Path to the plot | 
| -names | list(string) (optional) | names for each boxplot |
#### Example

```console
foo@bar:~$ python boxplot.py nb_orthologs /path/2/phosphorylation_prediction/data_for_test/csv/T/phospho_sites_T.csv,/path/2/phosphorylation_prediction/data_for_test/csv/H/phospho_sites_H.csv /path/2/phosphorylation_prediction/data_for_test/figures/boxplot_example.html -names orthologs_T,orthologs_H
```

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
