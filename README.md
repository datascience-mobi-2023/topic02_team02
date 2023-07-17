# Topic 02, Team 02: Deep mutational scanning data analysis to reveal sequence function relationships for tumor suppressor protein p53

Contributors
----------
Dario Prifti, Enno Schäfer, Frido Petersen and Maximilian Fidlin 


Supervisor
----------
_Prof. Dominik Niopek_  ([dominik.niopek@uni-heidelberg.de](mailto:dominik.niopek@uni-heidelberg.de))  
_Jan Mathony_  ([jan.mathony@tu-darmstadt.de](mailto:jan.mathony@tu-darmstadt.de)) 

Tutor: _Benedict Wolf_ ([b.wolf@stud.uni-heidelberg.de](mailto:b.wolf@stud.uni-heidelberg.de))    


Research Question
----------
With our project, we wanted to investigate whether the DNA binding domain, as a very essential domain of p53, is less prone
to mutations that are more likely to occur by chance. To do so, we inspected the datasets on p53 for the full length and
domain-wise, among other things. Furthermore, we invented a basic method to determine each mutation's probability to answer
our research question.


Structure of this repository
-----------------------------
First of all: To run our code, one has to download the data from the link supplied at the bottom and put all the data into a folder called "DMS_data". If the folder is named differently, our code
will not be able to load the data. The "DMS_data" folder must then be a direct subdirectory of topic02_team02.

In our repository, the final notebook, that generates all important plots is found in the Documentation folder and is called [P53_DMS_Documentation.ipynb](Documentation/P53_DMS_Documentation.ipynb).
The [report (as a pdf)](Documentation/report_DMS_topic02_team02.pdf) can be found in the same folder. Within our Documentation, you will find **five sub-topics**. The first one looks at the Comparability of p53 Datasets.
The code generating the relevant plots can be found [within the Documentation folder](Documentation/backgrounddata.py) and 
for visualization purposes, that might also be used later on, code from the [Visualization folder](visualization) was used.
The code for plots on the other four topics can be found here: 
- [Data cleanup](data_cleanup) 
- [Data exploration](data_exploration)
- [Domain comparison](domain_comparison) 
- [Calculating severity scores](severity_score)

In each of these folders, the relevant functions are defined in the .py file. In most folders, there is an additional **exploratory**
folder containing all the experimental notebooks. Jupyter notebooks, that are within the sub-topic folders but not in the exploratory 
folders are mentioned in the [P53_DMS_Documentation.ipynb file](Documentation/P53_DMS_Documentation.ipynb) and contain further 
information, outlook or important additional information. 


Covering the mandatory aspects of the project
------------
Our project was supposed to contain the following elements. Here, we list which sub-topic covers which mandatory aspect: 
- **descriptive statistics** about the datasets: Comparability of p53 Datasets
- **graphical representations**: all sub-topics
- **dimension reduction** analysis (PCA, clustering or k-means): Data exploration
- **statistical tests** (t-test, proportion tests etc): Domain comparison
- We did not implement a **linear regression**, which we discussed with our tutor Benedict Wolf beforehand.


Additional files and folders
---------
- The [environment.yml file](environment.yml) contains all packages necessary to run the code  
- The [DMS_data folder](DMS_data) is a directory, each of us contributors added on their own device. It is listed in the
.gitignore file to prevent the data from being pushed to GitHub. It must be added manually by each person using our code.
The contents of that folder can be downloaded via the link supplied below
- A file containing the chemical properties of all amino acids must be named "aminoacids.csv" in order to guarantee that
all plots can be generated in the Documentation file. It can be downloaded by the link provided at the bottom


Download the datasets worked on
----------
- The DMS data belonging into "DMS_data" can be downloaded [here](https://heibox.uni-heidelberg.de/d/d8754d7929d145efb9be)
- To calculate the severity score for a protein of interest, one has to download the dna sequence of the POI and change 
the variable "dna_sequence" to the desired sequence [here](severity_score/aa_prob.py)
- The file on chemical properties of amino acids can be downloaded via [this link](https://www.kaggle.com/datasets/alejopaullier/aminoacids-physical-and-chemical-properties?resource=download) and must be renamed to "aminoacids.csv"
and must be put into the "DMS_data" folder as well