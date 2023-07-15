# Topic 02, Team 02: Deep mutational scanning data analysis to reveal sequence function relationships for tumor suppressor protein p53

Contributors
----------
Dario Prifti, Enno Schäfer, Frido Petersen and Maximilian Fidlin 


Supervisor
----------
_Prof. Dominik Niopek_  ([dominik.niopek@uni-heidelberg.de](mailto:dominik.niopek@uni-heidelberg.de))  
_Jan Mathony_  ([jan.mathony@tu-darmstadt.de](mailto:jan.mathony@tu-darmstadt.de)) 

Tutor: _Benedict Wolf_ ([b.wolf@stud.uni-heidelberg.de](mailto:b.wolf@stud.uni-heidelberg.de))    


Introduction
------------


Research Question
----------


Structure of this repository
-----------------------------
First of all: To run our code, one has to download the data from the link supplied at the bottom and put all the data into a folder called "DMS_data". If the folder is named differently, our code
will not be able to load the data. The "DMS_data" folder must then be a direct subdirectory of topic02_team02.

In our repository, the final notebook, that generates all important plots is found in the Documentation folder and is called [P53 - DMS Documentation.ipynb](Documentation/P53_DMS_Documentation.ipynb).
Within our Documentation, you will find **five sub-topics**. The first one looks at the Comparability of p53 Datasets.
The code generating the relevant plots can be found [within the Documentation folder](Documentation/backgrounddata.py) and 
for visualization purposes, that might also be used later on, code from the [Visualization folder](visualization) was used.
The code for plots on the other four topics can be found here: 
- [Data cleanup](data_cleanup) 
- [Data exploration](data_exploration)
- [Domain comparison](domain_comparison) 
- [Calculating severity scores](severity_score)

In each of these folders, the relevant functions are defined in the .py file. In most folders, there is an additional **exploratory**
folder containing all the experimental notebooks. Jupyter notebooks, that are within the sub-topic folders but not in the exploratory 
folders are mentioned in the [P53 - DMS Documentation.ipynb file](Documentation/P53_DMS_Documentation.ipynb) and contain further 
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
The contents of that folder can be downloaded via the link supplied below.


Download the datasets worked on
----------
https://heibox.uni-heidelberg.de/d/d8754d7929d145efb9be/