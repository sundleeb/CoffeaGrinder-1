# CoffeeGrinder
### Prerequisites

* Python 2.7
* numpy version 1.13.1 
* jupyter 4.3+   
* ipython version 5.7+  
* histbook     
* vega version 1.1 
* uproot 

```
pip install numpy==1.13
pip install jupyter --user
pip install ipython
pip install histbook --user
pip install vega==1.1 --user
pip install uproot --user
```
Note: you might need to upgrade pip to version 18 to get Jupyter 4.3 or above  (`pip install --upgrade pip`)

### Setup
Once you have installed the prerequisites, set up the striped client
```
  git clone http://cdcvs.fnal.gov/projects/nosql-ldrd striped     
  cd striped
  python setup.py install --user 
```

Now on a different directory, clone Coffea repository
```
  git clone https://github.com/CoffeaTeam/CoffeaGrinder.git
  
 ```
### Run
Launch Jypyter Notebook
```
   cd CoffeaGrinder/
   jupyter notebook
```
It should open a new page in your default browser. If not , you can follow the link displayed in the terminal.  At the page you will find the directories and files from the location you ran jupyter notebook. select one of the samples to run it.

