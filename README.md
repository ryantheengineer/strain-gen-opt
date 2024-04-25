# strain-gen-opt

Minimization of circuit board strain using MATLAB-based finite-element analysis
package and a genetic algorithm optimization.

## Installation

Clone the git repository and install in your local directory. Requires the
MATLAB Runtime and a proprietary `FEA.exe` to run.

## Usage

### First Run
For convenience, before running an optimization, open the script runFEA.py and
uncomment the following lines, under the main code:

```python
initialdirFEA = "E:/github/strain-gen-opt/FEA"
FEApath = chooseFEApath(initialdirFEA)
```

Change the directory path for `initialdirFEA` to a directory on your machine,
then run `runFEA.py` in the terminal. A GUI will open where you can select your
local installation of `FEA.exe`. Complete the code run, which will create a 
file called `FEApath.pk` in your local repository. Comment out the previously
mentioned lines and `runFEA.py` will remember the location of the FEA
executable for all future FEA cases, which will make the optimization run
smoothly in the background.

### Main Run
To perform a design optimization, run `optimize.py` in the terminal. This 
script will start by opening a file selection GUI where you can choose your 
input design. Typically this will be a file named Fea.xml.

