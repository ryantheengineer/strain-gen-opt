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

### Typical Usage
To perform a design optimization, run `optimize.py` in the terminal. This 
script will start by opening a file selection GUI where you can choose your 
input design. Typically this will be a file named `Fea.xml`.

Once the input XML file is selected, the terminal will print some statements
about reading in design constraints. The script will make an educated guess 
about the number of pressure rods the design can fit on the top side and prompt
the user to change the value of `nprods_top` in the terminal, if desired, based
on this information.

Once all the design parameters are set, the initial design population will be
created based on the values set in the main function of `optimize.py`. The
progress of the creation of the initial designs, or chromosomes, will print
to the terminal.

Next, the genetic algorithm loop will take over. The current generation number
will be printed, and status statements will be printed for the creation of
child designs for each of the three child creation methods: crossover,
mutation, and local search (a special case of mutation). Then `FEA.exe` will be
called for each of the potential designs. The code uses some multiprocessing
libraries to speed up this part of the process.

After the FEA cases are computed, design fitness is calculated for each of the
objectives. The designs from the current generation are scored against each
other and the best designs are selected as parents for the next generation.
This continues until the maximum generations are reached, or if the appropriate
flags are set to `True`, then the optimization continues until at least one
design of the current generation has values for all objectives of less than 500
microstrain.

Status plots are produced with each generation to show the characteristics of
the best overall design for each generation. These can help the user determine
if convergence issues exist.

## Parameters and Tuning
The genetic algorithm parameters can be found in the main function of 
`optimize.py`.

### Useful Parameters
| Variable Name | Description |
| ---- | ---- |
| `pop_size` | The number of designs that will be selected as parents out of each generation. Also the size of the initial random population. |
| `rate_crossover` | The number of child designs that will be created for each generation, using the crossover method. |
| `rate_mutation` | The number of child designs that will be created for each generation, using the mutation method. |
| `chance_mutation` | The normalized percent chance that an individual pressure rod will be mutated within a given child design. |
| `n_searched` | The number of child designs that will be created for each generation, using the local search method. |
| `chance_localsearch` | The normalized percent chance that an individual pressure rod will be mutated within a given child design. |
| `maximum_generation` | The number of generations that the optimization will be allowed to run to, unless options are selected to end early. |
| `end_early` | A boolean flag that determines whether the optimization will automatically stop early if an acceptable design is discovered. |
| `nprods_top` | The number of pressure rods to include on the top side of the UUT. |
| `nstandoffs` | The number of standoffs, or board stops, to include on the bottom side of the UUT. |
| `all_on` | A boolean flag to determine whether the design space should be simplified to a single pressure rod type. |
| `rod_type` | A string variable with the name of the single pressure rod type to be used, if desired. |
