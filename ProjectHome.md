# Introduction #
The G4Analysis is a program package write in python mostly base on MDAnalysis
Pacakge and 3DNA package. The aim of this package is to analyze the G-quartets structural parameters, mostly rise and twist from the molecular dynamics trajectories. MD softwares like Gromacs, Amber, and NAMD trajectories all can be read in this program.

# Usage #
**Files**

|Option |   Type   |   Filename  |   Description                |
|:------|:---------|:------------|:-----------------------------|
|-p     |    Input |  coor\_file  |  Structure fie: gro pdb etc. |
|-f     |   Input  |  traj\_file  |  Trajectory: xtc trr.        |
|-o     |   Input  | output\_file | xvgr/xmgr file.              |
|-i     |   Input  | para\_an.in  | input parmarter file.        |


**Other options**


|Option |  Type  |  Value   |     Description                              |
|:------|:-------|:---------|:---------------------------------------------|
|--rise |   bool |   False  |      Calculate the distance of DNA bases groups. |
|--twist|   bool |   False  |      Calculate the twist of DNA bases groups.    |
|--rmsd |   bool |   False  |      skip Calculate the RMSD of DNA bases groups.|
|--begin|   int  |   0      |      First frame (ps) to read from trajectory.   |
|--end  |   int  |   -1     |      Last frame (ps) to read from trajectory.    |
|--skip |   int  |   1      |      Get frames when frame MOD skip = 0          |
|-h     |   bool |   yes    |      Print help info and quit                    |

# Interactive Mode #
G4Analysis -p coor\_file -f traj\_file -o output\_file –rise/twist/rmsd [–begin/end/skip]

# Input File Mode #
G4Analysis -p coor\_file -f traj\_file -i input\_file [–begin/end/skip]