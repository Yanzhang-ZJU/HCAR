01.Install   
   environment configuration
    * numpy 1.25.1
    * mdtraj 1.9.8
    * python 3.11.4
    * gromacswrapper 0.8.4
    * biopandas 0.4.1
    * tqdm 4.65.0
    * typing-extension 4.7.1

02.Input&Output
   Input of comlig.py: The trajectory file, template file and index file of MD simulation system. These three files need to be specified, details will be found in comments of comlig.py.
   The default output of comlig.py is comlig.pdb. It contains the center of mass coordinates of the selected atoms in the trajectory.
   Output of comlig.py is used as input of grid.py. Meanwhile the step size for the grid also need to be specified.
   Output of grid.py is new_comlig.pdb. It contains the frequency distribution information of the centroid.

03.Usage
   For comlig.py: 
       set three parameters in comlig.py include, xtcName, topName, comIndName
       python comlig.py
   For grid.py: 
       python grid.py <input pdb file> <step size>

04.Example
   a. comlig.py
   All input will be found in example/, include, traj.xtc, temp.gro, index.ndx. The three parameters required by comlig.py have also been set in line 44-46.
   python comlig.py
   b. grid.py
   python grid.py comlig.pdb 1
   In this example, <step size> is set to 1
   All files have been checked, and output files can be found in output/.
