# ExportMDF
Converts Matlab/Simulink data saved in workspace or 'MAT' file to DAT format (MDF 3.0) which can be read by ETAS MDA(INCA) or CANape tool.
Currently, script works using Matlab libraries version 8.2 (R2013) or higher. 

# Usage
For use with Matlab App 
-
1) The script uses simulation data saved to MAT file.  First use either of the following ways to save simulation data to MAT file - 

 - Use 'to workspace' Simulink block for data logging respective signals (use option "Save format" as "Timeseries") & in the workspace select the variable names generated after simulation to create the MAT file.  
 - Use 'ToFile' Simulink block with 'bus' block to write multiple signals to same *.mat file (use option "Save format" as "Timeseries") 
 - Export simulation data using signal logging (check link: https://www.mathworks.com/help/simulink/ug/exporting-signal-data-using-signal-logging.html).  Use the simulation output variable of signal logging (eg. logsout) stored as "Simulink.SimulationData.Dataset" datatype and create a MAT file.    

2) Install the "ExportMDF" app and use the app to select the MAT file. 

3) Output file is stored in the same location as input mat file.

For use with inline simulation 
- 
First export simulation data using signal logging (check link: https://www.mathworks.com/help/simulink/ug/exporting-signal-data-using-signal-logging.html).  

Then use the simulation output variable of signal logging (eg. logsout) stored as "Simulink.SimulationData.Dataset" datatype with the script.  Eg.- mat2dat(logsout, 'target_filename.dat')

For inline, script uses two parameters; 
 - param1 - exported simulation data stored as "Simulink.SimulationData.Dataset" datatype
 - param2 - target filename
