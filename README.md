# SEAPy_Climate
These codes perform diagnostics of mean precip, winds, and SSTs over southeast Asia
It outputs plots and netcdf files of NDJF mean, and composites of different cold surge types.

Designed to run in Python version 3.+

Dependencies to install: 
Iris (https://scitools.org.uk/iris/docs/latest/index.html)
Matplotlib (https://matplotlib.org/)
pandas (https://pandas.pydata.org/)


Run the code:
1. Edit the main.py with information on the model runs/observations
2. Enable/disable sections/subsections in main.py   
3. Run 'python main.py'
   
   This will generate a number of plots and netcdf files that can be
   used for model evaluation
   
A basic workflow of the code is below:
![SEAPy_flowchart](https://user-images.githubusercontent.com/2921691/158994431-9ea0d0d1-e66c-46e7-8947-b5939d1718b5.png)
