# SEAPy_Climate
These codes perform diagnostics of mean precip, winds, and SSTs over southeast Asia
It outputs plots and netcdf files of NDJF mean, and composites of different cold surge types.

Designed to run in Python version 3.+

Dependencies to install: 
Iris (https://scitools.org.uk/iris/docs/latest/index.html)
Matplotlib (https://matplotlib.org/)
pandas (https://pandas.pydata.org/)


Run the combined SEA and MJO workflow (Python 3.10 or later):
1. Install the package and dependencies with `pip install -e .`.
2. Edit `sea_runs.yaml` with the control, experiment and observation settings.
3. Check the site-specific data/output paths in `main.py` and the retrieval settings.
4. Run `python main.py --help` for diagnostic selection and processing options.

```bash
python main.py --runs obs control experiment
python main.py --runs experiment --mjo-levels 1 2 --sea-levels none
python main.py --html-only
```
   
   This will generate a number of plots and netcdf files that can be
   used for model evaluation

## Saved-output assessment

Install the optional reporting dependencies with `pip install -e '.[assessment]'`.
The scripts in `analysis_gc6_mjo_sea/` compare existing diagnostic products without
rerunning the model diagnostics:

```bash
python analysis_gc6_mjo_sea/analyse_outputs.py
python analysis_gc6_mjo_sea/build_focused_report.py
```

Review the configured output root and run identifiers in `analyse_outputs.py`
before applying the scripts to other runs. The focused report creates a Word
assessment and two figures using surge composites, wave propagation, spectra
and fractional variance. Use `--figures-only` to rebuild figures and the Word
report from saved evidence, or `--document-only` to rebuild just the document.
Generated reports, evidence tables, plots and logs are not committed. These
assessments retain sampling, provenance and composite-generation caveats;
they do not establish statistical significance. LibreOffice can optionally
render the Word file to PDF for layout checking.
   
A basic workflow of the code is below:
![SEAPy_flowchart](https://user-images.githubusercontent.com/2921691/158994431-9ea0d0d1-e66c-46e7-8947-b5939d1718b5.png)
