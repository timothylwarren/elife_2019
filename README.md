 elife_2019
Analysis and plotting scripts for behavioral experiments in Sullivan, Warren, and Doe 2019
Author of Repository:
Timothy L. Warren, timlwarren AT gmail 

Files in this Repository: All to be implemented in Python 2.7


(1) write_elife_repo_csv.py
	
	This was used to write raw data files (in original acquisition format) to csv format for Dryad repository.

(2) raw_analysis.py 

	This does initial analysis on raw data.

	USAGE 
	python raw_analysis.py raw_data_range.txt

	COMMENTS
	Loads distributed raw data files within this date range, does analysis, saves .pck (pickle) files for each
	data file.
	Raw data are archived in csv files a Dryad repository: doi:10.5061/dryad.45177sc
	However, this script is not designed to run on those csv files but on original experimental files.
	I plan to write a new version of this script to run on Dryad files. Also, I can share raw experimental files.

(3) raw_data_range.txt

	Example of file format used by raw_analysis.py

(4) processed_analysis.py

	This does further analysis on data already analyzed by raw_analysis.py

	USAGE
	python processed_analysis.py processed_data_range.txt

	COMMENTS
	Data processing of .pck files created by raw_analysis.py. Calculates summary data for different experimental types. Saves a .pck file. A version of that .pck file is in repository: summary_data.pck

(5) processed_date_range.txt
	
	Example of file format used by processed_analysis.py


(6) plot_data.py 

	Makes plots used in paper. Loads summary .pck file (eg. 'summary_data.pck'). Saves plots as .svg files.

	USAGE
	python plot_data.py
