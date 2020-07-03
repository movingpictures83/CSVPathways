# CSVPathways
# Language: Python
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, Python 3.6
# Dependency: PythonCyc 1.0

PluMA plugin to take a multiomics CSV file and query the PathwayTools
(Karp et al, 2015) database to retrieve a list of all metabolic pathways
involving taxa and/or metabolites in the CSV. 

The plugin relies on PathwayTools being installed and running on a host.

The input is a parameters TXT file of tab-delimited keyword-value pairs:
hostname: The hostname on which PathwayTools is currently running
csvfile: Multiomics CSV file (typically abundances)
mapping: Maps metabolite names in the CSV to PathwayTools identifiers

Using the user-specified PREFIX, two output files will be generated:
prefix.txt: List of pathways involving taxa/metabolites in the CSV
prefix.noa: Table of CSV taxa/metabolites and their PathwayTools reaction identifiers
