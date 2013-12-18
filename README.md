MarkovFlowMatrices
==================

This is a Python script for use in the field of land change science. It performs Markov and Flow calculations on a given transition matrix.

MarkovAndFlowScript.py Read-me File
MarkovAndFlowScript.py is  a python script file which takes as an input a “TransitionMatrix.csv” file containing a land cover transition matrix between two time points (this file must be formatted correctly – see the TransitionMatrix.csv read me file.
The script file outputs a “MatrixOutput.csv” file with information related to Markov and Flow matrices for this land cover change.
IMPORTANT!!
1) You must change the folderpath variable near the top of the script to point to the file directory containing  your python script and the TransitionMatrix.csv file. 
Other important points:
1) keep the TransitionMatrix.csv file in the same directory as your python script. 
2) if necessary, uncomment and change the folderpath variable at the top of the script to the file path where your “TransitionMatrix.csv” file is located. 
3) you must have certain Python libraries installed for this script to work, notably sys, string, csv, and numpy. If necessary, uncomment the “sys.path.append…” line at the top of the script to point to where your downloaded libraries are.
Numpy library available here: http://www.numpy.org/
4) Youtube video available here: 

