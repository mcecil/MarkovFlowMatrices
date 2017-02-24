MarkovFlowMatrices
==================

This is a Python script for use in the field of land change science. It performs Markov and Flow calculations on a given transition matrix.

MarkovAndFlowScript.py is  a python script file which takes as an input a “TransitionMatrix.csv” file containing a land cover transition matrix between two time points (this file must be formatted correctly – see the TransitionMatrix.csv read me file.)
The script file outputs a “MatrixOutput.csv” file with information related to Markov and Flow matrices for this land cover change.


IMPORTANT!!
1) You must change the folderpath variable (in line 16) near the top of the script to point to the file directory containing  your python script and the TransitionMatrix.csv file. 

Other important points:

2) Keep the TransitionMatrix.csv file in the same directory as your python script. 

3) If necessary, uncomment and change the folderpath variable at the top of the script to the file path where your “TransitionMatrix.csv” file is located. 

4) You must have certain Python libraries installed for this script to work, notably sys, string, csv, and numpy. If necessary, uncomment the “sys.path.append…” line at the top of the script to point to where your downloaded libraries are.
Numpy library available here: http://www.numpy.org/

5) Youtube video available here: http://www.youtube.com/watch?v=UD9ZTW7rFE0

6) Download link for related files: https://drive.google.com/file/d/0B1Dh-GK4po1rempzbzVUUGlCRms/edit?usp=sharing 
This includes a chart template, other sample input transition matrix files, and read-me files for the script itself, the input TransitionMatrix.csv file, and the output MatrixOutput.csv file. 
