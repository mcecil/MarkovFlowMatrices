# find the limiting Markov distribution for different land cover classes in a two time-step analysis

# IMPORT MODULES
import sys
# sys.path.append("F:\GIS and Land Change")

import numpy
from numpy import matrix
from numpy import linalg

import string
import csv

# FILE MANAGEMENT
# here are our folder and file paths and names
folderPath = "C:\Users\mcecil\Desktop\GIS and Land Change Science\FinalFiles\\"
newFileName = "Output"
csvName = "TransitionMatrix.csv"
newFilePathAndName = folderPath + newFileName
csvPathAndName = folderPath + csvName 
f1 = open(csvPathAndName,'r')

# define and open csv output file
ofile  = open('MatrixOutput.csv', "wb")
writer = csv.writer(ofile, delimiter=',',)

# FUNCTION DEFINITIONS (can be added to when code is refactored)
# define operation to convert matrices to lists, so that we can write them properly to csv output
def matrixToList(matrix):
    returnList = numpy.array(matrix).flatten().tolist()
    return returnList

# INITIAL CALCULATIONS
# get the number of categories
line1 = f1.readline()
line1Split = line1.split(",")
numberCategories = int(line1Split[0])

# get the start year, end year, and length of time interval
line2 = f1.readline()
line2Split = line2.split(",")
startYear = int(line2Split [0])
endYear = int (line2Split [1])
intervalLength = endYear - startYear

# create header list to write to output
headerList = ["Number of Categories", numberCategories, "Start Year", startYear, "End Year", endYear, "Time Interval", intervalLength]
writer.writerow(headerList)


# get the category names
line3 = f1.readline()
line3Split = line3.split(",")
categoryNames = line3Split[0:numberCategories]
categoryNames[numberCategories - 1] = categoryNames[numberCategories - 1].split("\n")[0]
writer.writerow(categoryNames)



# MAKE TRANSITION MATRIX

matrixList = []
for i in range(0, numberCategories):
    tempLine = f1.readline()
    tempList = tempLine.split(",")
    tempList[numberCategories - 1] = tempList[numberCategories - 1].split("\n")[0]
    matrixList.append(tempList)


# make a copy of the transitionMatrix for later reference
transitionMatrix = matrix(matrixList)
transitionMatrixStatic = numpy.zeros((numberCategories, numberCategories))
for i in range(0, numberCategories):
    for j in range(0, numberCategories):
        transitionMatrixStatic[i,j] = float(transitionMatrix[i,j]) # transitionMatrixStatic is the matrix for later reference


# MAKE MARKOV MATRIX


markovMatrix = numpy.zeros((numberCategories, numberCategories)) # this step initializes a float Matrix of the appropriate size
# make everything a float
for i in range(0, numberCategories):
    for j in range(0, numberCategories):
        markovMatrix[i,j] = float(transitionMatrix[i,j])

transitionMatrixFloat = markovMatrix.copy()
# divide by row sum
# i is row, j is column
for i in range(0, numberCategories):
    rowSum = 0
    for j in range(0, numberCategories):
        rowSum = rowSum + markovMatrix[i,j]
    for j in range(0, numberCategories):
        markovMatrix[i,j] = markovMatrix[i,j] / rowSum
markovMatrixStatic = markovMatrix.copy() # this is a copy of the markov matrix for future reference

# MARKOV EQUILIBRIUM CALCULATIONS (for more details, consult 


# transpose the Markov matrix 
squareMatrix = markovMatrix.T

# change all the last row entries to 1

for j in range(0, numberCategories):
    squareMatrix[numberCategories - 1, j] = 1

# subtract one from all diagonal entries except for last row

for i in range (0, numberCategories - 1):
    squareMatrix[i,i] = squareMatrix[i,i] - 1

# initialize the column matrix that is the other side of the equation
columnMatrix = numpy.zeros((numberCategories, 1))
# set last row entry in column matrix equal to 1
columnMatrix [numberCategories -1, 0] = 1

# solve the equation
Y = linalg.solve(squareMatrix, columnMatrix)


# make a copy of the equilibrium matrix and take its transpose
equilibriumMatrixT = Y.copy()
equilibriumMatrix = equilibriumMatrixT.T

# calculate initial land cover distribution matrix
landCover = numpy.zeros((1, numberCategories))
for j in range(0, numberCategories):
    columnSum = 0
    for i in range(0, numberCategories):
        columnSum = columnSum + float(transitionMatrix[i,j])
    landCover[0,j] = columnSum # this makes the landCover matrix the sum of all cells for that category at the final time point.


# these steps divide each entry in the landCover matrix by the total number of cells. The result is the decimal proportion each category represents of total land cover at the ending time point.
cellSum = 0
for j in range(0, numberCategories):
    cellSum = cellSum + landCover[0,j]
for j in range (0, numberCategories):
    landCover[0,j] = landCover [0,j] / cellSum

# write initial land cover distribution to output
writer.writerow(["Initial land cover distribution", endYear])
writer.writerow(matrixToList(landCover))
                
landCoverM = numpy.matrix(landCover)
landCoverMList = numpy.array(landCoverM).flatten().tolist()
markovMatrixM = numpy.matrix(markovMatrixStatic)



# write header for Markov calculations
writer.writerow(["Markov matrix calculations"])
writer.writerow(["Markov matrix "])


# write Markov matrix to output
for i in range(0, numberCategories):
    tempMarkovMatrixRow = numpy.zeros((1, numberCategories))
    tempMarkovRow = matrixToList(tempMarkovMatrixRow)
    for j in range(0, numberCategories):
        tempMarkovRow[j]=markovMatrixStatic[i,j]
    writer.writerow(tempMarkovRow)


# calculate how many intervals it takes to get all categories < 1% to equilibrium
# the while loop below is rather long, but basically it continually applies the Markov probabilities for change to calculate new land cover distributions.
# This process stops when each category is < 0.01 from its equilibrium value, where '0.01' represents proportion of total land area


writer.writerow(["Markov Steps"])
writer.writerow(categoryNames) # writes category names so this copies nicely to Excel


tempMatrix = landCoverM.copy()
count = 0
check = 0
finalCount = -1
while count <= 1000:
    check = 0
    
    # check if equilbrium and initial land cover matrices are already within 0.1 of each other
    if count == 0:
        maxValue = 0
        differenceMatrix = numpy.zeros((1, numberCategories))
        for j in range(0, numberCategories):
            differenceMatrix [0,j] = abs(tempMatrix[0,j] - equilibriumMatrix[0,j])
        maxValue = max(differenceMatrix[0,j] for j in range(0, numberCategories))
        if maxValue < 0.01:
            finalCount = count
            count = 999999
            break
    count = count + 1
    tempMatrix = tempMatrix * markovMatrixM
    tempMatrix2 = numpy.zeros((1, numberCategories + 2))
    for j in range(0, numberCategories):
        tempMatrix2[0,j] = tempMatrix[0,j]
    templist = matrixToList(tempMatrix2)
    templist[numberCategories] = str(endYear + intervalLength * count)
    templist[numberCategories + 1] = "Markov Step: " + str(count)
    
    if count <= 10:
        writer.writerow(matrixToList(templist)) # the first 10 Markov steps are written to the csv file
    while check == 0:
        if count == 1000: # this stops the process after 1000 Markov steps
            check = 1
            writer.writerow(["No equilibrium found"])
        else:
            maxValue = 0
            differenceMatrix = numpy.zeros((1, numberCategories))
            for j in range(0, numberCategories):
                differenceMatrix [0,j] = abs(tempMatrix[0,j] - equilibriumMatrix[0,j])
            maxValue = max(differenceMatrix[0,j] for j in range(0, numberCategories))
            if maxValue < 0.01:
                check = 1
                finalCount = count
                if count > 10: # this prints to output the step where each category is less than 1 % from its equilibrium
                    templist[numberCategories + 1] = "Markov Step: " + str(count) + " (within 1% of equilibrium for each category)"
                    writer.writerow(matrixToList(templist)) # writes this step to csv
                count = 999999
            check = 1



# this list contains the equilibrium proportions
templist2 = list(templist) # makes a copy
for j in range(0, numberCategories):
    templist2[j] = matrixToList(Y)[j]
templist2[numberCategories] = "Equilbrium"
templist2[numberCategories + 1] = "Equilibrium" 
writer.writerow(templist2) # write equilibrium proportions to output


    
#FLOW MATRIX CALCULATIONS

# calculate Flow Matrix
writer.writerow(["Flow Matrix Calculations"])
# create copy of transition matrix
flowMatrix = transitionMatrixFloat.copy()
# convert from array type to matrix type
flowMatrixM = matrix(flowMatrix)


# calculate the sum of all transitions
flowSum = 0
for i in range(0, numberCategories):
    for j in range(0, numberCategories):
        flowSum = flowSum + flowMatrixM[i,j]
# divide by sum
for i in range(0, numberCategories):
    for j in range(0, numberCategories):
        flowMatrixM[i,j] = flowMatrixM[i,j] / flowSum
# make diagonal entries 0
for i in range (0, numberCategories):
    flowMatrixM[i,i] = 0


# annualize
for i in range(0, numberCategories):
    for j in range(0, numberCategories):
        flowMatrixM[i,j] = flowMatrixM[i,j] / intervalLength # this is the flow matrix



# write the flow matrix to output
writer.writerow(["Flow matrix "])
for i in range(0, numberCategories):
    tempFlowMatrixRow = numpy.zeros((1, numberCategories))
    tempFlowRow = matrixToList(tempFlowMatrixRow)
    for j in range(0, numberCategories):
        tempFlowRow[j]=flowMatrixM[i,j]
    writer.writerow(tempFlowRow)





# make a copy of land cover matrix
landCoverM1=landCoverM.copy()

# make net gains matrix
netGains = numpy.zeros((1, numberCategories))
for j in range(0, numberCategories):
    columnSum = 0
    for i in range(0, numberCategories):
        columnSum = columnSum + flowMatrixM[i,j]
    netGains[0,j] = columnSum



# make net losses matrix

netLosses = numpy.zeros((1, numberCategories))
for i in range(0, numberCategories):
    rowSum = 0
    for j in range(0, numberCategories):
        rowSum = rowSum + flowMatrixM[i,j]
    netLosses[0,i] = rowSum



# make net change matrix

netChange = numpy.zeros((1, numberCategories)) # netChange means net annual change
netChange = netGains - netLosses



# calculate extinction time.
# extinction time is initial percentage divided by negative net change entries (ignore positive ones)

extinctionYears = numpy.zeros((1, numberCategories))

for i in range(0, numberCategories):
    if netChange[0,i] < 0:
        extinctionYears[0,i] = abs((landCoverM1[0,i] / netChange [0,i]))

# calculate 1st extinction year
minExtinctionYears = 999999
for i in range(0, numberCategories):
    if extinctionYears[0,i] > 0 and extinctionYears[0,i] < minExtinctionYears:
        minExtinctionYears = extinctionYears[0,i]
minExtinctionDate = endYear + minExtinctionYears



# calculate net change until extinction date
totalNetChange = minExtinctionYears * netChange

# calculate proportions at extinction date
finalProportions = landCoverM + totalNetChange


writer.writerow(["Annual net change under flow matrix"])
writer.writerow(matrixToList(netChange))


# flow step csv output
writer.writerow(["Flow Steps"])
writer.writerow(categoryNames) # writes category names so this copies nicely to Excel

currentYear = endYear
flowStep = 0
flowProportions = landCover.copy()
while currentYear < minExtinctionDate:
    currentYear = endYear + flowStep * intervalLength
    for j in range(0, numberCategories):
        flowProportions[0,j] = landCover[0,j] + ((netChange[0,j]) * (currentYear - endYear))
    tempMatrix3 = numpy.zeros((1, numberCategories + 2))
    for j in range(0, numberCategories):
        tempMatrix3[0,j] = flowProportions[0,j]
    templist3 = matrixToList(tempMatrix3)
    templist3[numberCategories] = str(currentYear)
    templist3[numberCategories + 1] = "Flow Step: " + str(flowStep)
    
    writer.writerow(templist3)
    
    flowStep = flowStep + 1
    currentYear = currentYear + intervalLength
    

# this list will contain the final flow proportions
templist4 = list(templist3)
for j in range(0, numberCategories):
    templist4[j] = matrixToList(finalProportions)[j]
    templist4[numberCategories] = "End (" + str(int(minExtinctionDate)) + ")" # converts minExtinctionDate to int for the csv output
    templist4[numberCategories + 1] = "End"
    
writer.writerow(templist4)
writer.writerow(["Final Proportions under flow matrix"])
writer.writerow(matrixToList(finalProportions))
writer.writerow(["Year flow matrix ends: ", minExtinctionDate])


ofile.close()
f1.close()













