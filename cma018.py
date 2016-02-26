#-------------------------------------------------------------------------------
# Author:      amir
# Created:     10/25/2015
#
# Instructions:
#
# 1) Make sure to rename the file (studentNetId.py) to your netId. (Do not include your first name, last name ... or any extra character)
# 2) To run the program type the following statement in the command line:  
#       -) python studentNetId.py DNASeq1FilePath DNASeq2FilePath OutputFilePath                                                                   
#    where  DNASeq1FilePath is the path to the file that contains First DNA sequence (e.g. DNASeq1.txt)
#           DNASeq2FilePath is the path to the file that contains Second DNA sequence (e.g. DNASeq2.txt)
#           OutputFilePath is the path that the output is goint to be saved (e.g. out.txt)
# 3) Make sure to replace FirstName, LastName, SectionNumber, NetId in studentInfo with your information
# 4) You may add as many functions as you want to this program
# 5) The core function in your program is DNASeqAlignment function, where it takes three arguments (DNASeq1,DNASeq2,outputPath) and 
#    computes the similarityScore, sequenceAlignment1 and sequenceAlignment2. At the end, the function writes the result to the output file (Do not make any changes to the output section).
# 6) sequenceAlignment1 and sequenceAlignment2 are strings and they are composed of following characters: 'A', 'T', 'G', 'C' and '-', Do not include extra space or any other character in the strings.
# 7) Make sure your program works with at least one of the following python versions: (2.7.9, 2.7.8, 2.6.6)
# 8) Once you have tested your program with one of the versions listed above, assign that version number to pythonVersion in studentInfo function
# 9) Make sure to write enough comments in order to make your code easy to understand. 
# 10) Describe your algorithm in ALGORITHM section below (you may add as many lines as you want).
# 11) To understand the flow of the program consider the following example:
#      0) Let say we have DNASeq1.txt file which contains AACCTGACATCTT and DNASeq2.txt file contains CCAGCGTCAACTT
#      1) If we execute the following command in the command line: -) python studentNetId.py DNASeq1.txt DNASeq2.txt out.txt
#      2) input arguments are parsed       
#      3) studentInfo() function will be executed and the output will be saved in out.txt file
#      4) DNASeqAlignment() function will be called
#      5) At the entry of the DNASeqAlignment function, DNASeq1='AACCTGACATCTT' and DNASeq2='CCAGCGTCAACTT'
#      6) You should compute the sequence alignment of DNASeq1 and DNASeq2. Let say the result is as follows:
#       A A C C T G A C - - - - A T C T T
#       | | | | | | | | | | | | | | | | |
#       - - C C A G - C G T C A A - C T T      
#      7) At the end of the DNASeqAlignment function sequenceAlignment1='AACCTGAC----ATCTT', sequenceAlignment2='--CCAG-CGTCAA-CTT', similarityScore=6.25
#      8) In the output section the result is going to be saved in out.txt file
#-------------------------------------------------------------------------------

# ALGORITHM: 
#
#
#
#
#


import os
import sys
import argparse

def studentInfo():
    pythonVersion = '2.7.6'
    
    student1FirstName = ""
    student1LastName = ""
    student1SectionNumber = ""
    student1NetId = ""
    
    student2FirstName = "FirstName2"
    student2LastName = "LastName2"
    student2SectionNumber = "SectionNumber2"
    student2NetId = "NetId2"
    
    info = 'Python version: ' + pythonVersion + '\n' + '\n'
    info = info + 'FirstName: ' + student1FirstName + '\n'
    info = info + 'LastName: ' + student1LastName + '\n'
    info = info + 'Section: ' + student1SectionNumber + '\n'
    info = info + 'NetId: ' + student1NetId + '\n' + '\n'
    
    info = info + 'FirstName: ' + student2FirstName + '\n'
    info = info + 'LastName: ' + student2LastName + '\n'
    info = info + 'Section: ' + student2SectionNumber + '\n'
    info = info + 'NetId: ' + student2NetId + '\n' + '\n'
    
    return info

def DNASeqAlignment(DNASeq1,DNASeq2,outputPath):
    similarityScore = -1
    sequenceAlignment1 = ''
    sequenceAlignment2 = ''
    
    #########################################################################################
    # Compute new values for similarityScore and sequenceAlignment1 and sequenceAlignment2  #                                                                  #
    #########################################################################################
    # the similarity scores can be computed as  #
    # S(0,j) = -0.2j #
    # S(i,0) = -0.2i #
    # the recurrence relation when i and j are positive is #
    # S(i,j) = max{S(i - 1,j - i)+ (i, j); S(i,j - 1) - 0.2; S(i - 1,j) - 0.2} #
    
    # the similirity scores are stored in a 2d matrix. 
    # each row i represents an element from DNA2 (y)
    # each column j represents an element from DNA1 (x)
    
    #########################################################################################
    # create a dictionary of dictionary to store the scores from delta table #
    
    #########################################################################################
    # Create two lists, one contains the similarity scores, the othere contains how to 
    # the index the scores are computed from 
    
    # in similarity score list:
    # Since 1d matrix is used to represent 2d matrix,
    # the element in (yj,xi) = yj * n + xi, where n is the size of sequenece yj
    
    # in direction list, where the index of the previous element is stored:
    # 'down' means S(i - 1,j) - 0.2 
    # 'slant' means max{S( i - 1,j - i)+ (i, j) 
    # 'right' means S(i,j - 1) - 0.2 
    # 'none' means it is the first element in the 2d matrix #
    #########################################################################################
 
    n = len(DNASeq1) + 1
    m = len(DNASeq2) + 1
    
    A = dict([('G', -0.1), ('C', -0.1), ('T', -0.15), ('A', 1.0)])
    G = dict([('A', -0.1), ('C', - 0.15), ('T', -0.1), ('G', 1.0)])
    C = dict([('T', -0.1), ('A', -0.1), ('G', - 0.15), ('C', 1.0)])
    T = dict([('A', -0.15), ('C', -0.1), ('G', -0.1), ('T', 1.0)])

    delta = dict([('A', A), ('G', G), ('C', C), ('T', T)])

    score = []
    direction = []
        
    j = 0
    i = 1
    
    for y in range(m):                # y represents columns ---- DNA 2#
        for x in range(n):            # x represents rows ---- DNA 1#
            
            if y == 0:
                score.append(-0.2 * j)
                if x == 0:
                    direction.append('none')
                else:
                    direction.append(j-1)
               
                j = j + 1
               # print y, x
                continue
                
            if x == 0:
                if y == 0:
                    continue
                
                else:
                    score.append(-0.2 * i)
                    direction.append((i-1)*n)
                    i = i + 1

#                    print y, x, score
                    
                continue
            
            score1 = score[(y - 1) * n + x] - 0.2 # down #
            score2 = score[(y - 1) * n + x - 1] + delta[DNASeq2[y-1]][DNASeq1[x-1]] # slant #
            score3 = score[y * n + x - 1] - 0.2 # right #
            
            tem = max(score1, score2, score3)
            score.append(tem)
            
            if tem == score1:
                direction.append(( y - 1) * n + x)
                #direction.append(('-', DNASeq2[y-1]))
                
            elif tem == score2:
                direction.append((y - 1) * n + x - 1)
                #direction.append((DNASeq1[x-1], DNASeq2[y-1]))
                
            elif tem == score3:
                direction.append(y * n + x - 1)
                #direction.append((DNASeq1[x-1], '-'))
                
            
           # print x, y, n
            #print score1, score2, score3, 'max', tem
            #print sequenceAlignment1, sequenceAlignment2
            #print (y-1)*n + x, (y-1)*n + x -1 , y*n + x - 1
            #print '\n'
            
            
    similarityScore = score[len(score)-1]
    
#    print direction
    
#    print score
    
    index = len(direction) - 1
    xi = index % n
    yj = index // n
    #print index, xi, yj
    #print DNASeq1, DNASeq2,'\n'
    
    str1 = ''
    str2 = ''
    while 1:
        if(direction[index] == (index-1)): # right
            str1 += DNASeq1[xi-1]
            str2 += '-'
            #print 'right',DNASeq1[xi-1],'-'
            #print str1, str2
            
        elif (direction[index] == (index - n)): # down
            str1 += '-'
            str2 += DNASeq2[yj-1]
            #print 'down','-',DNASeq2[yj-1]
            #print str1, str2
            
            
        else:
            str1 += DNASeq1[xi - 1]
            str2 += DNASeq2[yj - 1]
            #print 'slant',DNASeq1[xi - 1],DNASeq2[yj-1]
            #print str1, str2
        
        index = direction[index]
        xi = index % n
        yj = index // n
        #print index, xi, yj
        if(xi == 0):
            if(yj == 0):
                break
        
    
    
    
    
    #s1 = str1[::-1]
    #s2 = str2[::-1]
    #print s1, s2
    
    sequenceAlignment1 = str1[::-1]
    sequenceAlignment2 = str2[::-1]
    #print str1,str2
    
    
    
    #################################  Output Section  ######################################
    result = "Similarity score: " + str(similarityScore) + '\n'
    result = result + "Sequence alignment1: " + sequenceAlignment1 + '\n'
    result = result + "Sequence alignment2: " + sequenceAlignment2 + '\n'
    
    writeToFile(outputPath,result)
    
def writeToFile(filePath, content):
    with open(filePath,'a') as file:
        file.writelines(content)

def readFile(filePath):
    logLines = ''
    with open(filePath,'r') as file:
        for logText in file:
            logLines = logLines + logText

    uniqueChars = ''.join(set(logLines))
    for ch in uniqueChars:
        if ch not in ('A','a','C','c','G','g','T','t'):
            logLines = logLines.replace(ch,'')
    logLines = logLines.upper()
    return logLines

def removeFile(filePath):
    if os.path.isfile(filePath):
        os.remove(filePath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA sequence alignment')
    parser.add_argument('DNASeq1FilePath', type=str, help='Path to the file that contains First DNA sequence')
    parser.add_argument('DNASeq2FilePath', type=str, help='Path to the file that contains Second DNA sequence')
    parser.add_argument('OutputFilePath', type=str, help='Path to the output file')
    args = parser.parse_args()
    DNASeq1 = readFile(args.DNASeq1FilePath)
    DNASeq2 = readFile(args.DNASeq2FilePath)
    outputPath = args.OutputFilePath
    removeFile(outputPath)
    writeToFile(outputPath,studentInfo())
    DNASeqAlignment(DNASeq1,DNASeq2,outputPath)