"""
Created on Apr 14, 2014

@author: Jorge Munoz, Gabriel Alsheikh, and James Gomez
"""

import random
import sys
import re

nucleo_base_list = ["C", "T", "G", "A"]


def remove_from_string(arr, index):
    arr1 = arr[0:index]
    arr2 = arr[index + 1:]
    return arr1 + arr2


def remove_range_from_string(arr, index, size):
    arr1 = arr[0:index]
    arr2 = arr[index + size:]
    return arr1 + arr2


def insert_to_string(arr, index, char):
    arr1 = arr[0:index]
    arr2 = arr[index:]
    return arr1 + char + arr2


def copy_to_string(arr, index, size):
    arr1 = arr[index[0]:index[0] + size]
    #copy the substring into the array
    for i in range(1, len(index)):
        arr = remove_range_from_string(arr, index[i], len(arr1))
        arr = insert_to_string(arr, index[i], arr1)
    return arr


def generate_STR(length):    
    STR=[]
    for i in range(int(length*0.0001)):
        seqLen=random.randint(2,5)
        copies=random.randint(10,20)        
        temp=''
        for i in range(seqLen):
            temp += random.choice(nucleo_base_list)        
        STR.append([temp,temp*copies])
    return STR


def generate_ref_genome(genome_id, num_chromosomes, length_chromosome):
    """
    Generates a random reference genome with the specified number of chromosomes,
    each of length length_chromosome
    """
    print "Generating reference genome..."
    ref_file = open("ref_" + genome_id + ".txt", "w")
    ref_file.write(">" + str(genome_id))

    genome = []
    STR=generate_STR(length_chromosome)
        
    STRpos = []
    genome = ''
    #
    #    
    #Generate the string, then write it
    for i in range(1, num_chromosomes + 1):
        ref_file.write("\n>chr" + str(i) + "\n")
        # Generate the string
        for j in range(0, length_chromosome):
            genome += random.choice(nucleo_base_list)
        # Add STRs to the genome
        # for j in range(len(STR)):
        #     tmp=random.randint(0, length_chromosome - len(STR[j]))
        #
        #     genome = remove_range_from_string(genome, tmp, len(STR[j][1]))
        #     genome = insert_to_string(genome, tmp, str(STR[j][1]))
        #     STRpos.append([tmp, STR[j][0], STR[j][1]])

        for j in range(0, length_chromosome):
            # write a maximum of 80 alleles per line
            if j != 0 and j % 80 == 0:
                ref_file.write("\n")
            ref_file.write(genome[j])             
    
            
    print "Reference genome complete"
    ref_file.close()

    return (ref_file, STRpos)


def invert_str(arr, index, size):
    arr1 = arr[index:index + size]
    arr = remove_range_from_string(arr, index, size)
    arr = insert_to_string(arr, index, arr1[::-1])
    return arr, arr1


def usage(name):
    print "USAGE: python " + str(
        name) + " <genome_id> <#chromosomes> <chromosome_size x 1K> <coverage>"


def modSTR(genome, STR):
    fullSTR = []
    INS = []
    DEL = []

    for j in range(len(STR)):
        tmp=random.randint(-2,2)
        tmpStr='' 

        if(tmp>0):
            genome = remove_range_from_string(genome, int(STR[j][0])+len(STR[j][2]), tmp*len(STR[j][1]))
            genome = insert_to_string(genome, int(STR[j][0]), tmp*str(STR[j][1]))
            
            DEL.append([tmp*str(STR[j][1]), str(int(STR[j][0])+len(STR[j][2])-tmp*len(STR[j][1]))])
            
            fullSTR.append([str(STR[j][0]), str(STR[j][1]), str(str(STR[j][2])+tmp*str(STR[j][1]))])     
        if(tmp<0):
            for k in range(len(STR[j][1])):
                tmpStr+=random.choice(nucleo_base_list)    

            genome = remove_range_from_string(genome, int(STR[j][0]), -tmp * len(STR[j][1]))
            genome = insert_to_string(genome, int(STR[j][0]), -tmp*tmpStr)
            
            INS.append([tmpStr,str(STR[j][0])])
            
            fullSTR.append([str(STR[j][0]), str(STR[j][1]), str(STR[j][2][-tmp*len(STR[j][1]):])])
            
        if(tmp == 0):
            fullSTR.append([str(STR[j][0]), str(STR[j][1]), str(STR[j][2])])
            
    fullSTR = sorted(fullSTR, key=lambda fullSTR: int(fullSTR[:][0])) 
    
    return genome, INS, DEL, fullSTR
        
        
################################# START OF SCRIPT ###################################
if len(sys.argv) < 5:
    usage(sys.argv[0])
    sys.exit()

# get parameters from console
genome_id = str(sys.argv[1])
num_chromosomes = int(sys.argv[2])
# chromosome_size = int(sys.argv[3]) * 1000000
chromosome_size = int(sys.argv[3]) * 1000
coverage = float(sys.argv[4])

print "\ngenome-ID:\t" + genome_id
print "num-chroms:\t" + str(num_chromosomes)
print "chrom-size:\t" + str(chromosome_size)
print "coverage:\t" + str(coverage)

fullCOPYseq = []
fullCOPY = []
fullINVseq = []
fullINV = []
fullINS = []
fullDEL = []
fullSNP = []
fullSTR = []

# create unaltered reference genome
baseFile, STR = generate_ref_genome(genome_id, num_chromosomes, chromosome_size)

private = True
baseFile = open("ref_" + genome_id + ".txt", "r")
if private:
    baseFile2 = open("private_" + genome_id + ".txt", "w")
    baseFile2.write(">" + genome_id + "\n")

readsFile = open("reads_" + genome_id + ".txt", "w")
readsFile.write(">" + genome_id + "\n")
# skip the first two '>' labels in the ref genome file
baseFile.readline()
baseFile.readline()
for chromosome in range(1, num_chromosomes + 1):
    counter = 0
    print "Generating mutations in chromosome: " + str(chromosome)

    # Lists of indices at which deletions and inserts were made
    insList = []
    delList = []
    strList = []
    # get reference genome as single string
    baseFileList = ""
    for line in baseFile:
        if ">" in line:
            break
        baseFileList += str(line)
    
    baseFileList=baseFileList.replace("\n", "")
    
    #Modify the STRs
    # baseFileList, insList, delList, strList=modSTR(baseFileList,STR)
    # fullSTR.append(strList)

    #Copy Numbers
    #Sequence of random length between 20-50
    copyIndex = []
    copyLen = []
    # print "\tGenerating copies..."
    # for copyLoop in range(0, 1):
    #     copyLen.append(random.randint(20, 50))
    #     copyIndex.append(
    #         int(random.random() * (chromosome_size - copyLen[copyLoop])))
    #
    #     #0.001% of the time the string is copied
    #     for i in range(0, int(chromosome_size * 0.00001)):
    #         temp = (int(random.random() * (chromosome_size - copyLen[copyLoop])))
    #         write = 1
    #         for j in range(0, len(copyIndex) - 1):
    #             if copyIndex[j] < temp < copyIndex[j] + copyLen[copyLoop]:
    #                 write = 0
    #                 break
    #         if write == 1:
    #             copyLen.append(copyLen[copyLoop])
    #             copyIndex.append(temp)
    #     fullCOPY.append(copyIndex)
    #     fullCOPYseq.append(
    #         str(baseFileList[copyIndex[0]:copyIndex[0] + copyLen[copyLoop]]))
    #     baseFileList = copy_to_string(baseFileList, copyIndex, copyLen[copyLoop])

    #Inversions
    #0.001% of string is inverted
    actualInv = 0
    invIndex = []
    orig = []
    # print "\tGenerating inversions..."
    # for invLoop in range(0, int(0.00001 * chromosome_size)):
    #     write = 1
    #     temp = int(random.random() * (chromosome_size - copyLen[copyLoop]))
    #     invLen = random.randint(20, 50)
    #
    #     for i in range(0, len(copyIndex)):
    #         if ((temp > copyIndex[i] and temp < copyIndex[i] + copyLen[i]) or (temp +invLen > copyIndex[i] and temp+invLen < copyIndex[i] + copyLen[i])):
    #             write = 0
    #             invLoop-=1
    #             break
    #     if write == 1:
    #         actualInv += 1
    #         # baseFileList=
    #         # orig =
    #
    #         baseFileList, orig= invert_str(baseFileList, temp, invLen)
    #         invIndex.append([orig, temp])
    #
    #
    # # Sort invIndex by the index of the inversion
    # invIndex = sorted(invIndex, key=lambda invIndex: invIndex[:][1])
    #
    # fullINV.append(invIndex)
    
    #INS/DELS



    # Insertions/Deletions, split into sections of 2,000 (0.5% ins/del)
    # 500 below comes from Sequence length / (Seq. length * 0.1% * 2) = 1 / (0.1% * 2)
    print "\tGenerating indels..."
    sectionLen = 200 # int(chromosome_size * 0.001 * 2)
    for i in range(0, int(chromosome_size / sectionLen)):

        # Make deletions in i-th section
        # Get index to delete from

        # Number of nucleotides to insert and delete
        numToInsDel = random.randint(1, 5)

        # Sequence of nucleotides deleted
        delSeq = []
        delStartIndex = random.randint(i * sectionLen, ((i + 1) * sectionLen) - numToInsDel)
        # Delete the nucleotides starting at the index
        for j in range(0, numToInsDel):     
            nucleoDeleted = baseFileList[delStartIndex]
            baseFileList = remove_from_string(baseFileList, delStartIndex)
            delSeq.append(nucleoDeleted)

        delList.append([str(''.join(delSeq)), str(delStartIndex)])
        #delList.append(''.join(delSeq) + "," + str(delStartIndex))

        # Make insertions in i-th section
        # Get index to insert at

        # Sequence of nucleotides inserted
        insSeq = []
        insStartIndex = random.randint(i * sectionLen,
                                       ((i + 1) * sectionLen) - numToInsDel)
        # Insert the nucleotides starting at the index
        for j in range(0, numToInsDel):
            randomNucleo = random.choice(nucleo_base_list)
            baseFileList = insert_to_string(baseFileList, insStartIndex + j,
                                            str(randomNucleo))
            insSeq.append(randomNucleo)

        insList.append([str(''.join(insSeq)), str(insStartIndex)])
        #insList.append(''.join(insSeq) + "," + str(insStartIndex))
    insList = sorted(insList, key=lambda insList: int(insList[:][1])) 
    delList = sorted(delList, key=lambda delList: int(delList[:][1])) 
    fullINS.append(insList)
    fullDEL.append(delList)

    # END INSERTS AND DELETIONS

    #SNPS
    """
        Generates snps at random indeces in the array and returns a list
        of the indeces of snps in sorted ascending order
    """

    count = 0.0
    snp_list = []
    tempList=baseFileList
    import difflib
    
    print "\tGenerating SNPs..."
    while (count / len(baseFileList)) < .003:
        index = random.randint(0, len(baseFileList))
        if baseFileList[index] == "A":
            tp="A"
            s = random.choice(["C", "G", "T"])
        elif baseFileList[index] == "C":
            tp="C"
            s = random.choice(["A", "G", "T"])
        elif baseFileList[index] == "G":
            tp="G"
            s = random.choice(["A", "C", "T"])
        elif baseFileList[index] == "T":
            tp="T"
            s = random.choice(["A", "C", "G"])

        baseFileList=remove_from_string(baseFileList, index)
        baseFileList=insert_to_string(baseFileList, index, s)
        snp_list.append([tp, s, index])
        count += 1

    # Sort SNPs by increasing index
    snp_list = sorted(snp_list, key=lambda snps: snps[:][2]) 
    snps = snp_list

    fullSNP.append(snps)
    if private:
        print "\tWriting to private genome..."
        baseFile2.write(">chr" + str(chromosome) + "\n")
        counter = 0
        write_newline = False
        for c in baseFileList:
            if write_newline:
                baseFile2.write("\n")
                write_newline = False
            baseFile2.write(str(c))
            counter += 1
            if counter % 80 is 0:
                write_newline = True
        
        baseFile2.write("\n")
    #READS
    # File to hold reads from 1 million char sequence
    print "\tGenerating reads..."
    readsFile.write(">chr" + str(chromosome) + "\n")
    coverage /= float(100)
    for i in range(0, int(chromosome_size * coverage)):
        # First read
        startIndexPart1 = (int(random.random() * (chromosome_size - 210)))
        randomGap = random.randint(90, 110)
        # Second read
        startIndexPart2 = startIndexPart1 + randomGap + 1
        readList = baseFileList[startIndexPart1:startIndexPart1 + 210]

        # 1% error in reads 10% reads are garbage
        # Get value between 0 and 1
        randomReadCondition = random.random()

        #Throw in an error in 1% of the read length 200
        # errors = int(0.01 * 200)
        # for i in range(0, errors):
        #     #pick an index
        #     randomIndex = random.randint(0, 99)
        #     if randomIndex < 50:
        #         pass
        #     else:
        #         randomIndex += randomGap - 50
        #
        #     randomReadError = random.choice(nucleo_base_list)
        #     while randomReadError == readList[randomIndex]:
        #         randomReadError = random.choice(nucleo_base_list)

        #Full Garbage Read
        # if 0.01 < randomReadCondition < 0.11:
        #     readList = ""
        #     for i in range(0, 50):
        #         readList += (random.choice(nucleo_base_list) )
        #     for i in range(0, randomGap):
        #         readList += '-'
        #     for i in range(50, 100):
        #         readList += random.choice(nucleo_base_list)
        
        #for i in range(0, 10):
        #    readList += random.choice(nucleo_base_list)    
            
            
        readsFile.write(str(readList[:50]))
        readsFile.write(',')
        readsFile.write(str(readList[50 + randomGap:randomGap + 100]))
        readsFile.write("\n")


# create answer key
print "Creating answer key..."
baseAnswerFile = open("ans_" + genome_id + ".txt", "w")
baseAnswerFile.write(">")
baseAnswerFile.write(str(genome_id))
baseAnswerFile.write("\n")

#copies
try:
    fullCOPYseq[0]
    baseAnswerFile.write(">COPY: \n")
    for i in range(0, num_chromosomes):
        baseAnswerFile.write(str(i + 1) + "," + str(fullCOPYseq[i]))
        for j in range(len(fullCOPY[i])):
            baseAnswerFile.write("," + str(fullCOPY[i][j]))
        baseAnswerFile.write("\n")
except IndexError:
    pass

#inversions
try:
    fullINV[0]
    baseAnswerFile.write(">INVERSION: \n")
    for i in range(0, num_chromosomes):
        for j in range(0, len(fullINV[i])):
            baseAnswerFile.write(str(i + 1) + "," + str(fullINV[i][j][0]) + "," + (str)(
                fullINV[i][j][1]) + "\n")
except IndexError:
    pass

#STR
try:
    fullSTR[0]
    STRFile=open("STRtest.txt", "w")
    baseAnswerFile.write(">STR: \n")
    for i in range(0, num_chromosomes):
        for j in range(0, len(fullSTR[i])):
            baseAnswerFile.write(str(i + 1) +","+str(
                        fullSTR[i][j][1])+ ","+str(len(fullSTR[i][j][2])/len(fullSTR[i][j][1])) + ","+str(fullSTR[i][j][0])+"\n")

            STRFile.write(str(i + 1) +"," + (str)(
                        fullSTR[i][j][0]) + ","+str(fullSTR[i][j][2])+ "\n")
    STRFile.close()
except IndexError:
    pass

#inserts
try:
    fullINS[0]
    baseAnswerFile.write(">INSERT:\n")
    for i in range(0, num_chromosomes):
        for j in range(0, len(fullINS[i])):
            baseAnswerFile.write(str(i + 1) + "," + str(fullINS[i][j][0]) + "," + str(fullINS[i][j][1]))
            baseAnswerFile.write("\n")
except IndexError:
    pass

#deletes
try:
    fullDEL[0]
    baseAnswerFile.write(">DELETE:\n")
    for i in range(0, num_chromosomes):
        for j in range(0, len(fullDEL[i])):
            baseAnswerFile.write(
                str(i + 1) + "," + str(fullDEL[i][j][0]) + "," + str(fullDEL[i][j][1]))
            baseAnswerFile.write("\n")
except IndexError:
    pass
        #snps

try:
    fullSNP[0]
    baseAnswerFile.write(">SNP")
    for i in range(0, num_chromosomes):
        for j in range(0, len(fullSNP[i])):
            baseAnswerFile.write("\n" + str(i + 1) + "," + str(fullSNP[i][j][0]) + ',')
            baseAnswerFile.write(str(fullSNP[i][j][1]) + ',')
            baseAnswerFile.write(str(fullSNP[i][j][2]))
except IndexError:
    pass

baseAnswerFile.close()
readsFile.close()
baseFile2.close()
print "DONE"
