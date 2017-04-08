# -- coding: utf-8 -*-

from collections import defaultdict
import os
import os.path
import re

class Spectrum():
    def __init__(self):
        self.Charge=2
        self.Name=''
        self.Uni_Pet=1
        self.Samples=1
        self.Score=0
        self.Evalue=0
        self.Condition='Distinct'
        self.Spec='#,'
        self.Peptide=''
        self.SampleID='1.pFind'
        self.Mod_Sites=''
        self.Calc_M=0
        self.Delta_M=0
        self.ppm=0
        self.Proteins=''

def getEvalue(list):
    return float(list[14])
def buildProteinList(proteinID,totalProteinList,f):
    proList=[]
    for i in range(6):
        proList.append(f.readline())
    totalProteinList.append((proList[0].split('='))[1].strip('\n'))
    if 'SameSet=0' not in proList[4]:
        sameSet=[]
        sameSet=((proList[4].split('='))[1]).split(',')
        for j in range(int(sameSet[0])):
            totalProteinList.append(sameSet[j+1])
    if 'SubSet=0' not in proList[5]:
        SubSet=[]
        SubSet=((proList[5].split('='))[1]).split(',')
        for j in range(int(SubSet[0])):
            totalProteinList.append(SubSet[j+1])
def buildSpectrumList(spectrumNum,totalList,f,count):
    specList=[]
    for i in range(7):
        specList.append(f.readline())  
    ValidCandidate=int(filter(str.isdigit,specList[6]))
    if ValidCandidate != 0:
        spectrumData=Spectrum() 
        spectrumData.Name=(specList[0].split('='))[1].strip('\n')
        spectrumData.Charge=int(filter(str.isdigit,specList[1]))
        noScore=0
        for i in range(ValidCandidate):
            noList=[]
            for j in range(15):
                noList.append(f.readline())
            noScoreList=noList[0].split('=')
            scoreList2=[]
            if float(noScoreList[1])>=noScore:
                scoreList2.append(noScoreList[1])
                noScore=float(noScoreList[1])
                spectrumData.Evalue=format(float((noList[1].split('='))[1].strip('\n')),'.2e')
                spectrumData.Peptide=(noList[5].split('='))[1].strip('\n')
                spectrumData.Mod_Sites=(noList[12].split('='))[1].strip('\n')
                spectrumData.Score=noScore
                spectrumData.Calc_M=float((noList[2].split('='))[1].strip('\n'))
                spectrumData.Proteins=(noList[7].split('='))[1].strip('\n')
            #for k in range(len(scoreList2)):
            #    if noScoreList[1] == k: 
            #        count +=1
          #  elif float(noScoreList[1])==noScore and 'REVERSE_' in spectrumData.Proteins:
          #      scoreList2.append(noScoreList[1])
          #      noScore=float(noScoreList[1])
          #      spectrumData.Evalue=format(float((noList[1].split('='))[1].strip('\n')),'.2e')
          #      spectrumData.Peptide=(noList[5].split('='))[1].strip('\n')
          #      spectrumData.Mod_Sites=(noList[12].split('='))[1].strip('\n')
          #      spectrumData.Score=noScore
          #      spectrumData.Calc_M=float((noList[2].split('='))[1].strip('\n'))
          #      spectrumData.Proteins=(noList[7].split('='))[1].strip('\n')
        totalList.append((str(spectrumNum)+' '+spectrumData.Name+' '+str(spectrumData.Uni_Pet)+' '+str(spectrumData.Samples)+' '+str(spectrumData.Score)+' '+spectrumData.Condition+' '+spectrumData.Spec+' '+spectrumData.Peptide+' '+str(spectrumData.SampleID)+' '+spectrumData.Mod_Sites+' '+str(spectrumData.Calc_M)+' '+str(spectrumData.Delta_M)+' '+str(spectrumData.ppm)+' '+spectrumData.Proteins+' '+str(spectrumData.Evalue)+' '+str(spectrumData.Charge)).split())
        print '[Spectrum'+str(spectrumNum)+']\tE-value='+str(spectrumData.Evalue)
    else:
       print '[Spectrum'+str(spectrumNum)+']\tValidCandidate=0\t已经过滤'

writeFile=open('data\out.txt','w')
writeFile.write('#\tSpectrum\tUni_Pet\tSamples\tScore\tCondition\n*\t#,Spec\tPeptide\tSampleID\tMod_Sites\tScore\tCalc_M\tDelta_M\tppm\tProteins\n')
writeFile.close()

readFile=open('data\data1.txt','r')
print '...londing data'
dic=defaultdict(list)
totalList=[]
totalProteinList=[]
count =0
FDR=0.01
while True:
    line=readFile.readline()
    if not line:
        break
    if '[Spectrum' in line:
        spectrumNum=int(filter(str.isdigit,line))
        buildSpectrumList(spectrumNum,totalList,readFile,count)
    if '[Protein' in line:
        proteinID=int(filter(str.isdigit,line))
        buildProteinList(proteinID,totalProteinList,readFile)
for ti in range(len(totalList)):
    if int(totalList[ti][15]) == 2:
        dic['charge2'].append(totalList[ti])
    if int(totalList[ti][15]) == 3:
        dic['charge3'].append(totalList[ti])
    if int(totalList[ti][15]) == 4:
        dic['charge4'].append(totalList[ti])
    if int(totalList[ti][15]) == 5:
        dic['charge5'].append(totalList[ti])
    if int(totalList[ti][15]) == 6:
        dic['charge6'].append(totalList[ti])
    if int(totalList[ti][15]) == 7:
        dic['charge7'].append(totalList[ti])
print 'sort dict...'
for di in range(2,8):
    dic['charge'+str(di)].sort(key=getEvalue)
print 'complete'
print '\n\n\n\n'
print 'calculating FDR and writing file...'
writeCount=0
for m in range(2,8):
    decoyCount=0
    targetCount=0
    tempFDR=0.0
    for i in range(len(dic['charge'+str(m)])):
        flag=0
        #proFlag=0
        tempProteins=dic['charge'+str(m)][i][13]
        proteinList=tempProteins.split(',')
        for j in range(int(proteinList[0])):
            #proFlag=0
            protemp=proteinList[j+1]
            if 'REVERSE_' in protemp:
                decoyCount += 1
                flag=1
            else:
                targetCount += 1
                #for pro in totalProteinList:
                    #if protemp == pro:
                        #proFlag=1
                        #break
        #if proFlag ==1:
        tempFDR=round(float(decoyCount)/(targetCount),2)
        if tempFDR<=FDR:
            if flag == 0:
                writeCount +=1
                print 'charge'+str(m)+'  FDR='+str(tempFDR)+'\twrite '+str(writeCount)+' line in file\tspectrumID='+dic['charge'+str(m)][i][0]
                writeFile=open('data\out.txt','a')
                writeStr=str(writeCount)+'\t'+dic['charge'+str(m)][i][1]+'\t'+dic['charge'+str(m)][i][2]+'\t'+dic['charge'+str(m)][i][3]+'\t'+dic['charge'+str(m)][i][14]+'\t'+dic['charge'+str(m)][i][5]+'\n*\t'+dic['charge'+str(m)][i][6]+'\t'+dic['charge'+str(m)][i][7]+'\t'+dic['charge'+str(m)][i][8]+'\t'+dic['charge'+str(m)][i][9]+'\t'+dic['charge'+str(m)][i][14]+'\t'+dic['charge'+str(m)][i][10]+'\t'+dic['charge'+str(m)][i][11]+'\t'+dic['charge'+str(m)][i][12]+'\t'+dic['charge'+str(m)][i][13]+'\n'
                writeFile.write(writeStr)
                writeFile.close()
        else:
            continue
print count

