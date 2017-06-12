# -- coding: utf-8 -*-

from collections import defaultdict
import os
import os.path
import re
from fdr import *
from exclude_pif import *
class Spectrum():
    def __init__(self):
        self.Charge=2
        self.Spectrum_Name=''
        self.Evalue=0
        self.Peptide=''
        self.Mod_Sites=''
        self.Proteins=''

def getEvalue(list):
    return float(list[5])

def build_spectrum_list(spectrum_num,total_list,f):
    spec_base_list=[]
    for i in range(7):
        spec_base_list.append(f.readline())  
    ValidCandidate=int(re.sub("\D", "", spec_base_list[6]))
    if ValidCandidate != 0:
        spectrum_data=Spectrum() 
        spectrum_data.Spectrum_Name=(spec_base_list[0].split('='))[1].strip('\n')
        spectrum_data.Charge=int(re.sub("\D", "", spec_base_list[1]))
        number_score=0
        number_list=[]
        for i in range(ValidCandidate):
            number_temp_list=[]
            for j in range(15):
                number_temp_list.append(f.readline())
            num_score_list=number_temp_list[0].split('=')
            spectrum_data.Evalue=format(float((number_temp_list[1].split('='))[1].strip('\n')),'.2e')
            spectrum_data.Peptide=(number_temp_list[5].split('='))[1].strip('\n')
            spectrum_data.Mod_Sites=(number_temp_list[13].split('='))[1].strip('\n')
            spectrum_data.Proteins=(number_temp_list[7].split('='))[1].strip('\n')
            number_list.append((str(num_score_list[1])+'\t'+str(spectrum_data.Evalue)+'\t'+spectrum_data.Peptide+'\t'+str(spectrum_data.Mod_Sites)+'\t'+spectrum_data.Proteins).split())
            if float(num_score_list[1])>=number_score:
                number_score=float(num_score_list[1])
        select_list=[]
        for i in range(ValidCandidate):
            if number_score == float(number_list[i][0]):
                select_list.append(number_list[i])
        if len(select_list) > 1:
            for i in range(len(select_list)):
                proteinList=select_list[i][4].split(',')
                t_cunt=0
                proten_count=int(proteinList[0])
                for j in range(proten_count):
                    protemp=proteinList[j+1]
                    if 'REVERSE_' in protemp:
                        t_cunt += 1
                if t_cunt == proten_count:
                    continue
                else:
                    total_list.append((str(spectrum_num)+'\t'+spectrum_data.Spectrum_Name+'\t'+select_list[i][2]+'\t'+select_list[i][3]+'\t'+select_list[i][4]+' '+str(select_list[i][1])+' '+str(spectrum_data.Charge)).split())
                    #print('[Spectrum'+str(spectrum_num)+']\tE-value='+str(spectrum_data.Evalue))
                    break
        else:
            total_list.append((str(spectrum_num)+'\t'+spectrum_data.Spectrum_Name+'\t'+select_list[0][2]+'\t'+select_list[0][3]+'\t'+select_list[0][4]+' '+str(select_list[0][1])+' '+str(spectrum_data.Charge)).split())
          #  print('[Spectrum'+str(spectrum_num)+']\tE-value='+str(spectrum_data.Evalue))

    #else:
       #print('[Spectrum'+str(spectrum_num)+']\tValidCandidate=0\tfilted')
    return total_list



def get_spectrum_dict(total_list):
    dic=defaultdict(list)
    for ti in range(len(total_list)):
        if int(total_list[ti][6]) == 2:
            dic['charge2'].append(total_list[ti])
        if int(total_list[ti][6]) == 3:
            dic['charge3'].append(total_list[ti])
    #print('sort dict...')
    for di in range(2,4):
        dic['charge'+str(di)].sort(key=getEvalue)
    #print('complete')
    return dic



def main(spath):
    
    total_list=[]
    rf=open(spath,'r')
    while True:
        line=rf.readline()
        if not line:
            break
        if '[Spectrum' in line:
            spectrum_num=int(re.sub("\D", "", line))
            total_list=build_spectrum_list(spectrum_num,total_list,rf)
    dic=get_spectrum_dict(total_list)
    filter_with_fdr(0.01,dic,'data/pfind_out/filter_result/filter_pfind_result.txt',1)

if __name__=='__main__':
    writeFile=open('data/pfind_out/filter_result/filter_pfind_result.txt','w')
    writeFile.write('Spectrum\tPeptide\tMod_Sites\tProteins\tEvalue\tCharge\n')
    writeFile.close()
    paths=folder('data\pfind_out\pfind_search_out\search_result')
    for spath in paths:
        main(spath)
        print(spath+' Complete filtration')