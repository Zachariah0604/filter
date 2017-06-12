import math
import re
from pfind_filter import *
from collections import defaultdict
import urllib.parse



def getEvalue(list):
    return float(list[5])
def get_peptide(str_):
    
    str_list=str_.split('=')
    str_peptide=str_list[1].split(';')
    str_proteins=str_peptide[1]
    peptide=[]
    peptide.append(str_list[0].split('_')[0])
    peptide.append(str_list[0].split('_')[1])
    peptide.append(str_peptide[0].split(',')[4])
    peptide.append(str_peptide[0].split(',')[6])
    peptide.append(str_proteins.strip('\n'))
    peptide.append(str_peptide[0].split(',')[7])

    return peptide
def get_merge_peptides(rf):
    temp_query_number = ''
    temp_list = []
    peptides = []
    line=rf.readline()
    while True:
        line=rf.readline()
        if '--' in line:
            break
        if [s for s in ['_terms','_subst','p1=-1'] if s in line]:
            continue
        else:
            peptide= get_peptide(line)
       
            query_num = peptide[0]
            if query_num != temp_query_number:
                if temp_query_number == '':
                    temp_list.append(peptide)
                else:
                    peptides.append(temp_list)
                    temp_list = []
                    temp_list.append(peptide)
                temp_query_number = query_num
            else:
                temp_list.append(peptide)
    peptides.append(temp_list)
    return peptides
def get_peptides(merge_petide,peptides):
    for i in range(len(merge_petide)):
        score_list=[]
        select_list=[]
        for j in range(len(merge_petide[i])):
            score_list.append(float(merge_petide[i][j][5]))
        max_socre=max(score_list)
        for j in range(len(merge_petide[i])):
            if max_socre == float(merge_petide[i][j][5]):
                select_list.append(merge_petide[i][j])
        if len(select_list) > 1:
            mmm=0
            for m in range(len(select_list)):
                proteinList=select_list[m][4].split(',')
                t_cunt=0
                proten_count=int(len(proteinList))
                for n in range(proten_count):
                    protemp=proteinList[n]
                    if 'REVERSE_' in protemp:
                        t_cunt += 1
                if t_cunt == proten_count:
                    continue
                else:
                    peptides.append((select_list[m][0]+'\n'+select_list[m][2]+'\t'+select_list[m][3]+'\t'+select_list[m][4]+'\t'+str(select_list[m][5])).split())
                    break
        else:
            peptides.append((select_list[0][0]+'\n'+select_list[0][2]+'\t'+select_list[0][3]+'\t'+select_list[0][4]+'\t'+str(select_list[0][5])).split())
    return peptides
def get_spectrum(str_,rf):
    spectrum=[]
    spectrum.append(re.sub("\D","",str_))
    list_=[]
    for i in range(3):
        line=rf.readline()
        list_.append(line)
    spectrum.append(urllib.parse.unquote(list_[0].split('=')[1].strip('\n')))
    spectrum.append(list_[2].split('=')[1].strip('\n'))
    return spectrum
def get_spectrum_dict(peptides,spectrums):
    dic=defaultdict(list)
    for ti in range(len(peptides)):
        spec_index=int(re.sub("\D","",peptides[ti][0]))-1
        charge=int(re.sub("\D","",spectrums[spec_index][2]))
        if charge== 2:
            dic['charge2'].append([peptides[ti][0],spectrums[spec_index][1],peptides[ti][1],peptides[ti][2],peptides[ti][3],float(peptides[ti][4]),charge])
        if charge == 3:
            dic['charge3'].append([peptides[ti][0],spectrums[spec_index][1],peptides[ti][1],peptides[ti][2],peptides[ti][3],float(peptides[ti][4]),charge])
    #print('sort dict...')
    for di in range(2,4):
        dic['charge'+str(di)].sort(key=getEvalue,reverse = True)
    #print('complete')
    return dic
def main(spath):
    writeFile=open('data/mascot_out.txt','w')
    writeFile.write('Spectrum\tPeptide\tMod_Sites\tProteins\tEvalue\tCharge\n')
    writeFile.close()  
    peptides=[]
    spectrums=[]
    with open(spath,'r') as rf:
        while True:
            line=rf.readline()
            if not line:
                break
            if 'Content-Type: application/x-Mascot; name="peptides"' in line:
                merge_peptides=get_merge_peptides(rf)
                peptides=get_peptides(merge_peptides,peptides)

            if 'Content-Type: application/x-Mascot; name="query' in line:
                rf.readline()
                spectrums.append(get_spectrum(line,rf))
    rf.close()
    dic=get_spectrum_dict(peptides,spectrums)
    filter_with_fdr(0.01,dic,'data/mscot_out/filter_result/filter_mascot_result.txt',0)

if __name__=='__main__':
    writeFile=open('data/mscot_out/filter_result/filter_mascot_result.txt','w')
    writeFile.write('Spectrum\tPeptide\tMod_Sites\tProteins\tEvalue\tCharge\n')
    writeFile.close()
    paths=folder('data\mscot_out\mscot_search_out')
    for spath in paths:
        main(spath)
        print(spath+' Complete filtration')