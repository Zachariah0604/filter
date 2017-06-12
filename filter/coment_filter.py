import re
import gc
import sys   
from pfind_filter import *
from fdr import *
sys.setrecursionlimit(1000000)
writeFile=open('data/coment_out.txt','w')
writeFile.write('Spectrum\tPeptide\tMod_Sites\tProteins\tEvalue\tCharge\n')
writeFile.close()
def select_spectrum(peptide_list,max_score,total_list,spectrum_num,spectrum_name,charge):
    select_list=[]
    for i in range(len(peptide_list)):
        if float(max_score) == float(peptide_list[i][0]):
            select_list.append(peptide_list[i])
    if len(select_list) > 1:
        for i in range(len(select_list)):
            proteinList=select_list[i][4].split(';')
            t_cunt=0
            proten_count=len(proteinList)
            for j in range(proten_count):
                protemp=proteinList[j]
                if 'REVERSE_' in protemp:
                    t_cunt += 1
            if t_cunt == proten_count:
                continue
            else:
                total_list.append((str(spectrum_num)+'\t'+spectrum_name+'\t'+select_list[i][2]+'\t'+select_list[i][3]+'\t'+select_list[i][4]+' '+str(select_list[i][1])+' '+str(charge)).split())
                #print('[Spectrum'+str(spectrum_num)+']\tE-value='+str(select_list[i][1]))
                break
    else:
        total_list.append((str(spectrum_num)+'\t'+spectrum_name+'\t'+select_list[0][2]+'\t'+str(select_list[0][3])+'\t'+select_list[0][4]+' '+str(select_list[0][1])+' '+str(charge)).split())
        #print('[Spectrum'+str(spectrum_num)+']\tE-value='+str(select_list[0][1]))
    return total_list
def get_spectrum(rf,spectrum_num,spectrum_name,charge,total_list):
    peptide_list=[]
    max_score=0.0
    while True:
        index=rf.tell()
        line=rf.readline()
        if not line:
            return rf.tell(),total_list
        if '.dta' in line:

            if len(peptide_list)>0:
                total_list=select_spectrum(peptide_list,max_score,total_list,spectrum_num,spectrum_name,charge)
            return index,total_list
        else:
            temp_list=[]
            
            temp_list.append(line.split()[5])
            if float(temp_list[0])>float(max_score):
                max_score=float(temp_list[0])
            temp_list.append(float(line.split()[6]))
            temp_list.append(line.split()[9].split('.')[1])
            temp_list.append(line.split()[10])
            line=rf.readline()
            temp_list.append(line.split()[1])
        peptide_list.append(temp_list)

def main(spath):
    total_list=[]
    rf=open(spath,'r')
    while True:
        line=rf.readline()
        if not line:
            break
        if '.dta' in line:
            spectrum_num=int(re.sub("\D", "",line.split()[1]))
            
            spectrum_name=line.split()[10]
            charge=line.split()[3]
            index,total_list=get_spectrum(rf,spectrum_num,spectrum_name,charge,total_list)
            rf.seek(index,0)
    rf.close()
    dic=get_spectrum_dict(total_list)
    filter_with_fdr(0.01,dic,'data/coment_out/filter_result/filter_coment_result.txt',0)
if __name__=='__main__':
    writeFile=open('data/coment_out/filter_result/filter_coment_result.txt','w')
    writeFile.write('Spectrum\tPeptide\tMod_Sites\tProteins\tEvalue\tCharge\n')
    writeFile.close()
    paths=folder('data\coment_out\coment_search_out')
    for spath in paths:
        main(spath)
        print(spath+' Complete filtration')