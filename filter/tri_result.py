import pandas as pd
from collections import defaultdict
from merge_file import *
def getEvalue(list):
    return float(list[2])
def get_feature_list(file):
    spectrums=[]
    peptides=[]
    evalues=[]
    charges=[]
    with open(file,'r') as rf:
        rf.readline()
        while True:
            line = rf.readline()
            if not line:
                break
            
            spectrums.append(line.split()[0])
            peptides.append(line.split()[1])
            if len(line.split())>2:
                evalues.append(line.split()[2])
                charges.append(line.split()[3])
            else:
                evalues.append('-1')
                charges.append('-1')
        rf.close()
    return spectrums,peptides,evalues,charges
def write_result(file,spectrums,peptides,evalues,charges):
    wf=open(file,'w')
    result_list=[]
    spec_len=len(spectrums)
    print('number of spectrums:' +str(spec_len))
    cunt =0
    pecs =0.01
    for spectrum in spectrums:
        cunt+=1
        if spectrums.count(spectrum)==3 and spectrum != '-1':
            index=[u for u,v in enumerate(spectrums) if v==spectrum]
            if peptides[index[0]] ==peptides[index[1]] and peptides[index[1]]==peptides[index[2]]:
                wf.write(spectrum+'\t'+peptides[index[0]]+'\t'+evalues[index[0]]+'\t'+charges[index[0]]+'\n')
                result_list.append((spectrum+'\t'+peptides[index[0]]+'\t'+evalues[index[0]]+'\t'+charges[index[0]]).split())
                #print(spectrum)
                for m in index:
                    spectrums[m] = '-1'
                    peptides[m]= '-1'
        if cunt == int(spec_len*pecs):
            print('\t..'+str(round(pecs*100,2))+'% complete')
            pecs += 0.01
    wf.close()
    return result_list 

def get_spectrum_dict(total_list):
    dic=defaultdict(list)
    for ti in range(len(total_list)):
        if int(total_list[ti][3]) == 2:
            dic['charge2'].append(total_list[ti])
        if int(total_list[ti][3]) == 3:
            dic['charge3'].append(total_list[ti])
    for di in range(2,4):
        dic['charge'+str(di)].sort(key=getEvalue)
   
    return dic

def get_redundancy(dic,file): 
    wf=open(file,'w')
    redundancy=[]
    for m in range(2,4):
        peptides_le=[]
        for i in range(len(dic['charge'+str(m)])):
            peptides_le.append(dic['charge'+str(m)][i][1])
        for i in range(len(dic['charge'+str(m)])):
            if dic['charge'+str(m)][i][1]!='-1':
                if peptides_le.count(dic['charge'+str(m)][i][1]) == 1:
                    redundancy.append(dic['charge'+str(m)][i])
                    wf.write(dic['charge'+str(m)][i][0]+'\t'+dic['charge'+str(m)][i][1]+'\t'+dic['charge'+str(m)][i][2]+'\t'+dic['charge'+str(m)][i][3]+'\n')
                else:
                    index=[u for u,v in enumerate(peptides_le) if v==dic['charge'+str(m)][i][1]]
                    redundancy.append(dic['charge'+str(m)][index[0]])
                    wf.write(dic['charge'+str(m)][index[0]][0]+'\t'+dic['charge'+str(m)][index[0]][1]+'\t'+dic['charge'+str(m)][index[0]][2]+'\t'+dic['charge'+str(m)][index[0]][3]+'\n')
                    for n in index[1:]:
                        dic['charge'+str(m)][n][1]='-1' 
                        peptides_le[n] ='-1'
    wf.close()
def main():
    print('..merge file')
    merge('data/merge_file.txt','data/pfind_out/filter_result/filter_pfind_result.txt','data/mscot_out/filter_result/filter_mascot_result.txt','data/coment_out/filter_result/filter_coment_result.txt')
    print('..find the common spectrum with same peptide of the three files')
    spectrums,peptides,evalues,charges=get_feature_list('data/merge_file.txt')
    result_list=write_result('data/result.txt',spectrums,peptides,evalues,charges)
    print('..redundancy')
    dic=get_spectrum_dict(result_list)
    get_redundancy(dic,'data//redundancy_result.txt')
    print('complete')
if __name__ == '__main__':
    main()