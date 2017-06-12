from exclude_pif import *
import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

dicM={'A':71.03711,'C':160.03065,'D':115.02694,'E':129.04259,'F':147.06841,\
      'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,\
      'M':131.04048,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,\
      'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06332}

def spectrums_from_mgf():
    files=folder('data/LTQ-mgf')
    mgf_spectrums=[]
    mgf_charges=[]
    mgf_pepmass=[]
    ions=[]
    for mgf_file in files:
        with open(mgf_file,'r') as rf:
            while True:
                line = rf.readline()
                if not line:
                    break
                if 'TITLE=' in line:
                    mgf_spectrums.append(line.split('=')[1].strip('\n'))
                if 'CHARGE=' in line:
                    mgf_charges.append(line.split('=')[1].strip('\n').strip('+'))
                if 'PEPMASS=' in line:
                    mgf_pepmass.append(float(line.split('=')[1].strip('\n')))
                    str_=''
                    while True:
                        line=rf.readline()
                        if 'END IONS' in line:
                            break
                        str_+=line.split()[0]+','+line.split()[1].strip('\n')+';'
                    ions.append(str_)
        print(mgf_file+'\t complete')
    mgf_pd=pd.DataFrame({"Charge":mgf_charges,"Pepmass":mgf_pepmass,"Ions":ions},index=mgf_spectrums)
    return mgf_pd
def spectrum_from_filter():
    data=pd.read_csv('data\pfind_out\\filter_result\\filter_pfind_result.txt',delimiter='\t',header=0)
    dic = defaultdict(list)
    spectrums=[]
    b_ions=[]
    y_ions=[]
    for row in data.itertuples():
        charge=row.Charge
        if int(charge) == 2:
            spectrums.append(row.Spectrum)
            peptide=row.Peptide
            len_peptide=len(peptide)
            b=[]
            y=[]
            mass_of_peptide=0
            for i in range(len_peptide):
                mass_of_peptide+=dicM[peptide[i]]
            for i in range(1,len_peptide):
                if i==1:
                    b.append(dicM[peptide[i-1]])
                    y.append(mass_of_peptide-b[-1])
                else:
                    b.append(dicM[peptide[i-1]]+b[i-2])
                    y.append(mass_of_peptide-b[-1])
            b_ions.append(np.array(b)+1)
            y_ions.append(np.array(y)+19)
    allen=len(spectrums)
    filter_spectrums=pd.DataFrame({"Spectrum":spectrums,"Bions":b_ions,"Yions":y_ions})
    return filter_spectrums,allen
def get_scatter_list(scatter_list,Ions,ion_mass,bory):
    for index in range(len(Ions)):
            min_dvalue=0.5
            flag=0
            for i_mass in ion_mass:
                devalue=float(Ions[index])-float(i_mass)
                if abs(devalue)<0.5 and devalue < min_dvalue:
                    flag=1
                    min_dvalue=devalue
                elif devalue < -0.5:
                    break
            if flag:
                scatter_list.append([bory*(index+1),min_dvalue])
    return scatter_list
def get_dvalue(mgf_pd,fileter_spectrums,allen):
    scatter_list=[]
    pecs=0.1
    index=0
    for row in fileter_spectrums.itertuples():
        index+=1
        str_ions=mgf_pd.loc[row.Spectrum,'Ions']
        ions=str_ions.split(';')
        ion_mass=[]
        for ion in ions:
            if len(ion)>0:
                ion_mass.append(ion.split(',')[0])
        
        scatter_list=get_scatter_list(scatter_list,row.Bions.tolist(),ion_mass,1)
        scatter_list=get_scatter_list(scatter_list,row.Yions.tolist(),ion_mass,-1)
       
        if int(allen*pecs) ==index:
            print('\t..'+str(round(pecs*100,2))+'% complete')
            pecs+=0.1
    return scatter_list
def plot_scatter(scatter_list):
    #f1 = plt.figure(1)
    #plt.subplot(211)
    x1=scatter_list[:,1]
    x2=scatter_list[:,0]
    plt.scatter(x1,x2,s=0.01)
    plt.show()

if __name__=='__main__':
    print('get spectrums from source file..')
    mgf_pd=spectrums_from_mgf()
    print('get sprctrums from filtered file..')
    filter_spectrums,allen=spectrum_from_filter()
    print('get all the coordinate list')
    scatter_list=get_dvalue(mgf_pd,filter_spectrums,allen)
    print('plot..')
    plot_scatter(np.array(scatter_list))
    i=0