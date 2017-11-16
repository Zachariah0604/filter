from main import *
class Spectrum():
    def __init__(self):
        self.Charge=2
        self.Spectrum_Name=''
        self.Evalue=0
        self.Peptide=''
        self.Mod_Sites=''
        self.Proteins=''
writeFile=open('data/simple_out.txt','w')
writeFile.write('Spectrum\tPeptide\tMod_Sites\tProteins\tEvalue\tCharge\n')
writeFile.close()
total_list=[]
with open('data//2017_05_19_16_52_12.spectra.txt','r') as rf:
    line=rf.readline()
    line=rf.readline()
    while True:
        line1=rf.readline()
        line2=rf.readline()
        if not line2:
            break
        spectrum_num=line1.split()[0]
        spectrum_data=Spectrum()
        spectrum_data.Spectrum_Name=line1.split()[1]
        spectrum_data.Evalue=line1.split()[4]
        spectrum_data.Peptide=line2.split()[2]
        spectrum_data.Mod_Sites=line2.split()[4]
        spectrum_data.Proteins=line2.split()[9]
        spectrum_data.Charge=spectrum_data.Spectrum_Name.split('.')[3]

        ##
        spec_str=str(spectrum_num)+'\t'+spectrum_data.Spectrum_Name+'\t'+spectrum_data.Peptide+'\t'+spectrum_data.Mod_Sites+'\t'+spectrum_data.Proteins+'\t'+str(spectrum_data.Evalue)+'\t'+str(spectrum_data.Charge)+'\n'
        total_list.append(spec_str.split())
    rf.close()
dic=get_spectrum_dict(total_list)
for m in range(2,8):
    for i in range(len(dic['charge'+str(m)])):
        writeFile=open('data\simple_out.txt','a')
        writeStr=dic['charge'+str(m)][i][1]+'\t'+dic['charge'+str(m)][i][2]+'\t'+dic['charge'+str(m)][i][3]+'\t'+dic['charge'+str(m)][i][4]+'\t'+dic['charge'+str(m)][i][5]+'\t'+dic['charge'+str(m)][i][6]+'\n'
        writeFile.write(writeStr)
        writeFile.close()

        
       