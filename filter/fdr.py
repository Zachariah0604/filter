import os
def write_file(write_count,dic,m,i,temp_fdr,file_name):
    #print('charge'+str(m)+'  FDR='+str(temp_fdr)+'\twrite '+str(write_count)+' line in file\tspectrumID='+dic['charge'+str(m)][i][0])
    writeFile=open(file_name,'a')
    writeStr=dic['charge'+str(m)][i][1]+'\t'+dic['charge'+str(m)][i][2]+'\t'+dic['charge'+str(m)][i][3]+'\t'+dic['charge'+str(m)][i][4]+'\t'+str(dic['charge'+str(m)][i][5])+'\t'+str(dic['charge'+str(m)][i][6])+'\n'
    writeFile.write(writeStr)
    writeFile.close()
def filter_with_fdr(fdr,dic,file_name,role):
    write_count=0
    for m in range(2,3):
        decoyCount=0
        targetCount=0
        temp_fdr=0.0
        for i in range(len(dic['charge'+str(m)])):
            try:
                flag=0
                tempProteins=dic['charge'+str(m)][i][4]
                proteinList=tempProteins.split(',')
                t_cunt=0
                if role==1:
                    proten_count=int(proteinList[0])
                elif role == 0:
                    proten_count=int(len(proteinList))
                if proten_count > 1:
                    for j in range(proten_count):
                        protemp=proteinList[j+role]
                        if 'REVERSE_' in protemp:
                            t_cunt += 1
                    if t_cunt == proten_count:
                        decoyCount+=1
                        flag=1
                    else:
                        targetCount+=1
                else:
                    protemp=proteinList[role]
                    if 'REVERSE_' in protemp:
                        decoyCount+=1
                        flag=1
                    else:
                        targetCount+=1
                #print(targetCount)
                #print(decoyCount)
                temp_fdr=float(decoyCount)/(targetCount)
                #print(temp_fdr)
               # print(temp_fdr)
                #os.system('pause')
            except:
                print(file_name+'\t'+str(i))
            if temp_fdr<=float(fdr):
                if flag == 0:
                    write_count +=1
                    write_file(write_count,dic,m,i,temp_fdr,file_name)
            else:
                continue
