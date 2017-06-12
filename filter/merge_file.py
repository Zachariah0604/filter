def merge(merge_file,pfind_source_file,mascot_source_file,coment_source_file):
    writeFile=open(merge_file,'w')

    with open(pfind_source_file,'r') as pf:
        pf.readline()
        while True:
            line = pf.readline()
            if not line :
                break
            
            write_str=line.split()[0]+'\t'+line.split()[1]+'\t'+line.split()[4]+'\t'+line.split()[5]+'\n'
            writeFile.write(write_str)
        pf.close()
    print('pfind compelte')
    with open(mascot_source_file,'r') as mf:
        mf.readline()
        while True:
            line=mf.readline()
            if not line:
                break
            write_str=(line.split()[0]+'\t'+line.split()[1])+'\n'
            writeFile.write(write_str)
        mf.close()
    print('mascot compelte')
    with open(coment_source_file,'r') as cf:
        cf.readline()
        while True:
            line=cf.readline()
            if not line:
                break
            write_str=(line.split()[0]+'\t'+line.split()[1])+'\n'
            writeFile.write(write_str)
        cf.close()
    print('coment compelte')
    writeFile.close()