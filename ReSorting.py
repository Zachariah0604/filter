def getValue(list):
    return float(list[2])

f=open("data/out.txt",'r')
list=[]
while True:
    line=f.readline()
    if not line:
        break
    if '20091211_' in line:
        tmpList=line.split('\t')
        list.append((tmpList[0]+' '+tmpList[1]+' '+tmpList[4]).split())
#        i=1
list.sort(key=getValue)

writeFile=open("data/out2.txt",'w')
for i in list: 
    writeStr=''
    for j in i:
        writeStr+=j+'\t'
    writeStr+='\n'
    writeFile.write(writeStr)
writeFile.close()
i=1
