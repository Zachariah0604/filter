import shutil
import os
def folder(dir):
    paths=[]
    for lists in os.listdir(dir):
        path =os.path.join(dir,lists)
        paths.append(path)
        if os.path.isdir(path):
            folder(path)
    return paths
if __name__ == '__main__':
    paths=folder('data/LTQ-mgf/')
    for path in paths:
        with open(path,'r') as rf:
            with open(path+'.new.mgf','w') as wf:
                for line in rf.readlines():
                    if 'PIF=' not in line:
                        wf.write(line)
        shutil.move(path+'.new.mgf', path)