from tri_result import *
if __name__ == '__main__':
    total_list=[]
    with open('data/pfind_out/filter_result/filter_pfind_result.txt','r') as rf:
        rf.readline()
        while True:
            line =rf.readline()
            if not line:
                break
            list_=line.split()
            temp=[]
            temp.append(list_[0])
            temp.append(list_[1])
            temp.append(list_[4])
            temp.append(list_[5].strip('\n'))
            total_list.append(temp)
    dic=get_spectrum_dict(total_list)
    get_redundancy(dic,'data/pfind_redundancy.txt')
    print('complete')