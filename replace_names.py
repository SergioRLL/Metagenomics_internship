import sys
import re

dic = {}
with open(sys.argv[1]) as f:
    for line in f:
        a = re.match(r'(\S*)\s*(\S*)',line)
        dic[a.group(1)] = a.group(2)


#Select only the ERR accession id in the character string and create a file with ERZ id and ERR id (prepare a file first)
'''
for i in dic:
    co = dic[i].find('ERR')
    fco = co+3
    for x in range(fco, len(dic[i])):
        if dic[i][fco+1].isdigit():
            fco = fco+1
        else:
            break
    print(i, dic[i][co:fco+1])
'''

#Print contents of dic:
'''
for i in dic:
    print(i, dic[i])
'''

with open(sys.argv[2]) as f:
    for line in f:
        for ele in dic:
            if ele in line:
                line = re.sub(ele, dic[ele], line)
        print(line)
