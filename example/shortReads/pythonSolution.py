import sys, operator

inFile = open(sys.argv[1],'r')

header = ''
data = ''

allData = []

i=1

for line in inFile:
    if i%4==1:
        if data != '':
            allData.append((header,data))
        header = line.strip()[1:]
        data = line
    else:
        data += line
    i+=1

allData.append((header,data))
allData.sort(key = operator.itemgetter(0))

for item in allData:
    print (item[1].strip())

