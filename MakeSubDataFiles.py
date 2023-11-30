
#
folder = '/sciclone/scr10/dchendrickson01/EFIT/Double4m/'
fileO = 'Anima.csv'
files = []
sets = []
for i in range(10):
    j = (i+1) * 1000
    sets.append(j)
    files.append('AnimaClean'+str(i).zfill(5)+'.csv')

FileO= open(folder+fileO, 'r')
Lines = FileO.readlines()
FileO.close()

filesClean=[]
for file in files:
    filesClean.append(open(folder+'frames/'+file,'w'))
count = 0
for line in Lines:
    for i in range(10):
        if count == 0:
            filesClean[i].writelines(line)
        elif line[:len(str(sets[i]))+1] == str(set[i])+',':
            filesClean[i].writelines(line)
    count +=1
print(count)
for file in filesClean:
    file.close()