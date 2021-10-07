import glob
import pandas as pd
import sys


#Get file names


pharmFiles = glob.glob(f'{sys.argv[1]}*.gninatypes')
fragFiles = glob.glob(f'{sys.argv[2]}*.gninatypes')

print(len(fragFiles))

pharmCol = []
fragCol = []
Col = []
for i in range(len(fragFiles)):
    pharmName = f'{sys.argv[1]}_{i}.gninatypes'
    fragName = f'{sys.argv[2]}_{i}.gninatypes'
    pharmCol.append(pharmName)
    fragCol.append(fragName)
    Col.append(1)

print(len(fragCol))
#Make pandas dataframe
dataOut = pd.DataFrame({'col':Col, 'fragCol':fragCol, 'pharmCol':pharmCol})

dataOut.to_csv(f'{sys.argv[3]}.types', sep = ' ', header = False, index = False)
