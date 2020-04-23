# Combine the two files to form one clinical data file

import pandas as pd
seperator="\t"
pathToCancer="C:/Users/tim.kuijpers/Desktop/EGM/MetaData_Cancer_matched_all.txt"
pathToSex="C:/Users/tim.kuijpers/Desktop/EGM/SubjectDemographicData.txt"

cancerData = pd.read_csv(pathToCancer,sep=seperator)
sexData = pd.read_csv(pathToSex,sep=seperator)

dataFrameMerged = pd.merge(cancerData,sexData,on="egm_id")