import pandas as pd

Output = pd.read_csv('NaturalHistoryOfSjgr_DATA_LABELS_2020-07-28_1317.csv')
Output = Output[Output['USA code'].notna()]
Output = Output[Output['Sex'] == 'Female']
Output = Output[Output['SjD Diagnosis'] != 'Unknown']
Output.index = Output['USA code']
meta = pd.read_csv('2023-06-20_MITO_RNASEQ-METADATA_forThomas.csv',index_col=8)
meta['ANA'] = (meta['ANA'] > 0).astype('str').replace('False','Negative').replace('True','Positive')
meta['RF'] = (meta['RF'] > 0).astype('str').replace('False','Negative').replace('True','Positive')
Output = pd.concat([Output,meta])[['Age','SjD Diagnosis','Focus Score','SSA','SSB','ANA','RF']]
Output['SjD Diagnosis'] = Output['SjD Diagnosis'].replace(['SjD','SS1','SS2'],'Positive').replace(['nonSjD','HV','NSS'],'Negative')
Output['SSA'] = Output['SSA'].replace(1,'Positive').replace(0,'Negative')
Output['SSB'] = Output['SSB'].replace(1,'Positive').replace(0,'Negative')
Output.to_csv('meta.csv')