import pandas as pd
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

Output = pd.read_csv('NaturalHistoryOfSjgr_DATA_LABELS_2020-07-28_1317.csv')
Output = Output[Output['USA code'].notna()]
Output = Output[Output['Sex'] == 'Female']
Output = Output[Output['SjD Diagnosis'] != 'Unknown']
Output.index = Output['USA code']
meta = pd.read_csv('2023-06-20_MITO_RNASEQ-METADATA_forThomas.csv',index_col=8)
meta = meta[meta['Sex'] == 'FEMALE']
meta['ANA'] = (meta['ANA'] > 0).astype('str').replace('False','Negative').replace('True','Positive')
meta['RF'] = (meta['RF'] > 0).astype('str').replace('False','Negative').replace('True','Positive')
Output = pd.concat([Output,meta])
Output['SjD Diagnosis'] = Output['SjD Diagnosis'].replace(['SjD','SS1','SS2'],'Positive').replace(['nonSjD','HV','NSS'],'Negative')
Output['SSA'] = Output['SSA'].replace(1,'Positive').replace(0,'Negative')
Output['SSB'] = Output['SSB'].replace(1,'Positive').replace(0,'Negative')
Output['Cohort'] = [x[:3] for x in Output.index]
samples = glob('/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/BRA*')+glob('/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/MSG_*')
samples = [s.split('/')[-1] for s in samples]
inter = list(set(samples).intersection(set(Output.index)))
Output = Output.loc[inter]
WUS = pd.read_csv('RNAseq_WUS-WSSF-EYE.csv',index_col=0)
WUS.index = ['MSG_' + str(i) for i in WUS.index]
WUS = WUS.loc[list(set(WUS.index).intersection(set(Output.index)))]

def min_or_avg(label1,label2):
	brz = Output[[label1+" (OS)",label1+" (OD)"]]
	amr = WUS[[label2+'_R',label2+'_L']]
	plt.clf()
	brz_min = brz.min(axis=1); amr_min = amr.min(axis=1)
	brz_avg = brz.mean(axis=1); amr_avg = amr.mean(axis=1)
	brz_min = brz_min / brz_min.max(); amr_min = amr_min / amr_min.max()
	brz_avg = brz_avg / brz_avg.max() ; amr_avg = amr_avg / amr_avg.max()
	DF = pd.DataFrame([amr_min,amr_avg,brz_min,brz_avg],index=['American Minimum','American Average','Brazilian Minimum','Brazilian Average']).T
	DF['Sample'] = DF.index
	DF = DF.melt(id_vars='Sample')
	sns.kdeplot(DF,hue='variable',x='value',common_norm=False)
	plt.savefig('kde_plots/'+label2+'.pdf',format='pdf',bbox_inches='tight')
	return pd.concat([brz_avg,amr_avg]).dropna()

Output['Schirmer'] = min_or_avg("Schirmer's Test",'SCH')
min_or_avg("Schirmer's Test",'SCHWA')
Output['Tear Film Breakup Time'] = min_or_avg("Tear Film Breakup Time",'TBUT')
Output['Rose Bengal/van Bijsterveld'] = min_or_avg("Rose bengal score",'VANB')

def compare_data(label,inp1,inp2):
	Output[label] = pd.concat([inp1,inp2]).dropna()
	plt.clf()
	sns.kdeplot(Output,hue='Cohort',x=label,common_norm=False)
	plt.savefig('kde_plots/'+label+'.pdf',format='pdf',bbox_inches='tight')

compare_data('Salivary Flow',Output['Salivary Flow (ml/min)'],WUS['WUS_int']/15)

Labs = pd.read_csv('Labs_For_Thomas.csv',encoding_errors='replace')
Labs = Labs[pd.to_numeric(Labs['Observation Value'],errors='coerce').notnull()]
Labs['Observation Value'] = Labs['Observation Value'].astype('float')
Labs = Labs.pivot_table(columns='Observation Name',values='Observation Value',index='SPIT',aggfunc='mean')
Labs.index = ['MSG_'+str(int(i)) for i in Labs.index]

compare_data('C3',Output['C3 G/L']*100,Labs['C3 Complement'])
compare_data('C4',Output['C4 G/L']*100,Labs['C4 Complement'])
compare_data('IgG',Output['IgG  (700-1600 mg/dl)'],Labs['IgG Total'])
compare_data('Lactate Dehydrogenase',Output['Lactato Desidrogenase - LDH']/2,Labs['Lactate Dehydrogenase'])
compare_data('Gamma Fraction',Output['Gamma Fraction EFP'],Labs['Gamma'])
compare_data('B2M Fraction',Output['BETA 2 MICROGLOBULINAS']/10000,Labs['Beta 2'])

Output = Output[['Age','SjD Diagnosis','Focus Score','SSA','SSB','ANA','RF','Schirmer','Tear Film Breakup Time','Rose Bengal/van Bijsterveld','Salivary Flow','C3','C4','IgG','Lactate Dehydrogenase','Gamma Fraction','B2M Fraction']]
Output['Cohort'] = [i[:3] for i in Output.index]
Output.to_csv('meta.csv')


