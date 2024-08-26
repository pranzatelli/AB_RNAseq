# coding=utf-8
from os import system

VARIABLES = ['Age','SjD Diagnosis','Focus Score','SSA','SSB','ANA','RF']
REGRESSIONS = [1,0,1,0,0,0,0]

with open('swarmfile.sh','w') as swarmfile:
	for e in range(2):
		for b in range(2):
			for t in range(2):
				for l in range(2):
					system('mkdir Summary/E'+str(e)+'B'+str(b)+'T'+str(t)+'L'+str(l))
					for i in range(len(VARIABLES)):
						label = VARIABLES[i].replace(' ','_').replace('/','+').replace('(','{').replace(')','}').replace("'",'?').replace('&','ö').replace('<','ø')
						swarmfile.write('source /data/ChioriniCompCor/tools/conda/etc/profile.d/conda.sh; conda activate python39; python swarMORFR.py '+label+' '+str(REGRESSIONS[i])+' '+str(e)+' '+str(b)+' '+str(t)+' '+str(l)+' > Summary/E'+str(e)+'B'+str(b)+'T'+str(t)+'L'+str(l)+'/'+label+'.txt\n')
