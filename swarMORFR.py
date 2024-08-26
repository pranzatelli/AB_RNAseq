# coding=utf-8
from os import system,devnull
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeRegressor,DecisionTreeClassifier
from sklearn.model_selection import ShuffleSplit
import sys
import json
import pickle
from warnings import simplefilter

def read_salmon(sample):
	with open('/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/'+sample+'/quant.sf') as openfile:
		lines = openfile.read().split('\n')[1:-1]
	lines = [line.split('\t') for line in lines]
	lines = {line[0].split('.')[0]:float(line[3]) for line in lines}
	return lines

def get_enst(sample_list):
	frameable = {}
	for sample in sample_list:
		try:
			frameable[sample] = read_salmon(sample)
		except:
			pass
	return pd.DataFrame(frameable)

def make_data(column,enst=False,lim=False):
	Output = pd.read_csv('meta.csv',index_col=0)
	Diagnosis = Output['SjD Diagnosis']
	Output = Output[column][Output[column].notna()]
	if enst:
		Input = get_enst(Output.index)
	else:
		Input = pd.read_csv('all_gene.csv',index_col=0)
		Input = Input.divide(Input.sum()) * 1000000
	Output = Output.loc[list(set(Input.columns).intersection(set(Output.index)))]
	Input = Input[Output.index].T
	Diagnosis = Diagnosis[Output.index]
	if lim:
		Input = Input.T.loc[((Input == 0).sum() == 0).values].T
	return Input,Output,Diagnosis

def tree_importances(model):
	node_depth = np.zeros(shape=model.tree_.node_count, dtype=np.int64)
	is_leaves = np.zeros(shape=model.tree_.node_count, dtype=bool)
	stack = [(0, -1)]
	while len(stack) > 0:
		node_id, parent_depth = stack.pop()
		node_depth[node_id] = parent_depth + 1
		if (model.tree_.children_left[node_id] != model.tree_.children_right[node_id]):
			stack.append((model.tree_.children_left[node_id], parent_depth + 1))
			stack.append((model.tree_.children_right[node_id], parent_depth + 1))
		else:
			is_leaves[node_id] = True
	impD = {}; occD = {}
	for i in range(model.tree_.node_count):
		if is_leaves[i]:
			pass
		else:
			if model.tree_.feature[i] not in impD:
				impD[model.tree_.feature[i]] = 0
				occD[model.tree_.feature[i]] = 0
			impD[model.tree_.feature[i]] += 0.5 ** node_depth[i]
			occD[model.tree_.feature[i]] += 1
	return impD,occD

def split(DF,Y,diag,frac=False,n=False,balance=False):
	diag = diag.loc[DF.index]
	ix1 = diag[diag == 'Positive']
	ix2 = diag[diag == 'Negative']
	if frac:
		if balance:
			ix1 = ix1.sample(frac=frac).index.tolist()
			ix2 = ix2.sample(frac=frac).index.tolist()
			ix = ix1+ix2
		else:
			ix = diag.sample(frac=frac).index.tolist()
	elif n:
		if balance:
			ix1 = ix1.sample(n=int(n*len(ix1)/len(diag)+0.5)).index.tolist()
			ix2 = ix2.sample(n=int(n*len(ix2)/len(diag)+0.5)).index.tolist()
			ix = ix1+ix2
		else:
			ix = diag.sample(n=n).index.tolist()
	nx = [x for x in DF.index if x not in ix]
	return DF.loc[ix],DF.loc[nx],Y.loc[ix],Y.loc[nx]

def build_model(DF,Y,diag,Correlations,reg=True,ntrees=10000,enriched=True,balance=True):
	print('Enriched: '+str(enriched)+' Balanced: '+str(balance))
	sys.stdout = open(devnull, "w"); simplefilter('ignore')
	R2 = []; Importances = np.zeros(len(Correlations)); Occurrences = np.zeros(len(Correlations))
	R = []; Rsum = 0
	testX,holdx,testY,holdy = split(DF,Y,diag,frac=0.5,balance=balance)
	for tree in range(ntrees):
		if enriched:
			Genes = np.random.choice(np.arange(len(Correlations)),1000,replace=False,p=Correlations/np.sum(Correlations))
		else:
			Genes = np.random.choice(np.arange(len(Correlations)),1000,replace=False)
		valx,trainx,valy,trainy = split(holdx[Correlations[Genes].index],holdy,diag,n=5,balance=balance)
		testx = testX[Correlations[Genes].index]
		if reg:
			model = DecisionTreeRegressor()
			model.fit(trainx,trainy)
			predy = model.predict(valx).transpose()
			score = np.nan_to_num(np.corrcoef(valy,predy))[1][0]
		else:
			model = DecisionTreeClassifier()
			model.fit(trainx,trainy)
			score = model.score(valx,valy)
		score = score * score
		impD,occD = tree_importances(model)
		for gene in impD.keys():
			Importances[Genes[gene]] += score * impD[gene]
			Occurrences[Genes[gene]] += occD[gene]
		R.append(model.predict(testx) * score); Rsum += score	
	R = np.array(R).sum(axis=0) / Rsum
	Score = np.corrcoef(R,testY)[1][0]
	R2.append(str(Score * Score))
	sys.stdout = sys.__stdout__
	print(' '.join(R2))
	return Importances, Occurrences


def __main__():
	label = sys.argv[1].replace('_',' ').replace('+','/').replace('{','(').replace('}',')').replace('?',"'").replace('ö','&').replace('ø','<')
	reg = bool(int(sys.argv[2])); enr = bool(int(sys.argv[3])); bal = bool(int(sys.argv[4])); enst = bool(int(sys.argv[5])); lim = bool(int(sys.argv[6]))
	I,O,D = make_data(label,enst,lim)
	if not reg:
		ix = O.index
		O = LabelEncoder().fit_transform(O)
		O = pd.Series(O)
		O.index = ix
	Correlations = I.corrwith(O,method='kendall').dropna()
	Correlations = Correlations * Correlations
	for i in range(100):
		Importances,Occurrences = build_model(I,O,D,Correlations,reg=reg,enriched=enr,balance=bal)
		if enr and bal and not lim:
			C = pd.DataFrame(Correlations,columns=['Correlations'])
			C['Importances'] = Importances
			C['Occurrences'] = Occurrences
			C['Importances'] = C['Importances'] / C['Occurrences']
			C.to_csv('Summary/'+sys.argv[1]+'_'+str(enst)+'_'+str(i)+'.csv')


if __name__ == '__main__':
	__main__()

