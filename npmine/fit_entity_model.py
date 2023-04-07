from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from collections import Counter
import os
import json

import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn import preprocessing
from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef
from sklearn import datasets, metrics, model_selection
import joblib

import matplotlib.pyplot as plt  # doctest: +SKIP
from matplotlib.backends.backend_pdf import PdfPages

def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)

def get_features(inchi, tp='descriptors'):
    mols = []
    for x in inchi:
        mols.append(Chem.MolFromInchi(x))

    if tp=='descriptors':
        nms=[x[0] for x in Descriptors._descList]

        calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
        descrs = [calc.CalcDescriptors(x) for x in mols]

        descrs = pd.DataFrame(descrs, columns=nms)
        return descrs
    elif tp=='fingerprint':
        # generate binary Morgan fingerprint with radius 2
        fp = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]
        descrs = rdkit_numpy_convert(fp)
        return descrs
    else:
        raise ValueError('Undefined type.')

def plot_roc(x_ts, y_ts, model, out):
    probs = model.predict_proba(x_ts)
    preds = probs[:,1]
    fpr, tpr, threshold = metrics.roc_curve(y_ts, preds)
    roc_auc = metrics.auc(fpr, tpr)

    # method I: plt
    #https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    with PdfPages(os.path.join(out, 'rf_roc.pdf')) as pdf:
        plt.title('Receiver Operating Characteristic')
        plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
        plt.legend(loc = 'lower right')
        plt.plot([0, 1], [0, 1],'r--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()


def fit_model(inchi, groups, out, tp='descriptors'):
    if not os.path.exists(out):
        os.mkdir(out)
    else:
        raise ValueError('Folder already exists.')

    descrs = get_features(inchi, tp='descriptors')
    if np.any(np.isnan(descrs)):
        # better way to replace?
        descrs[np.isnan(descrs)] = 0

    lb = preprocessing.LabelBinarizer()
    gr = np.unique(groups)
    lb.fit(gr)
    y = lb.transform(groups)
    keys = lb.transform(np.unique(groups))[:,0]
    labels = {}
    for k,v in zip(keys, gr):
        labels[int(k)] = v

    with open(os.path.join(out, "labels.json"), 'w+') as f:
        json.dump(labels, f)

    print(Counter(groups))

    seed = 42
    x_tr, x_ts, y_tr, y_ts = train_test_split(descrs, y, random_state=seed)

    cv = StratifiedKFold(n_splits=5, random_state=seed)

    # obtain scale object which can be further applied to scale any data to fit the training set
    # for binary?
    if tp=='fingerprint':
        scaler = preprocessing.StandardScaler()
    elif tp=='descriptors':
    # it is a good idea to save it for future use
        scaler = preprocessing.MinMaxScaler()

    scaler.fit(x_tr)
    x_tr = scaler.transform(x_tr)
    x_ts = scaler.transform(x_ts)
    joblib.dump(scaler, os.path.join(out, "scaler.pkl"),
                compress=3)

    # create grid search dictionary
    param_grid = {"max_features": [x_tr.shape[1] // 10, x_tr.shape[1] // 7, x_tr.shape[1] // 5, x_tr.shape[1] // 3],
                  "n_estimators": [100, 250, 500]}

    # setup model building
    m = GridSearchCV(RandomForestClassifier(), param_grid,
                     n_jobs=2, cv=cv, verbose=1)

    # run model building
    m.fit(x_tr, y_tr)

    print(m.best_params_)
    print(m.best_score_)

    # save model
    joblib.dump(m, os.path.join(out, "rf_model.pkl"), compress=3)

    pred_rf = m.predict(x_ts)

    print('Accuracy:', accuracy_score(y_ts, pred_rf))
    print('Matthews_Corr:', matthews_corrcoef(y_ts, pred_rf))
    print('Cohen_Kappa:', cohen_kappa_score(y_ts, pred_rf))

    pd.DataFrame(m.best_estimator_.feature_importances_,
                 index=descrs.columns).to_csv(os.path.join(out,
                 "feat_imp.tsv"), sep='\t', header=None)

    plot_roc(x_ts, y_ts, m, out)

def model_predict(input_model, inchi, tp='descriptors'):
    descrs = get_features(inchi, tp='descriptors')
    if np.any(np.isnan(descrs)):
        # better way to replace?
        descrs[np.isnan(descrs)] = 0
    scaler = joblib.load(os.path.join(input_model,"scaler.pkl"))
    clf = joblib.load(os.path.join(input_model,"rf_model.pkl"))
    x = scaler.transform(descrs)
    pred = clf.predict(x)
    with open(os.path.join(input_model, 'labels.json')) as f:
        labels = json.load(f)
    return [labels[str(x)] for x in pred]

def apply_model(acty, comp, aid):
    sel_acty = acty[acty['AID']==aid]
    print('Unique CID:', len(sel_acty['CID'].unique()))

    sel_pub = sel_acty[['CID', 'Bioactivity Outcome']]
    sel_pub = sel_pub[sel_pub['Bioactivity Outcome'].isin(['Inactive', 'Active'])]
    sel_pub = sel_pub[~sel_pub['CID'].duplicated()]


    sel_pub = pd.merge(sel_pub,
                       comp[~comp['pubchem'].duplicated()],
                       left_on='CID', right_on='pubchem', how='left')

    sel_pub.fillna('', inplace=True)
    print('Missing InChI before processing',
          sum(sel_pub['standardInChI']==''))

    for i in sel_pub.index:
        if sel_pub.loc[i, 'standardInChI']=='':
            sel_pub.loc[i, 'standardInChI'] = Chem.MolToInchi(Chem.MolFromSmiles(sel_pub.loc[i, 'smiles']))

    print('Missing InChI before processing',
          sum(sel_pub['standardInChI']==''))
    print('Fitting the model')
    fit_model(sel_pub['standardInChI'],
              sel_pub['Bioactivity Outcome'],
              out=aid)


