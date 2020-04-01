from rdkit import DataStructs
from rdkit import Chem
import pandas as pd

hits_df = pd.read_csv('/home/alejandro/Escritorio/covid19/hits.csv')
all_act_smiles = list(hits_df.loc[:,'Compound SMILES'].values)

hits = {}
for act in all_act_smiles:
    rdkitmol = Chem.MolFromSmiles(act)
    fp = Chem.RDKFingerprint(rdkitmol)
    hits[act] = (rdkitmol, fp)

ok_drugs = {}
fda_drugs = Chem.SDMolSupplier('FDA_approved_Drugs_Collection_682cmpds.sdf')
for drug in fda_drugs:
   if drug is not None:
       fdaid = drug.GetProp('BRC_ID')
       fp = Chem.RDKFingerprint(drug)
       ok_drugs[fdaid] = (drug, fp) 


results = {}
for hit in hits:
    for drug in ok_drugs:
        sim = DataStructs.FingerprintSimilarity(hits[hit][1],ok_drugs[drug][1])
        try:
            results[hit].append((drug,sim))
        except KeyError:
            results[hit] = [(drug,sim)]


#show best FDAdrug for fragment hit:
for hit in hits:
    print('Fragment hit:',hit)
    ranking = sorted(results[hit],key=lambda x:x[1],reverse=True)
    for rank in ranking[:5]:
        print(rank)

