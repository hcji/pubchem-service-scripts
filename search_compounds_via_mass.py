import json
import requests
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

def get_cid_from_mass(mass, tol=0.5):
    min_mass = mass - tol
    max_mass = mass + tol
    
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={}[MIMass]:{}[MIMass]&retmode=json&retmax=100000'.format(min_mass, max_mass)
    res = requests.get(url)
    string = res.content.decode(encoding='gbk')
    dicts = json.loads(string)
    cids = dicts['esearchresult']['idlist']
    
    idstring = ''
    smiles = []
    inchikey = []
    all_cids = []
    
    for i, cid in enumerate(cids):
        idstring += ',' + str(cid)
        if ((i%100==99) or (i==len(cids)-1)):
            url_i = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + idstring[1:(len(idstring))] + "/property/InChIKey,CanonicalSMILES/JSON"
            res_i = requests.get(url_i, timeout=999)
            soup_i = BeautifulSoup(res_i.content, "html.parser")
            str_i = str(soup_i)
            properties_i = json.loads(str_i)['PropertyTable']['Properties']
            idstring = ''
            for properties_ij in properties_i:
                smiles_ij = properties_ij['CanonicalSMILES']
                if smiles_ij not in smiles:
                    smiles.append(smiles_ij)
                    inchikey.append(properties_ij['InChIKey'])
                    all_cids.append(str(properties_ij['CID']))
                else:
                    wh = np.where(np.array(smiles)==smiles_ij)[0][0]
                    all_cids[wh] = all_cids[wh] + ', ' + str(properties_ij['CID'])
    
    result = pd.DataFrame({'InChIKey': inchikey, 'SMILES': smiles, 'PubChem': all_cids})
    return result