import re
import time
import requests
from rdkit import Chem

"""
PubChem Smiles Standardization:
author: Jeff van Santen
date: Nov. 1, 2018
Usage:
To use this tool from a script, import it as follows
`from pubchem_smiles_standardizer import get_standardized_smiles`
`get_standarided_smiles` is the main function to use from this script.
INPUT:
in_smiles - a SMILES string representing a compound to standardize the structure of
OPTIONAL INPUT:
max_retry - DEFAULT = 3 - Can increase the number of times this function tries to get
                          results from PubChem service
OUTPUT:
out_smiles - A SMILES string representing a standardized compound
EXCEPTIONS:
TypeError       - Uses rdkit check if SMILES string is a compound raise TypeError if not
                  rdkit may also raise a C++ exception which can be also be caught with TypeError
ValueError      - If unable to to reach PubChem for some unknown reason
ConnectionError - If requests library is unable to reach PubChem
(Note the last two can be considered as redundant)
"""

def get_structural_similar_smiles(in_smiles, max_retry=3):
    print("Input SMILES:\t%s" % in_smiles)
    # Check that smiles is smile string like
    if not is_smiles(in_smiles):
        raise TypeError

    request1 = requests.post(url=urlSend, data=pubchem_structural_similar_string.format(in_smiles))
    if request1.status_code == requests.codes.ok:
        reqid = get_PCT_reqid(request1.text)
        smiles = None
        counter = 0
        while not smiles and counter < max_retry:
            request2 = poll_PCT(reqid)
            # Check connection was made
            if request2.status_code != requests.codes.ok:
                print("There was a failure contacting PUG gateway.")
                break

            # Check status code
            # This will be "success" if queued or done otherwise break
            if not check_PCT_status(request2.text):
                print("PUG was unable to search the given structure")
                break

            # Try and get smiles string from request
            # smiles will either be retrieved or set to None, continuing loop
            smiles = get_PCT_smiles(request2.text)

            # Sleep for a second if didn't get smiles
            time.sleep(1)
            counter += 1
        if not smiles:
          smiles = in_smiles
        print("Output SMILES:\t%s" % smiles)
        return smiles
    else:
        print("There was a failure contacting PUG gateway.")
        raise ValueError

 # DO NOT TOUCH BELOW THIIS
pubchem_structural_similar_string = \
"""<?xml version="1.0"?>
<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "https://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
<PCT-Data>
  <PCT-Data_input>
    <PCT-InputData>
      <PCT-InputData_query>
        <PCT-Query>
          <PCT-Query_type>
            <PCT-QueryType>
              <PCT-QueryType_css>
                <PCT-QueryCompoundCS>
                  <PCT-QueryCompoundCS_query>
                    <PCT-QueryCompoundCS_query_data>{0}</PCT-QueryCompoundCS_query_data>
                  </PCT-QueryCompoundCS_query>
                  <PCT-QueryCompoundCS_type>
                    <PCT-QueryCompoundCS_type_similar>
                      <PCT-CSSimilarity>
                        <PCT-CSSimilarity_threshold>80</PCT-CSSimilarity_threshold>
                      </PCT-CSSimilarity>
                    </PCT-QueryCompoundCS_type_similar>
                  </PCT-QueryCompoundCS_type>
                  <PCT-QueryCompoundCS_results>2000000</PCT-QueryCompoundCS_results>
                </PCT-QueryCompoundCS>
              </PCT-QueryType_css>
            </PCT-QueryType>
          </PCT-Query_type>
        </PCT-Query>
      </PCT-InputData_query>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data>
"""

pubchem_poll_string = \
"""<?xml version="1.0"?>
<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
<PCT-Data>
  <PCT-Data_input>
    <PCT-InputData>
      <PCT-InputData_request>
        <PCT-Request>
          <PCT-Request_reqid>{0}</PCT-Request_reqid>
          <PCT-Request_type value="status"/>
        </PCT-Request>
      </PCT-InputData_request>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data>"""

urlSend = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"

def get_PCT_reqid(request_text):
    reqid = None
    for l in request_text.split('\n'):
        if "<PCT-Waiting_reqid>" in l:
            reqid = re.sub("</?PCT-Waiting_reqid>", "", l).strip()
            break
    return reqid

def check_PCT_status(request_text):
    for l in request_text.split('\n'):
        if "<PCT-Status value=" in l:
            try:
                code = re.search('<PCT-Status value="([a-z]{4,})"/>', l).group(1)
            except AttributeError:
                code = "none"
            break
    return True if code == "success" else False

def get_PCT_smiles(request_text):
    smiles = None
    for l in request_text.split('\n'):
        if "<PCT-Structure_structure_string>" in l:
            smiles = re.sub("</?PCT-Structure_structure_string>", "", l).strip().replace("&#xa;", "")
            break
    return smiles

def poll_PCT(request_id):
    return requests.post(url=urlSend, data=pubchem_poll_string.format(request_id))

def is_smiles(smiles):
    return bool(Chem.MolFromSmiles(smiles)) and bool(smiles)