import pandas as pd
import re
import yaml
import sdrf_pipelines
import os.path
from sdrf_pipelines.zooma.zooma import OlsClient
from sdrf_pipelines.openms.unimod import UnimodDatabase
from sdrf_pipelines.sdrf.sdrf import SdrfDataFrame

## Accessing ontologies and CVs
unimod = UnimodDatabase()
olsclient = OlsClient()
#print(ols_out)

# modifications have the same column name, not working with pandas
# therefore separated
mod_columns = pd.DataFrame()

# variable and fixed modifications are joined and added as one (each) parameter in the parameter yaml file
fixed_mods = []
variable_mods = []

with open(r'param2sdrf.yml') as file:
   param_mapping = yaml.safe_load(file)
   mapping = param_mapping["parameters"]

# Output dictionary that will become a yaml file
out_yaml = dict()
   
## WE NEED AN SDRF FILE FOR THE EXPERIMENTAL DESIGN, CONTAINING FILE LOCATIONS
sdrf_content = pd.DataFrame()
has_sdrf = os.path.isfile("./sdrf.tsv")
if has_sdrf :
   sdrf_content = pd.read_csv("sdrf.tsv", sep="\t")
   mod_columns = sdrf_content.filter(like="comment[modification parameters]")
   sdrf_content = sdrf_content.drop(columns=mod_columns.columns)
else:
   ## THROW ERROR FOR MISSING SDRF
   exit("ERROR: No SDRF file given. Add an at least minimal version\nFor more details, see https://github.com/bigbio/proteomics-metadata-standard/tree/master/sdrf-proteomics")


## FIRST STANDARD PARAMETERS
# FOR GIVEN PARAMETERS
# CHECK WHETHER COLUMN IN SDRF TO PUT WARNING AND OVERWRITE
# IF NOT GIVEN, WRITE COLUMN
for p in mapping:
   pname = p["name"]
   ptype = p["type"]
   psdrf = "comment[" + p["sdrf"] + "]"
   print("---- Reading parameter: " + pname + " ----")

   if (psdrf) in sdrf_content.keys() :
     pvalue = set(sdrf_content[psdrf])
     if (len(set(sdrf_content[psdrf])) > 1) :
         exit("ERROR: multiple values for parameter " + pname + " in sdrf file\n We recommend separating the sdrf file into parts with the same data analysis parameters")
     # for each type: check consistency
     #print(type(pvalue))
     pvalue = pvalue.pop()
     if ptype == "boolean":
         if not isinstance(pvalue, bool) :
           exit("ERROR: " + pname + " needs to be either \"true\" or \"false\"!!")
     elif ptype == "str":
         if not isinstance(pvalue, str) :
           exit("ERROR: " + pname + " needs to be a string!!")
     elif ptype == "integer":
         if not isinstance(pvalue, int) :
           exit("ERROR: " + pname + " needs to be a string!!")
     elif ptype == "float":
         if not isinstance(pvalue, (float, int)):
           exit("ERROR: " + pname + " needs to be a numeric value!!")
     elif ptype == "class":
         not_matching = [x for x in pvalue.split(",") if x not in p["value"]]
         if not_matching != [] :
            exit("ERROR: " + pname + " needs to have one of these values: " + ' '.join(p["value"]) + "!!\n" + ' '.join(not_matching) + " did not match")



     # Mass tolerances: do they include Da or ppm exclusively?
     if pname == "fragment_mass_tolerance" or pname == "precursor_mass_tolerance":
        print(pvalue)
        unit = pvalue.split(" ")[1]
        pvalue = pvalue.split(" ")[0]
        if unit != "Da" and unit != "ppm" :
          exit("ERROR: " + pname + " allows only units of \"Da\" and \"ppm\", separated by space from the value!!\nWe found " + unit)


     ## ENZYME AND MODIFICATIONS: LOOK UP ONTOLOGY VALUES
     elif pname == "enzyme":
       pvalue = pvalue.split(";")
       pvalue =  [s for s in pvalue if "NT=" in s][0]
       pvalue = pvalue.replace("NT=","")
       ols_out = olsclient.search(pvalue, ontology="MS", exact=True)
       if ols_out == None or len(ols_out) > 1 :
         exit("ERROR: enzyme " + pvalue + " not found in the MS ontology, see https://bioportal.bioontology.org/ontologies/MS/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMS_1001045 for available terms")

     ## Now finally writing the value
     print("Wrote " + str(pvalue))
     out_yaml[pname] = pvalue


## Modifications: look up in Unimod
for m in mod_columns.columns: 
        pvalue = set(mod_columns[m])
        if len(pvalue) > 1:
           exit("ERROR: multiple PTMs in the same column " + m)
        pvalue = pvalue.pop()
        pvalue = pvalue.split(";")
        mmod =  [s for s in pvalue if "NT=" in s][0]
        mmod = mmod.replace("NT=","")
        mresidue =  [s for s in pvalue if "PP=" in s or "TA=" in s][0]
        mresidue = mresidue.replace("PP=","").replace("TA=","")
        mtype =  [s for s in pvalue if "MT=" in s][0]
        mtype = mtype.replace("MT=","")
        print("Found modification " + mmod + " at " + mresidue + " with type " + mtype)

        pvalue = mmod + " of " + mresidue

        if mtype == "fixed":
           fixed_mods.append(pvalue)
        else:
           variable_mods.append(pvalue)
           
## summarizing modifications
out_yaml["fixed_mods"] = ",".join(fixed_mods)
out_yaml["variable_mods"] = ",".join(variable_mods)



print("--- Writing sdrf file into params.yml ---")
with open('params.yml', 'w') as outfile:
    yaml.dump(out_yaml, outfile, default_flow_style=False)
    
