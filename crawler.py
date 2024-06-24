import requests
import xml.etree.ElementTree as ET
import pandas as pd
import os

class PDBNotFoundException(Exception):
    pass
class DownloadNotSuccessfulException(Exception):
    pass    

def get_pdb_ids(uniprot_id):
    # UniProt API endpoint for retrieving entry in XML format
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    
    if response.status_code != 200:
        raise Exception(f"Failed to retrieve data for UniProt ID {uniprot_id}")

    # Parse the XML response
    root = ET.fromstring(response.content)
    pdb_ids = []
    
    # Find all PDB cross-references
    for ref in root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='PDB']"):
        pdb_ids.append(ref.attrib['id'])
    
    
    return pdb_ids

def download_pdb(pdb_id, output_file):
    # RCSB PDB API endpoint for downloading PDB files
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    
    if response.status_code != 200:
        raise DownloadNotSuccessfulException(f"Failed to download PDB file for PDB ID {pdb_id}")

    # Write the PDB file content to a local file
    with open(output_file, 'wb') as file:
        file.write(response.content)

def get_pdb(uniprot_id):
    pdb_ids = get_pdb_ids(uniprot_id)
    path = os.getcwd() + "/tmp/"
    if pdb_ids:
        for pdb_id in pdb_ids:
            try:
                output_file = f"{path}{uniprot_id}.pdb"
                download_pdb(pdb_id, output_file)
                print(f"PDB file {output_file} downloaded successfully.")
                break
            except DownloadNotSuccessfulException:
                continue
    else:

        raise PDBNotFoundException(f"No PDB IDs found for UniProt ID {uniprot_id}")

def download_alphafold_structure(uniprot_id):
    path = os.getcwd() + "/tmp/"
    output_file = f"{path}{uniprot_id}.pdb"
    # AlphaFold Protein Structure Database URL template
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    
    response = requests.get(url)
    
    if response.status_code != 200:
        raise Exception(f"Failed to download AlphaFold structure for UniProt ID {uniprot_id}")

    # Write the PDB file content to a local file
    with open(output_file, 'wb') as file:
        file.write(response.content)


if __name__ == "__main__":
    
    train_data = pd.read_csv("./scaffold_split/protein_train.csv")
    train_protein_ID = set(train_data["protein_ID"].values.tolist())
    val_data = pd.read_csv("./scaffold_split/protein_semi_inductive.csv")
    val_protein_ID = set(val_data["protein_ID"].values.tolist())
    test_data = pd.read_csv("./scaffold_split/protein_transductive.csv")
    test_protein_ID = set(test_data["protein_ID"].values.tolist())
    for protein_id in train_protein_ID | val_protein_ID | test_protein_ID:
        try:
            path = os.getcwd() + "/tmp/" + f"{protein_id}.pdb"
            if not os.path.isfile(path):
                get_pdb(protein_id)
            else:
                print(f"There exists {protein_id}. Skip!")
        except PDBNotFoundException:
            download_alphafold_structure(protein_id)