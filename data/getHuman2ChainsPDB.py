import requests
import os
import sys
from Bio.PDB import PDBList, PDBParser, PDBIO, Select


OUTPUT_DIR = "Human_2Chains_pdbs"
DOWNLOAD_LIMIT = None


class ChainSelect(Select):
    def accept_residue(self, residue):
        return residue.get_id()[0] == " "

def normalize_chains(structure):

    model = structure[0]
    chains = list(model.get_chains())
    

    if len(chains) != 2:
        return False
    
    chain_ids = {c.id for c in chains}
    

    if chain_ids == {'A', 'B'}:
        return True
        

    elif chain_ids == {'H', 'L'}:
        
        for c in chains:
            if c.id == 'H': c.id = 'temp_A'
            if c.id == 'L': c.id = 'temp_B'
        
        for c in chains:
            if c.id == 'temp_A': c.id = 'A'
            if c.id == 'temp_B': c.id = 'B'
        return True
        

    else:

        chains.sort(key=lambda c: len(list(c.get_residues())), reverse=True)
        chains[0].id = 'A'
        chains[1].id = 'B'
        return True

def get_human_dimers():
    print("Querying RCSB PDB for Homo sapiens (Strict 2 chains)...")
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entity_source_organism.taxonomy_lineage.name", "operator": "exact_match", "value": "Homo sapiens"}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.polymer_entity_count_protein", "operator": "equals", "value": 2}}
            ]
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    
    try:
        r = requests.post(url, json=query)
        r.raise_for_status()
        ids = [item['identifier'] for item in r.json().get('result_set', [])]
        print(f"Found {len(ids)} structures.")
        return ids
    except Exception as e:
        print(f"Search failed: {e}")
        return []

def process_pdbs(pdb_ids):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pdbl = PDBList(verbose=False)
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    
    targets = pdb_ids[:DOWNLOAD_LIMIT] if DOWNLOAD_LIMIT else pdb_ids
    
    success_count = 0
    
    for i, pdb_id in enumerate(targets):
        sys.stdout.write(f"\rProcessing [{i+1}/{len(targets)}]: {pdb_id} ... ")
        sys.stdout.flush()
        
        try:

            raw_file = pdbl.retrieve_pdb_file(pdb_id, pdir=OUTPUT_DIR, file_format="pdb")
            if not os.path.exists(raw_file): continue
            

            structure = parser.get_structure(pdb_id, raw_file)

            is_valid = normalize_chains(structure)
            
            if is_valid:

                clean_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_AB.pdb")
                io.set_structure(structure)
                io.save(clean_path, select=ChainSelect())
                success_count += 1

            os.remove(raw_file)
            
        except Exception as e:

            continue
            
    print(f"\n\nDone! Saved {success_count} normalized files to '{OUTPUT_DIR}/'")

if __name__ == "__main__":
    ids = get_human_dimers()
    if ids:
        process_pdbs(ids)
