from typing import List, Tuple

import pandas as pd
import requests
#from rcsbsearchapi.search import AttributeQuery

def GetPDBdata(pdbID_list: List[str], savePATH: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fetch PDB structure validation data via RCSB GraphQL API.

    Args:
        pdbID_list: List of PDB IDs to query (e.g., ['1a3i', '2a3i']).
        savePATH: Directory path to save lipid info CSV.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: A tuple of three dataframes:
            - df_date: Release dates indexed by PDB ID.
            - df_protein: Protein clash validation scores.
            - df_lipid: Ligand (nonpolymer entity) validation data including
              completeness, clashes, stereochemistry, and covalent bond info.
    """
    COLUMNSLIST = ['BioDolphinID', 'Release_Date', 'Covalent', 'Protein_Clash', 'Lipid_Clash', 'Lipid_Stereo']
    
    query = GetQuery(pdbID_list)
    url = f'https://data.rcsb.org/graphql'
    response = requests.post(url, json={'query': query})

    # create dataframe for protein info
    COLUMNSLIST = ['Protein_rcsb_id', 'Protein_Num_Clash'] #rcsb_id: format is in: PDBID.asymChainID
    df_protein = pd.DataFrame(columns=COLUMNSLIST)

    # create dataframe for lipid info
    COLUMNSLIST = ['PDBID', 'Ligand_resNum', 'Ligand_Chain', 'Ligand_CompID', 'Ligand_Covalent', 'Ligand_Num_Clash', 'Ligand_Stereo', 'Ligand_Completeness']
    df_lipid = pd.DataFrame(columns=COLUMNSLIST)

    # create dataframe for PDB release date
    COLUMNSLIST = ['PDBID', 'PDB_release_date']
    df_date = pd.DataFrame(columns=COLUMNSLIST)


    if response.status_code == 200:
        data = response.json()
        for entry in data["data"]["entries"]:
            
            PDB_ID = entry['rcsb_id'].lower()
            date = entry['rcsb_accession_info']['initial_release_date']

            df_new_date = pd.DataFrame({'PDBID': [PDB_ID], 'PDB_release_date': [date]})
            df_date = pd.concat([df_date, df_new_date], ignore_index=True)
            # parse nonpolymer information
            nonpolymer = entry.get('nonpolymer_entities') or []
            
            for entity in nonpolymer:
                for np in entity['nonpolymer_entity_instances']:
                    #print(np)
                    # Get identity information
                    #print('---------------------np-------------------------------')
                    IDs = np['rcsb_nonpolymer_entity_instance_container_identifiers']
                    resNum, ligandChain, compID = IDs['auth_seq_id'], IDs['auth_asym_id'], IDs['comp_id']


                    # Get validation scores information
                    try:
                        validation_scores = np['rcsb_nonpolymer_instance_validation_score'][0]
                        lignad_completeness, ligand_clash, ligand_stereo_outliers = validation_scores['completeness'], validation_scores['intermolecular_clashes'], validation_scores['stereo_outliers']

                    except:
                        pass

                    # Get structure connection information

                    try:
                        connections = np['rcsb_nonpolymer_struct_conn']
                        ligand_covalent = False
                        for conn in connections:
                            connect_type, connect_partner = conn['connect_type'], conn['connect_partner']
                            if connect_type == "covalent bond":
                                ligand_covalent = True
                    except:
                        ligand_covalent = None

                    df_new_lipid = pd.DataFrame({'PDBID':[PDB_ID], 'Ligand_resNum':[resNum], 'Ligand_Chain':[ligandChain], 'Ligand_CompID':[compID], 'Ligand_Covalent':[ligand_covalent], 'Ligand_Num_Clash':[ligand_clash], 'Ligand_Stereo':[ligand_stereo_outliers], 'Ligand_Completeness':[lignad_completeness]})
                    df_lipid = pd.concat([df_lipid, df_new_lipid], ignore_index=True)
        df_lipid.to_csv(f'{savePATH}/PDB_lipid_info.csv')
        return df_date, df_protein, df_lipid
    else:
        print(f'unable to fetch data from PDB. error code is {response.status_code}')
            
            
        



def GetQuery(pdbID_list: List[str]) -> str:
    """Construct RCSB GraphQL query string for PDB structures.

    Args:
        pdbID_list: List of PDB IDs to include in the query.

    Returns:
        str: GraphQL query string ready for POST request to RCSB API.
    """
                
    pdb_ids = "["
    for iL, itm in enumerate(pdbID_list):
        if iL == len(pdbID_list)-1:
            pdb_ids += '"%s"'%itm
        else:
            pdb_ids += '"%s",'%itm
    pdb_ids += "]"
    
    
    query = """{
    entries(entry_ids: """ + pdb_ids +""") {
        rcsb_id
        rcsb_accession_info {
        initial_release_date
        }
        polymer_entities {
        polymer_entity_instances {
            rcsb_id
            rcsb_polymer_instance_feature_summary {
            count
            coverage
            type
            }
        }
        }
        nonpolymer_entities {
        nonpolymer_entity_instances {
            rcsb_nonpolymer_entity_instance_container_identifiers {
            auth_seq_id
            auth_asym_id
            asym_id
            comp_id
            entity_id
            entry_id
            }
            rcsb_nonpolymer_instance_validation_score {
            completeness
            intermolecular_clashes
            stereo_outliers

            }
            rcsb_nonpolymer_struct_conn {
            role
            connect_type
            connect_partner {
                label_asym_id
            }
            connect_target {
                label_asym_id
            }
            }
        }
        }
    }
    }"""

    return query
