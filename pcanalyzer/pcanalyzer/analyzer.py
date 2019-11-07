#   Copyright 2018 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import os
import warnings
import pandas as pd
import cptac

class Analyzer:
    """This class loads CPTAC datasets and provides functions for studying the behavior of protein complexes within them."""

    def __init__(self, load="all"):
        """Load datasets and protein complex lists.

        Parameters:
        load (str or list of str): The name(s) of the datasets you want to load, in the cptac package format. Default "all" loads all datasets.
        """

        # Initialize the datasets variable
        self._datasets = {}        

        # Parse the datasets parameter and load the specified datasets
        dataset_names_to_inits = {
            "brca": cptac.Brca,
            "ccrcc": cptac.Ccrcc,
            "colon": cptac.Colon,
            "endometrial": cptac.Endometrial,
            "gbm": cptac.Gbm,
            "hnscc": cptac.Hnscc,
            "luad": cptac.Luad,
            "ovarian": cptac.Ovarian,
        }

        if load == "all": # Load all datasets
            for name, init in dataset_names_to_inits.items():
                self._datasets[name] = init()

        else: # They only want specific datasets
            if isinstance(load, str): # If it's a single string other than "all", make it a list of length one
                load = [load]
                
            for name in load: # Load the datasets corresponding to the strings they passed, and save in our self._datasets dictionary
                if name in dataset_names_to_inits.keys():
                    self._datasets[name] = dataset_names_to_inits[name]()
                else:
                    warnings.warn(f"Requested dataset '{name}' is invalid. Dataset not loaded.")

        # Create a dictionary with protein complex names as keys and lists of included proteins as values
        # Data downloaded from CORUM 3.0 (released 2018, http://mips.helmholtz-muenchen.de/corum/#download) 
        member_proteins = pd.read_csv(os.path.join(data_files_path, 'allComplexes.txt'), sep='\t')
        member_proteins = member_proteins.loc[member_proteins['Organism'] == 'Human']
        member_proteins = member_proteins.set_index("ComplexName")

        # Select the member proteins column and split the proteins in each complex into values of a list
        member_proteins = member_proteins['subunits(Gene name)'].str.split(';')

        # For complexes with multiple entries (due to different samples), combine the lists
        member_proteins = member_proteins.groupby(member_proteins.index).agg(sum) # Sum will concatenate lists
        member_proteins = member_proteins.apply(set).apply(sorted) # Get rid of duplicates by converting to set. Then go back to list.
        self.mp = member_proteins
