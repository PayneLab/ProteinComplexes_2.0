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
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Stats functions

def permutation_test_means(data, num_permutations, paired=False):
    """Use permutation testing to calculate a P value for the difference between the means of two groups.

    Parameters:
    data (pandas Series or Dataframe): A series or single column dataframe containing all the data values, with a patient ID index (which indicates the tumor/normal grouping).
    num_permutations (int): The number of permutations to perform
    paired (bool, optional): Whether to do a paired test. Default is False.

    Returns:
    float: The P value for the null hypothesis that the two groups have the same mean.
    """
    # If input was a dataframe, check the shape
    if isinstance(data, pd.DataFrame):
        if data.shape[1] != 1:
            raise ValueError(f"Input was a dataframe, so expected 1 column. Found {data.shape[1]}.")

    # If the input was a series, make it a dataframe
    if isinstance(data, pd.Series):
        data = data.to_frame()

    # Drop NaN values
    data = data.dropna()

    # Split into tumor/normal and extract values
    tumor_selector = ~data.index.str.endswith(".N")
    normal_selector = data.index.str.endswith(".N")
    tumor = data[tumor_selector].iloc[:, 0]
    normal = data[normal_selector].iloc[:, 0]

    # Create an independent pseudo-random number generator
    generator = np.random.RandomState(0)

    # Calculate the actual correlation coefficient
    actual_diff = abs(np.mean(tumor) - np.mean(normal))

    null_dist = []
    extreme_count = 0

    for i in range(num_permutations):
        # Permute values
        perm_array = generator.permutation(data.iloc[:, 0])

        # Split into tumor/normal and extract values
        perm_tumor = perm_array[tumor_selector]
        perm_normal = perm_array[normal_selector]

        # Calculate the actual correlation coefficient
        perm_diff = abs(np.mean(perm_tumor) - np.mean(perm_normal))

        # Add it to our null distribution
        null_dist.append(perm_diff)

        # Keep count of how many are as or more extreme than our coefficient
        if perm_diff >= actual_diff: # We compare the absolute values for a two-tailed test
            extreme_count += 1

    # Calculate the P value
    P_val = extreme_count / num_permutations # Don't need to multiply by 2 because we compared the absolute values of difference between means.

    return actual_diff, P_val, null_dist

def permutation_test_corr(data, num_permutations):
    """Use permutation testing to calculate a P value for the linear correlation coefficient between two variables in several samples.

    Parameters:
    data (pandas DataFrame): A dataframe where the rows are samples, and the columns are the two variables we're testing correlation between.

    Returns:        
    float: The linear correlation coefficient for the two variables.
    float: The P value for the null hypothesis that the correlation coefficient is zero.
    """

    # Check the table dimensions
    if data.shape[1] != 2:
        raise ValueError(f"Expected 2 columns in dataframe. Found {data.shape[1]}.")

    # Drop NaN values
    data = data.dropna()

    # Extract the values
    var1 = data.iloc[:, 0].values
    var2 = data.iloc[:, 1].values

    # Create an independent pseudo-random number generator
    generator = np.random.RandomState(0)

    # Calculate the actual correlation coefficient
    actual_coef = np.corrcoef(var1, var2)[0, 1]

    null_dist = []
    extreme_count = 0

    for i in range(num_permutations):
        var1_perm = generator.permutation(var1)
        perm_coef = np.corrcoef(var1_perm, var2)[0, 1]

        # Add it to our null distribution
        null_dist.append(perm_coef)

        # Keep count of how many are as or more extreme than our coefficient
        if abs(perm_coef) >= abs(actual_coef): # We compare the absolute values for a two-tailed test
            extreme_count += 1

    # Calculate the P value
    P_val = extreme_count / num_permutations # Don't need to multiply by 2 because we compared the absolute values of coefficients.

    return actual_coef, P_val, null_dist

# Plotting functions

def plot_linear_regression(data):
    """Plot two variables against each other and show the linear regression line.

    Parameters:
    data (pandas DataFrame): A dataframe where the rows are samples, and the columns are the two variables we're testing correlation between.

    Returns:        

    """
    # Check the table dimensions
    if data.shape[1] != 2:
        raise ValueError(f"Expected 2 columns in dataframe. Found {data.shape[1]}.")

    # Drop NaN values
    data = data.dropna()

    # Flatten the columns if needed.
    if data.columns.nlevels > 1:
        tuples = data.columns.to_flat_index() # Converts multiindex to an index of tuples
        no_nan = tuples.map(lambda x: [item for item in x if pd.notnull(item)]) # Cut any NaNs out of tuples
        joined = no_nan.map(lambda x: '_'.join(x)) # Join each tuple
        data.columns = joined

    # Get the names
    var1 = data.columns[0]
    var2 = data.columns[1]

    # Create a column that marks sample type
    sample_type_col = np.where(data.index.str.endswith(".N"), "Normal", "Tumor")
    data.insert(2, "Sample_Type", sample_type_col)

    # Create the plot
    sns.set(style="darkgrid")
    plot = sns.lmplot(x=var1, y=var2, hue="Sample_Type", data=data, fit_reg=False) 
    sns.regplot(x=var1, y=var2, data=data, scatter=False)
    plot.set(xlabel=var1, ylabel=var2, title=f"{var1} vs. {var2}")
    plt.show()


# Data loading functions

def get_corum_protein_lists():
    """Reads file from CORUM and returns a dictionary where the keys are protein complex names, and the values are lists of proteins that are members of those complexes. Data downloaded from CORUM 3.0 (released 2018, http://mips.helmholtz-muenchen.de/corum/#download)"""
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_files_path = os.path.join(path_here, "data_files")

    member_proteins = pd.read_csv(os.path.join(data_files_path, 'corum_protein_complexes.tsv'), sep='\t')
    member_proteins = member_proteins.loc[member_proteins['Organism'] == 'Human']
    member_proteins = member_proteins.set_index("ComplexName")

    # Select the member proteins column and split the proteins in each complex into values of a list
    member_proteins = member_proteins['subunits(Gene name)'].str.split(';')

    # For complexes with multiple entries (due to different samples), combine the lists
    member_proteins = member_proteins.groupby(member_proteins.index).agg(sum) # Sum will concatenate lists
    member_proteins = member_proteins.apply(set).apply(sorted) # Get rid of duplicates by converting to set. Then go back to list.
    member_proteins = member_proteins.to_dict()

    return member_proteins

def get_hgnc_protein_lists():
    """Reads file from HGNC gene family dataset and returns a dictionary where the keys are protein complex names, and the values are lists of proteins that are members of those complexes. Data downloaded from HGNC BioMart server on 12 Mar 2020 (https://biomart.genenames.org/)."""
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_files_path = os.path.join(path_here, "data_files")

    member_proteins = pd.read_csv(os.path.join(data_files_path, 'hgnc_protein_families.tsv'), sep='\t')
    member_proteins = member_proteins.set_index("Family name")
    member_proteins = member_proteins["Approved symbol"]

    # Combine multiple rows per family into one row with a list
    member_proteins = member_proteins.groupby(member_proteins.index).agg(list)
    member_proteins = member_proteins.apply(set).apply(sorted) # Get rid of duplicates by converting to set. Then go back to list.
    member_proteins = member_proteins.to_dict()

    return member_proteins



def get_ubiquitination():
    pass
