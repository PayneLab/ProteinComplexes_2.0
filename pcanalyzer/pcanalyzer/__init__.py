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
    float: The difference between the means
    float: The P value for the null hypothesis that the two groups have the same mean
    list of float: The generated null distribution for the difference between the means
    """
    # If input was a dataframe, check the shape
    if isinstance(data, pd.DataFrame):
        if data.shape[1] != 1:
            raise ValueError(f"Input was a dataframe, so expected 1 column. Found {data.shape[1]}:\n{data}.")

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
    actual_diff = np.mean(tumor) - np.mean(normal)
    abs_actual_diff = abs(actual_diff)

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
        if perm_diff >= abs_actual_diff: # We compare the absolute values for a two-tailed test
            extreme_count += 1

    # Calculate the P value
    P_val = extreme_count / num_permutations # Don't need to multiply by 2 because we compared the absolute values of difference between means.

    return actual_diff, P_val, null_dist

def perm_test_omics_pancancer(datasets, id_list, data_type, num_permutations, paired = False):
    """Use permutation testing to calculate a P value for the difference between the means of tumor and normal sample omics for multiple datasets.

    Parameters:
    datasets (list of cptac.DataSet): A list of instantiations of the datasets we want to do the testing for
    id_list (list of str): A list of the IDs (e.g. protein names) to select for the test
    data_type (str): The type of omics data to do the testing on. Currently supported: proteomics
    num_permutations (int): The number of permutations to perform
    paired (bool, optional): Whether to do a paired test. Default is False.

    Returns:
    dict of str: pandas.DataFrame: A dictionary where:
        the keys are the names of the datasets analyzed, and
        the values are dataframes where
            the index contains the IDs for the data arrays that were tested
            the first column is the difference between the means of tumor and normal for that array, and
            the second column is the raw P value for the difference between the means of the tumor and normal groups in that array.
        The P values are not adjusted. Make sure that you do multiple testing correction.
    int: The number of tests performed, to help with multiple testing correction calculations.
    """
    results = {}
    num_tests = 0

    for dataset in datasets:
        
        ids = []
        diffs = []
        P_vals = []

        if data_type == "proteomics":
            omics = dataset.get_proteomics()
        else:
            raise ValueError(f"Unsupported data type {data_type}. See docstring for supported types.")

        included_cols = (omics.columns.get_level_values(0) & id_list).drop_duplicates()
        selected_omics = omics[included_cols]
        
        for id in selected_omics.columns:
            
            data = selected_omics[id]
            diff, P_val, null_dist = permutation_test_means(data, num_permutations)
            num_tests += 1

            ids.append(id)
            diffs.append(diff)
            P_vals.append(P_val)

        dataset_results = pd.DataFrame({"diff": diffs, "P_val": P_vals}, index=ids)
        results[dataset.get_cancer_type()] = dataset_results

    return results, num_tests

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

def linear_regplot(data):
    """Plot two variables against each other and show the linear regression line.

    Parameters:
    data (pandas DataFrame): A dataframe where the rows are samples, and the columns are the two variables we're testing correlation between.

    Returns:        
    None
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
    data.insert(2, "Sample_Tumor_Normal", sample_type_col)

    # Create the plot
    sns.set(style="darkgrid")
    plot = sns.lmplot(x=var1, y=var2, hue="Sample_Tumor_Normal", data=data, fit_reg=False) 
    sns.regplot(x=var1, y=var2, data=data, scatter=False)
    plot.set(xlabel=var1, ylabel=var2, title=f"{var1} vs. {var2}")
    plt.show()

def boxplot_omics_tumor_normal(dataset, data_type, ids=None, separated_ids=None, dims=(12, 8)):
    """Plot omics data boxplots for tumor vs. normal.

    Parameters:
    dataset (cptac.DataSet): The dataset object to plot the data for
    data_type (str): The type of omics data to plot. Currently supported: proteomics
    ids (list of str, optional): A list of the identifiers (e.g. proteins) to plot data for. Either ids or separated_ids parameter, but not both, must be provided.
    separated_ids (list of list of str, optional): List of groups of identifiers to plot data for. Groups will be plotted with whitespace between groups. Supports exactly two groups. Either ids or separated_ids parameter, but not both, must be provided.
    dims (tuple of int, int, optional): A tuple containing the dimensions for the plot, passed to figsize argument in plt.figure(). Defaults to (12, 8)

    Returns:
    None
    """
    if separated_ids is not None:
        if ids is not None:
            raise ValueError("ids and separated_ids parameters cannot both not be None.")
        if len(separated_ids) > 2:
            raise ValueError(f"boxplot_omics_tumor_normal does not support more than two separated groups. You passed {len(separated_ids)}:\n{separated_ids}")

        ids = separated_ids[0] + separated_ids[1]
    elif ids is not None:
        pass
    else:
        raise ValueError("ids and separated_ids cannot both be not None.")

    if data_type == "proteomics":
        omics = dataset.get_proteomics()
    else:
        raise ValueError(f"Unsupported data type {data_type}. See docstring for supported types.")

    included_cols = (pd.Series(ids)[pd.Series(ids).isin(omics.columns)]).drop_duplicates()
    selected_omics = omics[included_cols] # Note that if the ids list was constructed from separated_ids lists, we will have selected the columns in the desired order

    # Flatten the columns if needed
    if selected_omics.columns.nlevels > 1:
        #selected_omics.columns = flatten_idx_differentiate_duplicates(selected_omics.columns)
        tuples = selected_omics.columns.to_flat_index() # Converts multiindex to an index of tuples
        no_nan = tuples.map(lambda x: [item for item in x if pd.notnull(item)]) # Cut any NaNs out of tuples
        joined = no_nan.map(lambda x: '_'.join(x)) # Join each tuple
        selected_omics.columns = joined

    if separated_ids is not None:
        selected_omics.insert(len(separated_ids[0]), "", np.nan) # Plotting the column of nothing between the two groups will separate them visually

    # Grab the order of the columns so we can pass it to order them that way in the plot
    plot_order = selected_omics.columns

    # Create a column that marks sample type
    sample_type_col = np.where(selected_omics.index.str.endswith(".N"), "Normal", "Tumor")
    selected_omics = selected_omics.assign(**{"Sample_Tumor_Normal": sample_type_col})

    # Unpivot the selected_omicsframe (similar to a tidyr gather) so it works with Seaborn
    selected_omics = pd.melt(selected_omics, id_vars="Sample_Tumor_Normal", var_name="id", value_name="omics")

    # Create the plot
    sns.set(style="darkgrid")

    plt.figure(figsize=dims)
    boxplt = sns.boxplot(data=selected_omics, x='id', y='omics', hue='Sample_Tumor_Normal', order=plot_order, showfliers=False)
    boxplt.set_title(dataset.get_cancer_type())
    boxplt.set_ylabel('omics')
    boxplt.set_xlabel('')
    boxplt.set_xticklabels(boxplt.get_xticklabels(), rotation=45)
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

# Helper functions

def flatten_idx_differentiate_duplicates(idx):
    """Flattens a pandas.MultiIndex, dropping level values unless they're needed to differentiate duplicates, in which case it concatenates level values into a string with underscore separators.

    Parameters:
    idx (pandas.MultiIndex): The index to flatten.

    Returns:
    pandas.Index: The flattened index.
    """
