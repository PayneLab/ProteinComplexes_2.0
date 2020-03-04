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

def permutation_test_means(data, group_col, value_col, num_permutations):
    """Use permutation testing to calculate a P value for the difference between the means of two groups."""
    pass

def permutation_test_corr(data, num_permutations):
    """Use permutation testing to calculate a P value for the linear correlation coefficient between two variables in several samples.

    Parameters:
    data (pandas DataFrame): A dataframe where the rows are samples, and the columns are the two variables we're testing correlation between.

    Returns:        
    float: The linear correlation coefficient for the two variables.
    float: The P value for the coefficient.
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
        if perm_coef >= actual_coef:
            extreme_count += 1

    # Calculate the P value
    P_val = extreme_count / num_permutations * 2

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

    # Extract the values
    var1 = data.iloc[:, 0]
    var2 = data.iloc[:, 1]

    # Get the names
    var1_name = data.columns[0]
    var2_name = data.columns[1]

    # Create the plot
    sns.set(style="darkgrid")
    plot = sns.regplot(x = var1, y = var2)
    plot.set(xlabel = var1_name, ylabel = var2_name, title = f"{var1_name} vs. {var2_name}")
    plt.show()


# Data loading functions

def get_member_proteins(complex_name):
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
    self._member_proteins = member_proteins

def get_ubiquitination():
    pass
