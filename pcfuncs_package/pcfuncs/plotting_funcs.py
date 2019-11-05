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

from .analysis_funcs import *

import contextlib

"""
plot_ratios
-----------
Create a seaborn plot for the ratios in a dataframe

Parameters:
    prot1, prot2 = the two proteins to plot (ratio of prot1 / prot2)
    by_patient (optional) = T/F whether or not to visualize matched samples
    mutation_list (optional) = plot will distinguish patients with a mutation in any of the proteins in this list
    
Returns:
    Displays a plot of the data
    
"""
def plot_ratios(prot1, prot2, by_patient = False, mutation_list = None):
    
    # Set up dataframe for plotting
    tumor_ratio_df, normal_ratio_df = get_ratio_df(prot1, prot2)
    plot_data = tumor_ratio_df.append(normal_ratio_df)
    
    # Reformat for visualizing matched samples if necessary
    if by_patient or mutation_list:
        # Create a new column for matched status
        tumor_ratio_df['Matched_Status'] = 'Unmatched'
        normal_ratio_df['Matched_Status'] = 'Unmatched'
        # Classify samples as matched/unmatched
        normal_ratio_df.loc[((normal_ratio_df['Patient_ID']).str.replace('_NM', '')).isin(tumor_ratio_df['Patient_ID']), 'Matched_Status'] = 'Matched'
        tumor_ratio_df.loc[((tumor_ratio_df['Patient_ID']) + '_NM').isin(normal_ratio_df['Patient_ID']), 'Matched_Status'] = 'Matched'
        # Label samples accordingly in plot_data
        plot_data = tumor_ratio_df.append(normal_ratio_df)
        if mutation_list:
            # Find patients with mutations in the given proteins
            patient_list = find_unique_mutations(mutation_list)
            plot_data.loc[plot_data['Patient_ID'].isin(patient_list), 'Matched_Status'] = 'Has_Mutation'
        plot_data.loc[plot_data['Matched_Status'] == 'Unmatched', 'Patient_ID'] = 'Unmatched_Sample'
        plot_data['Patient_ID'] = plot_data['Patient_ID'].str.replace('_NM','')
    
    # Reformat for visualizing a certain list of patients if necessary
    if mutation_list and not by_patient:
        plot_data.loc[~plot_data['Patient_ID'].isin(patient_list), 'Patient_ID'] = 'Other Mutation'
        plot_data.loc[plot_data['Patient_ID'].isin(patient_list), 'Patient_ID'] = 'Given Mutation'
    elif mutation_list:
        plot_data.loc[plot_data['Patient_ID'].isin(patient_list), 'Patient_ID'] = 'Given Mutation'
        
    # Print results of statistical tests
    print('T-test p-value: ' + str(scipy.stats.ttest_ind(tumor_ratio_df['Ratio'], normal_ratio_df['Ratio'])[1]))
    print('Levene p-value: ' + str(scipy.stats.levene(tumor_ratio_df['Ratio'], normal_ratio_df['Ratio'])[1]))
    
    a4_dims = (10, 10)
    fig, ax = plt.subplots(figsize=a4_dims)

    # Create the plot
    if by_patient or patient_list:
        boxplt = sns.boxplot(data=plot_data, x='Sample_Type', y='Ratio', color='w', showfliers=False)
        boxplt = sns.stripplot(data=plot_data, x='Sample_Type', y='Ratio', hue='Patient_ID', size=10, dodge=True, jitter=True)
        boxplt.get_legend().set_bbox_to_anchor((1, 1, 0, 0))
    else:
        boxplt = sns.boxplot(data=plot_data, x='Sample_Type', y='Ratio', showfliers=False)
        boxplt = sns.stripplot(data=plot_data, x='Sample_Type', y='Ratio', dodge=True, jitter=True, color='.3')

    # Add styling
    boxplt.set_title('Ratio of ' + prot1 + ' / ' + prot2, fontsize='25')
    boxplt.set_xlabel('')
    boxplt.set_ylabel('Protein Expression Ratio', fontsize='20')
    boxplt.tick_params(labelsize='15')

# https://stackoverflow.com/questions/24685012/pandas-dataframe-renaming-multiple-identically-named-columns
# Rename duplicated columns in a dataframe so that each column has a unique name
def df_column_uniquify(df):
    df_columns = df.columns
    new_columns = []
    for item in df_columns:
        counter = 1
        newitem = item
        while newitem in new_columns:
            counter += 1
            newitem = "{}_{}".format(item, counter)
        new_columns.append(newitem)
    df.columns = new_columns
    return df

def plot_phosphoproteomics(protein, by_patient = False, print_pvals = True, remove_duplicates = False, cancer_type = 'ov'):
    """For a particular cancer type, plot the difference in phosphorylation levels between tumor and normal samples, for each phosphorylation site."""
    
    # Get data from the appropriate cancer type
    if cancer_type == 'ov':
        phos_data = ov.join_metadata_to_omics('clinical', 'phosphoproteomics', 
                            metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                            omics_genes = protein)
        phos_data = ov.reduce_multiindex(phos_data, flatten=True)
    elif cancer_type == 'colon': 
        phos_data = colon.join_metadata_to_omics('clinical', 'phosphoproteomics', 
                            metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                            omics_genes = protein)
        phos_data = colon.reduce_multiindex(phos_data, flatten=True)
    elif cancer_type == 'en': 
        phos_data = en.join_metadata_to_omics('clinical', 'phosphoproteomics', 
                            metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                            omics_genes = protein)
        phos_data = en.reduce_multiindex(phos_data, flatten=True)
    elif cancer_type == 'renal': 
        phos_data = renal.join_metadata_to_omics('clinical', 'phosphoproteomics', 
                            metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                            omics_genes = protein)
        phos_data = renal.reduce_multiindex(phos_data, flatten=True)
    else: print('Error: cancer_type must be "ov", "colon", "renal", or "en"')
    
    if remove_duplicates: phos_data = phos_data.loc[:, ~phos_data.columns.duplicated()]
    else: phos_data = df_column_uniquify(phos_data)
    plot_data = pd.melt(phos_data, id_vars = ['Patient_ID', 'Sample_Tumor_Normal'], var_name = 'Location', value_name = 'Reading').dropna(axis = 0)
    plot_data['Location'] = plot_data['Location'].str.replace('_phosphoproteomics', '')
    
    # Perform t-tests on tumor vs normal phosphoproteomics for each site
    if print_pvals:
        for column in phos_data.columns:
            if column != 'Patient_ID' and column != 'Sample_Tumor_Normal':
                sitedf = phos_data[['Patient_ID', 'Sample_Tumor_Normal', column]]
                tumordf = sitedf.loc[sitedf['Sample_Tumor_Normal'] == 'Tumor'].dropna(axis = 0)
                normaldf = sitedf.loc[sitedf['Sample_Tumor_Normal'] == 'Normal'].dropna(axis = 0)
                if len(tumordf) > 2 and len(normaldf) > 2:
                    pval = scipy.stats.ttest_ind(tumordf[column], normaldf[column])[1]
                    print(column + ' t-test: ' + str(pval))
    
    # Plot the data
    a4_dims = (20, 20)
    fig, ax = plt.subplots(figsize=a4_dims)

    boxplt = sns.boxplot(data=plot_data, x='Location', y='Reading', hue='Sample_Tumor_Normal', showfliers=False, palette = 'Blues')
    boxplt = sns.stripplot(data=plot_data, x='Location', y='Reading', hue='Sample_Tumor_Normal', dodge=True, jitter=True, color='.3')

    # Add styling
    boxplt.set_title(protein + ' Phosphorylation', fontsize='25')
    boxplt.set_xlabel('')
    boxplt.set_ylabel('Phosphorylation Level', fontsize='20')
    boxplt.tick_params(labelsize='10')

"""
plot_proteomics
---------------
Create a seaborn plot of the proteomics data for a single protein

Parameters:
    protein = the protein to plot
    use_cptac (optional) = T/F whether to use the CPTAC proteomics data (default) or the non-normalized data
    by_patient (optional) = T/F whether or not to visualize matched samples
    print_pvals (optional) = T/F whether or not to print the p-values of the statistical tests
    
Returns:
    Displays a plot of the data

"""

def plot_proteomics(protein, use_cptac = True, by_patient = False, print_pvals = True, cancer_type = 'ov', show_subtype = False):
    
    # This uses the normalized proteomics data from the CPTAC package
    if use_cptac:
        # Get data from the appropriate cancer type
        if cancer_type == 'ov':
            plot_data = ov.join_metadata_to_omics('clinical', 'proteomics', 
                                                metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                                                omics_genes = protein)
            plot_data = ov.reduce_multiindex(plot_data, levels_to_drop="Database_ID")
        elif cancer_type == 'colon': 
            plot_data = colon.join_metadata_to_omics('clinical', 'proteomics', 
                                                metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                                                omics_genes = protein)
        elif cancer_type == 'en': 
            plot_data = en.join_metadata_to_omics('clinical', 'proteomics', 
                                                metadata_cols = ['Patient_ID', 'Proteomics_Tumor_Normal'], 
                                                omics_genes = protein)
        elif cancer_type == 'renal':
            plot_data = renal.join_metadata_to_omics('clinical', 'proteomics', 
                                                metadata_cols = ['Patient_ID', 'Sample_Tumor_Normal'], 
                                                omics_genes = protein)
            plot_data = renal.reduce_multiindex(plot_data, levels_to_drop="Database_ID")
        else: 
            print('Error: cancer_type must be "ov", "colon", "en", or "renal"')
            return
            
        plot_data = plot_data.loc[:, ~plot_data.columns.duplicated()]
        if cancer_type == 'en':
            plot_data.rename(columns={protein+'_proteomics': protein, 'Proteomics_Tumor_Normal': 'Sample_Type'}, inplace = True)
            plot_data.loc[plot_data['Sample_Type'] != 'Tumor', 'Sample_Type'] = 'Normal'
        else:
            plot_data.rename(columns={protein+'_proteomics': protein, 'Sample_Tumor_Normal': 'Sample_Type'}, inplace = True)
        plot_data = plot_data.dropna(axis = 0)
        plot_data['Matched_Status'] = 'Unmatched'
        tumor_df = plot_data.loc[plot_data['Sample_Type'] == 'Tumor']
        normal_df = plot_data.loc[plot_data['Sample_Type'] == 'Normal']
    
    # Otherwise use the non-normalized data to make the plot
    else:
        plot_data = pd.DataFrame(data.loc[data.index == protein].transpose())
        plot_data['Sample_Type'] = 'Tumor'
        plot_data.loc[plot_data.index.str.contains('_NM'), 'Sample_Type'] = 'Normal'

        plot_data['Patient_ID'] = plot_data.index
        plot_data['Matched_Status'] = 'Unmatched'
        tumor_df = plot_data.loc[plot_data['Sample_Type'] == 'Tumor']
        normal_df = plot_data.loc[plot_data['Sample_Type'] == 'Normal']
    
    # Format to show matched patients if necessary
    if by_patient:
        pd.options.mode.chained_assignment = None 
        
        # Classify samples as matched/unmatched
        normal_df.loc[((normal_df['Patient_ID']).str.replace('_NM|N', '')).isin(tumor_df['Patient_ID']), 'Matched_Status'] = 'Matched'
        tumor_df.loc[((tumor_df['Patient_ID']) + '_NM').isin(normal_df['Patient_ID']), 'Matched_Status'] = 'Matched'
        tumor_df.loc[((tumor_df['Patient_ID']) + 'N').isin(normal_df['Patient_ID']), 'Matched_Status'] = 'Matched'
        tumor_df.loc[('N' + (tumor_df['Patient_ID'])).isin(normal_df['Patient_ID']), 'Matched_Status'] = 'Matched'
        
        # Label samples accordingly in plot_data
        plot_data = tumor_df.append(normal_df)
        plot_data.loc[plot_data['Matched_Status'] == 'Unmatched', 'Patient_ID'] = 'Unmatched_Sample'
        plot_data['Patient_ID'] = plot_data['Patient_ID'].str.replace('_NM|N','')
        
        if show_subtype:
            plot_data = pd.merge(plot_data, colon.get_clinical()[['Mutation_Phenotype', 'Patient_ID']], on = 'Patient_ID')
             
    if print_pvals:
        # Print results of statistical tests
        print('T-test p-value: ' + str(scipy.stats.ttest_ind(tumor_df[protein], normal_df[protein])[1]))
        print('Levene p-value: ' + str(scipy.stats.levene(tumor_df[protein], normal_df[protein])[1]))
    
    # Create the plot
    a4_dims = (10, 10)
    fig, ax = plt.subplots(figsize=a4_dims)
    
    if by_patient and show_subtype:
        boxplt = sns.boxplot(data=plot_data, x='Sample_Type', y=protein, color='w', showfliers=False)
        boxplt = sns.stripplot(data=plot_data, x='Sample_Type', y=protein, hue = 'Mutation_Phenotype', size=10, dodge=True, jitter=True)
        boxplt.get_legend().set_bbox_to_anchor((1, 1, 0, 0))
    elif by_patient:
        boxplt = sns.boxplot(data=plot_data, x='Sample_Type', y=protein, color='w', showfliers=False)
        boxplt = sns.stripplot(data=plot_data, x='Sample_Type', y=protein, hue='Patient_ID', size=10, dodge=True, jitter=True)
        boxplt.get_legend().set_bbox_to_anchor((1, 1, 0, 0))
    else:
        boxplt = sns.boxplot(data=plot_data, x='Sample_Type', y=protein, showfliers=False)
        boxplt = sns.stripplot(data=plot_data, x='Sample_Type', y=protein, dodge=True, jitter=True, color='.3')
        
    # Add styling
    boxplt.set_title(protein + ' Proteomics', fontsize='25')
    boxplt.set_xlabel('')
    boxplt.set_ylabel('Protein Expression', fontsize='20')
    boxplt.tick_params(labelsize='15')

"""
plot_complex_clinical
---------------------
Creates a countplot of the number of proteins in a complex that are differentially expressed in tumor/normal matched samples
Also displays the given clinical feature 

Parameters:
    cancer_type = can be colon, en (endometrial)
    clinical_feature = the name of the column from the clinical dataframe that you want to visualize
    protein_list = list of proteins to analyze (most likely the proteins in a complex)
    expression_change = increased, decreased, or nochange; will determine if the plot shows the number
        of proteins with increased expression in tumor samples, decreased in tumor samples, or no change between
        tumor/normal
    bin_size = size of the bins to put patients in as number of proteins on the x-axis (default is 1)
    complex_name (optional) = The title for the plot
    show_proportions = whether to display proportions of each group in each category or raw numbers
    
Returns:
    Displays a Seaborn countplot
    x-axis: Number of proteins with differential expression
    y-axis: Number of patients in this category
    hue: clinical feature
    
"""

def plot_complex_clinical(cancer_type, clinical_feature, raw_protein_list, expression_change = 'increased', bin_size = 1, complex_name = '', show_proportions = True, boxplot = False):
    # Remove duplicates from protein list
    raw_protein_list = list(set(raw_protein_list))
    protein_list = []
    
    if cancer_type == 'colon':
        # Make sure all the proteins in protein_list are in the proteomics df
        for protein in raw_protein_list:
            if protein in colon.get_proteomics().columns:
                protein_list.append(protein)
        # Start out with a df including the patient IDs and the given clinical feature
        plot_data = pd.DataFrame()
        plot_data = colon.get_clinical()[['Patient_ID', clinical_feature]]
        plot_data['Sample_ID'] = plot_data.index
        
        # Add protein expression values for all proteins in the given list
        proteomics_data = colon.get_proteomics()[protein_list]
        plot_data = plot_data.join(proteomics_data)

        # Extract all the normal matched samples into a new df and remove the 'N' to match them with tumor samples
        normal = plot_data.copy()
        normal = normal.loc[normal['Patient_ID'].str.contains('N')]
        normal['Patient_ID'] = normal['Patient_ID'].str.replace('N', '')
        normal.index = normal['Patient_ID']

        # Create new column to count the number of proteins with increased/decreased/nochange levels in tumor samples
        # Remove the normal samples from plot_data
        plot_data[expression_change] = 0 
        plot_data = plot_data.loc[~plot_data['Patient_ID'].str.contains('N')]
        plot_data.index = plot_data['Patient_ID']

        # Remove unmatched samples from each dataframe
        plot_data = plot_data.loc[plot_data['Patient_ID'].isin(normal['Patient_ID'])]
        normal = normal.loc[normal['Patient_ID'].isin(plot_data['Patient_ID'])] 
    
    elif cancer_type == 'en':
        # Make sure all the proteins in protein_list are in the proteomics df
        for protein in raw_protein_list:
            if protein in en.get_proteomics().columns:
                protein_list.append(protein)
        
        # Start out with a df including the patient IDs and the given clinical feature
        plot_data = pd.DataFrame()
        plot_data = en.get_clinical()[['Patient_ID', 'Proteomics_Tumor_Normal']]
        plot_data['Sample_ID'] = plot_data.index
        
        # Add protein expression values for all proteins in the given list
        #proteomics_data = en.get_proteomics()[protein_list]
        proteomics_data = en.join_metadata_to_omics('derived_molecular', 'proteomics', clinical_feature, protein_list)
        plot_data = plot_data.join(proteomics_data)
        plot_data.columns = plot_data.columns.str.replace('_proteomics', '')

        # Extract all the normal matched samples into a new df and remove the 'N' to match them with tumor samples
        normal = plot_data.copy()
        normal = normal.loc[normal['Proteomics_Tumor_Normal'].str.contains('Adjacent_normal')]
        normal.index = normal['Patient_ID']

        # Create new column to count the number of proteins with increased/decreased/nochange levels in tumor samples
        # Remove the normal samples from plot_data
        plot_data[expression_change] = 0 
        plot_data = plot_data.loc[plot_data['Proteomics_Tumor_Normal'] == 'Tumor']
        plot_data.index = plot_data['Patient_ID']

        # Remove unmatched samples from each dataframe
        plot_data = plot_data.loc[plot_data['Patient_ID'].isin(normal['Patient_ID'])]
        normal = normal.loc[normal['Patient_ID'].isin(plot_data['Patient_ID'])] 
    
    elif cancer_type == 'renal':
        # Make sure all the proteins in protein_list are in the proteomics df
        for protein in raw_protein_list:
            if protein in renal.get_proteomics().columns:
                protein_list.append(protein)
        
        # Start out with a df including the patient IDs and the given clinical feature
        plot_data = pd.DataFrame()
        plot_data = renal.get_clinical()[['Patient_ID', clinical_feature, 'Sample_Tumor_Normal']]
        plot_data['Sample_ID'] = plot_data.index
        
        # Add protein expression values for all proteins in the given list
        proteomics_data = renal.get_proteomics()[protein_list]
        # The renal dataset has duplicate columns for the same protein in some cases
        # We'll keep only the first column with that protein name
        proteomics_data = proteomics_data.loc[:, ~proteomics_data.columns.duplicated()]
        plot_data = plot_data.join(proteomics_data)

        # Extract all the normal matched samples into a new df and remove the 'N' to match them with tumor samples
        normal = plot_data.copy()
        normal = normal.loc[normal['Sample_Tumor_Normal'] == 'Normal']
        normal['Patient_ID'] = normal['Patient_ID'].str.replace('NC', 'C')
        normal.index = normal['Patient_ID']

        # Create new column to count the number of proteins with increased/decreased/nochange levels in tumor samples
        # Remove the normal samples from plot_data
        plot_data[expression_change] = 0 
        plot_data = plot_data.loc[plot_data['Sample_Tumor_Normal'] == 'Tumor']
        plot_data.index = plot_data['Patient_ID']
        
        # Remove unmatched samples from each dataframe
        plot_data = plot_data.loc[plot_data['Patient_ID'].isin(normal['Patient_ID'])]
        normal = normal.loc[normal['Patient_ID'].isin(plot_data['Patient_ID'])] 
    
    else: 
        print('Error: invalid or unimplemented cancer type')
        return
    
    # Fill out the Num_* column by comparing the protein expression for matched patients in plot_data/normal dfs
    for protein in protein_list:
        mean_diff = np.abs(np.mean(plot_data[protein]) - np.mean(normal[protein]))
        if expression_change == 'increased':
            plot_data.loc[((plot_data[protein] - normal[protein]) >= mean_diff), expression_change] += 1
        elif expression_change == 'decreased':
            plot_data.loc[((plot_data[protein] - normal[protein]) <= (-1*mean_diff)), expression_change] += 1
        else:
            plot_data.loc[((plot_data[protein] - normal[protein]) <= mean_diff)
                      & ((plot_data[protein] - normal[protein]) >= (-1*mean_diff)), expression_change] += 1

    # Call another function to create the plot
    if boxplot:
        create_complex_clinical_boxplot(plot_data[[clinical_feature, expression_change]], clinical_feature, expression_change, complex_name)
    else:
        create_complex_clinical_plot(plot_data, clinical_feature, protein_list, bin_size, expression_change, complex_name, show_proportions)

def create_complex_clinical_boxplot(plot_data, clinical_feature, expression_change, complex_name):
    clinical_vals = list(set(plot_data[clinical_feature].dropna()))
    for i in range(0, len(clinical_vals)):
        for j in range(i+1, len(clinical_vals)):
            val_1 = plot_data.loc[plot_data[clinical_feature] == clinical_vals[i]]
            val_2 = plot_data.loc[plot_data[clinical_feature] == clinical_vals[j]]
            pval = scipy.stats.ttest_ind(val_1[expression_change], val_2[expression_change])[1]
            print('t-test result ' + str(clinical_vals[i]) + ' vs ' + str(clinical_vals[j]) + ': ' + str(pval))
    a4_dims = (10, 10)
    fig, ax = plt.subplots(figsize=a4_dims)
    sns.set(font_scale = 1.5, style = 'white')
    boxplt = sns.boxplot(data = plot_data,
                           x = clinical_feature,
                           y = expression_change,
                           showfliers = False)
    boxplt = sns.stripplot(data=plot_data, 
                           x = clinical_feature, 
                           y = expression_change, 
                           dodge=True, 
                           jitter=True, 
                           color='.3')
    
    boxplt.set_title(complex_name, fontsize = '20')
    boxplt.set_ylabel('Number of proteins ' + expression_change + ' in ' + complex_name)

"""
create_complex_clinical_plot
----------------------------
Function called by plot_complex_clinical; this prepares the dataframe for plotting and creates the actual plot

Parameters are the same as those passed to plot_complex_clinical

"""
def create_complex_clinical_plot(plot_data, clinical_feature, protein_list, bin_size, expression_change, complex_name, show_proportions):
    # Drop the proteomics data (not necessary at this point)
    plot_data.drop(protein_list, axis = 1, inplace=True)
    category_col_title = expression_change + ' category'
    sns.set(font_scale = 1, style='white')
    # If indicated, bin the patients into groups based on number of proteins increased/decreased/nochange
    # This makes it easier to view, otherwise there are lots of ticks on the x-axis
    if bin_size > 1:
        # Create the list of cutoff values for each bin based on the given bin_size
        bin_vals = list(range(0, len(protein_list), bin_size))
        # Create the list of strings to use as category names
        category_names = []
        for i in range(0, len(bin_vals) - 1):
            category_names.append(str(bin_vals[i]) + '-' + str(bin_vals[i+1] - 1))
        category_names.append(str(bin_vals[len(bin_vals)-1]) + '+')
        
        # Place each sample into the correct category
        plot_data[category_col_title] = category_names[0]
        for i in range(1, len(bin_vals) - 1):
            plot_data.loc[((plot_data[expression_change] > bin_vals[i]) & (plot_data[expression_change] <= bin_vals[i+1])), category_col_title] = category_names[i]
            
        # Create the appropriate plot based on the desired differential expression type
        # increased means increased in tumor samples
        a4_dims = (10, 10)
        fig, ax = plt.subplots(figsize=a4_dims)
        x_order = category_names
        if show_proportions:
            percent_df = (plot_data[category_col_title]
                            .groupby(plot_data[clinical_feature])
                            .value_counts(normalize=True)
                            .rename('Proportion ' + clinical_feature)
                            .reset_index())
            countplot = sns.barplot(data=percent_df, x=category_col_title, y=('Proportion ' + clinical_feature), hue=clinical_feature, order = x_order)
            countplot.set_ylabel('Proportion of Individuals', fontsize = '15')
            tick_spacing = 0.05
            
        else:
            countplot = sns.countplot(x=category_col_title, hue=clinical_feature, data=plot_data, order = x_order)
            countplot.set_ylabel('Number of Individuals', fontsize = '15')
            tick_spacing = 2
    
    else:
        # Same idea as above, just without binning the patients
        plot_data[category_col_title] = plot_data[expression_change].astype(str)
        a4_dims = (10, 10)
        fig, ax = plt.subplots(figsize=a4_dims)
        x_axis_labels = list(set(plot_data[category_col_title].astype(int).sort_values()))
        x_axis_labels = [str(i) for i in x_axis_labels] 
        
        if show_proportions:
            percent_df = (plot_data[category_col_title]
                            .groupby(plot_data[clinical_feature])
                            .value_counts(normalize=True)
                            .rename('Proportion ' + clinical_feature)
                            .reset_index())
            countplot = sns.barplot(data=percent_df, x = category_col_title, y=('Proportion ' + clinical_feature), hue=clinical_feature, order = x_axis_labels)
            countplot.set_ylabel('Proportion of Individuals', fontsize = '15')
            tick_spacing = 0.05
        else:
            countplot = sns.countplot(x=category_col_title, hue=clinical_feature, data=plot_data, order = x_axis_labels)
            tick_spacing = 2
            countplot.set_ylabel('Number of Individuals', fontsize = '15')
        
    countplot.set_xlabel('Number of proteins with ' + expression_change + ' expression in tumor samples', fontsize = '15')      
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    countplot.set_title(complex_name, fontsize = '25')

def plot_subtype_complex(raw_protein_list, complex_name = '', cancer_type = 'en', simplify_subtype = True):
    """Plot differences in protein expression for proteasome proteins, for each sample, separated by MSI-H and MSS samples. Only implemented for colon and endometrial."""
    protein_list = []
    # Make sure all the proteins in protein_list are in the proteomics df
    for protein in raw_protein_list:
        if protein in en.get_proteomics().columns:
            protein_list.append(protein)
    if cancer_type == 'en':
        # Start out with a df including the patient IDs and molecular subtype
        source_data = en.get_derived_molecular()[['Genomics_subtype']]
        source_data['Sample_ID'] = source_data.index
        
        # Add protein expression values for all proteins in the given list
        proteomics_data = en.get_proteomics()[protein_list]
        patient_ids = en.get_clinical()[['Patient_ID', 'Proteomics_Tumor_Normal']]
        source_data = source_data.join(patient_ids)
        source_data = source_data.join(proteomics_data)

        # Extract all the normal matched samples into a new df and remove the 'N' to match them with tumor samples
        normal = source_data.copy()
        normal = normal.loc[normal['Proteomics_Tumor_Normal'].str.contains('Adjacent_normal')]
        normal.index = normal['Patient_ID']

        # Remove the normal samples from the tumor df
        tumor = source_data.loc[source_data['Proteomics_Tumor_Normal'] == 'Tumor']
        tumor.index = tumor['Patient_ID']

        # Remove unmatched samples from each dataframe
        tumor = tumor.loc[tumor['Patient_ID'].isin(normal['Patient_ID'])]
        normal = normal.loc[normal['Patient_ID'].isin(tumor['Patient_ID'])] 
        
        # Create the plot_data with ratios tumor/normal for each protein
        plot_data = tumor.copy()[['Genomics_subtype', 'Sample_ID']]
        for protein in protein_list:
            plot_data[protein + ' (T/N)'] = tumor[protein] - normal[protein]
            
        # Create only 3 categories of subtype: MSI-H, POLE, and Other
        if simplify_subtype:
            plot_data.loc[((plot_data['Genomics_subtype'] != 'MSI-H') 
                                       & (plot_data['Genomics_subtype'] != 'POLE')), 
                                      'Genomics_subtype'] = 'Other'
        plot_data = pd.melt(plot_data, id_vars = ['Genomics_subtype', 'Sample_ID'], var_name = 'Protein', value_name = 'Ratio')

        pole = plot_data.loc[plot_data['Genomics_subtype'] == 'POLE'].groupby('Sample_ID').mean()
        msi = plot_data.loc[plot_data['Genomics_subtype'] == 'MSI-H'].groupby('Sample_ID').mean()
        other = plot_data.loc[plot_data['Genomics_subtype'] == 'Other'].groupby('Sample_ID').mean()

        mwu_msi_pole = scipy.stats.mannwhitneyu(msi['Ratio'], pole['Ratio'])[1]
        mwu_msi_other = scipy.stats.mannwhitneyu(msi['Ratio'], other['Ratio'])[1]
        mwu_pole_other = scipy.stats.mannwhitneyu(pole['Ratio'], other['Ratio'])[1]
        print('p-value for Mann-Whitney rank test test MSI-H vs POLE ' + str(mwu_msi_pole))
        print('p-value for Mann-Whitney rank test test MSI-H vs Other ' + str(mwu_msi_other))
        print('p-value for Mann-Whitney rank test test POLE vs Other ' + str(mwu_pole_other))
        
    elif cancer_type == 'colon':
        source_data = colon.get_clinical()[['Mutation_Phenotype', 'Patient_ID']]
        source_data['Sample_ID'] = source_data.index
        
        proteomics_data = colon.get_proteomics()[protein_list]
        source_data = source_data.join(proteomics_data)

        # Extract all the normal matched samples into a new df and remove the 'N' to match them with tumor samples
        normal = source_data.copy()
        normal = normal.loc[normal['Patient_ID'].str.contains('N')]
        normal['Patient_ID'] = normal['Patient_ID'].str.replace('N', '')
        normal.index = normal['Patient_ID']

        # Remove the normal samples from the tumor df
        tumor = source_data.loc[~source_data['Patient_ID'].str.contains('N')]
        tumor.index = tumor['Patient_ID']

        # Remove unmatched samples from each dataframe
        tumor = tumor.loc[tumor['Patient_ID'].isin(normal['Patient_ID'])]
        normal = normal.loc[normal['Patient_ID'].isin(tumor['Patient_ID'])] 
        
        plot_data = tumor.copy()[['Sample_ID']]
        plot_data['Genomics_subtype'] = tumor['Mutation_Phenotype']
        print(tumor)
        print(normal)
        for protein in protein_list:
            plot_data[protein + ' (T/N)'] = (tumor[protein] + 10) - (normal[protein] + 10)
        print(plot_data)
        plot_data = pd.melt(plot_data, id_vars = ['Genomics_subtype', 'Sample_ID'], var_name = 'Protein', value_name = 'Ratio')
        
        mss = plot_data.loc[plot_data['Genomics_subtype'] == 'MSS'].groupby('Sample_ID').mean()
        msi = plot_data.loc[plot_data['Genomics_subtype'] == 'MSI-H'].groupby('Sample_ID').mean()

        whitneyutest = scipy.stats.mannwhitneyu(mss['Ratio'], msi['Ratio'])[1]
        print('p-value for Mann-Whitney rank test test MSS vs MSI-H: ' + str(whitneyutest))
        
    else:
        print('Error: invalid or unimplemented cancer type. Options are: en')
        return
    
    create_complex_boxplot(plot_data, complex_name, protein_list)

def create_complex_boxplot(plot_data, complex_name, protein_list):
    a4_dims = (100, 50)
    fig, ax = plt.subplots(figsize=a4_dims)
    sns.set(font_scale = 5, style = 'white')
    boxplt = sns.boxplot(data = plot_data,
                           x = 'Genomics_subtype',
                           y = 'Ratio',
                           hue = 'Sample_ID',
                           showfliers = False)
    boxplt = sns.stripplot(data=plot_data, 
                           x='Genomics_subtype', 
                           y='Ratio', 
                           hue = 'Sample_ID',
                           dodge=True, 
                           jitter=True, 
                           color='.3')
    
    boxplt.set_title(complex_name, fontsize = '300')
    boxplt.set_xlabel('Genomic Subtype', fontsize = '100')
    boxplt.set_ylabel('Ratio (Tumor/Normal)', fontsize = '100')
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:0], labels[0:0])
    boxplt.set_title(complex_name)

def ttest_all_complexes(alpha = 0.05, cancer_type = None):
    """For the specified cancer type, create a dataframe that tells you which protein complexes had significantly increased, significantly decreased, or no change in expression levels for all protein complexes."""

    results = pd.DataFrame(columns = ['complex_name', 'total_protein_num', 'num_increased_tumor', 
                                      'num_decreased_tumor', 'no_proteomics_data'])
    
    # Prepare dataframes for analysis (make them all look the same, remove duplicates, etc)
    if not cancer_type: 
        print('Error: cancer_type is a required parameter. Options are: colon, en, ov, and renal')
        return
    
    if cancer_type == 'colon':
        proteomics_data = colon.join_metadata_to_omics('clinical', 'proteomics', ['Sample_Tumor_Normal'])
    
    elif cancer_type == 'en':
        proteomics_data = en.join_metadata_to_omics('clinical', 'proteomics', ['Proteomics_Tumor_Normal'])
        proteomics_data.columns = proteomics_data.columns.str.replace('Proteomics_Tumor_Normal', 'Sample_Tumor_Normal')
        proteomics_data.loc[proteomics_data['Sample_Tumor_Normal'] == 'Adjacent_normal', 'Sample_Tumor_Normal'] = 'Normal'
    
    elif cancer_type == 'ov':
        proteomics_data = ov.join_metadata_to_omics('clinical', 'proteomics', ['Sample_Tumor_Normal'])
    
    elif cancer_type == 'renal':
        proteomics_data = renal.join_metadata_to_omics('clinical', 'proteomics', ['Sample_Tumor_Normal'])
    
    else:
        print('Error: cancer_type must be colon, en, ov, or renal')
        return
    
    proteomics_data.columns = proteomics_data.columns.str.replace('_proteomics', '')
    proteomics_data = proteomics_data.loc[:,~proteomics_data.columns.duplicated()]
    
    for complex_name, proteins in subunitNames.items():
        num_sig_increased = 0
        num_sig_decreased = 0
        num_no_data = 0
        tumor = proteomics_data.loc[proteomics_data['Sample_Tumor_Normal'] == 'Tumor']
        normal = proteomics_data.loc[proteomics_data['Sample_Tumor_Normal'] == 'Normal']
        
        for protein in proteins:
            if protein in proteomics_data.columns:
                tumor_vals = tumor[protein].dropna()
                normal_vals = normal[protein].dropna()
                if len(tumor_vals) > 2 and len(normal_vals) > 2:
                    p_val = scipy.stats.ttest_ind(tumor_vals, normal_vals)[1]
                    if p_val < alpha:
                        if (np.mean(tumor[protein]) > np.mean(normal[protein])):
                            num_sig_increased += 1
                        elif (np.mean(tumor[protein]) < np.mean(normal[protein])):
                            num_sig_decreased += 1
            else:
                num_no_data += 1
        results = results.append({'complex_name': complex_name,
                      'total_protein_num': len(proteins),
                      'num_increased_tumor': num_sig_increased,
                      'num_decreased_tumor': num_sig_decreased,
                      'no_proteomics_data': num_no_data}, ignore_index = True)
        
    results['percent_increased'] = (results['num_increased_tumor'] / results['total_protein_num'])*100
    results['percent_decreased'] = (results['num_decreased_tumor'] / results['total_protein_num'])*100
    return results

def make_histogram(results, show_change = 'percent_increased', size_cutoffs = [10, 20], norm = False, show_complex = None):
    """Create a stacked/overlaid histogram that displays the output of ttest_all_complexes."""
    a4_dims = (10, 10)
    fig, ax = plt.subplots(figsize=a4_dims)
    color_vals = ['skyblue', 'olive', 'gold', 'teal', 'red']
    if norm:
        if size_cutoffs is not None:
            # Get the column that shows the desired type of change, for the complexes with the minimum number of proteins
            plot_data = list(results.loc[results['total_protein_num'] < size_cutoffs[0]][show_change]) 

            # The following lines of code generate a different histogram for the complexes falling withing each range between the size_cutoffs
            # x axis is bins grouping by increments of the change type we're showing (e.g. percent increase in expression)
            # y axis is the proportion of proteins with that range of our chosen change type. But I'm not quite sure, because the proportions seem weird.
            distplot = sns.distplot(plot_data, color = color_vals[0], label = 'n < ' + str(size_cutoffs[0]), kde = False, norm_hist = True, bins = 20)
            for i in range(1, len(size_cutoffs)):
                plot_data = results.loc[(results['total_protein_num'] >= size_cutoffs[i-1]) & (results['total_protein_num'] < size_cutoffs[i])]
                plot_data = list(plot_data[show_change])
                distplot = sns.distplot(plot_data, color = color_vals[i], kde = False, norm_hist = True, label = (str(size_cutoffs[i-1]) + ' <= n < ' + str(size_cutoffs[i])), bins = 20)
            plot_data = list(results.loc[results['total_protein_num'] >= size_cutoffs[len(size_cutoffs) - 1]][show_change])
            distplot = sns.distplot(plot_data, color = color_vals[len(color_vals) - 1], kde = False, norm_hist = True, label = 'n >= ' + str(size_cutoffs[len(size_cutoffs)-1]), bins = 20)
            distplot.set(xlabel=show_change + "expression of complex's proteins", ylabel='Frequency') # I added this
            plt.legend()

        else:
            plot_data = list(results[show_change])
            distplot = sns.distplot(plot_data, color = color_vals[0], kde = False, norm_hist = True, label = 'No size cutoffs')
            plt.legend()
    else:
        if size_cutoffs is not None:
            plot_data = list(results.loc[results['total_protein_num'] < size_cutoffs[0]][show_change])
            distplot = sns.distplot(plot_data, color = color_vals[0], norm_hist = False, kde = False, label = 'n < ' + str(size_cutoffs[0]), bins = 20)
            for i in range(1, len(size_cutoffs)):
                plot_data = results.loc[(results['total_protein_num'] >= size_cutoffs[i-1]) & (results['total_protein_num'] < size_cutoffs[i])]
                plot_data = list(plot_data[show_change])
                distplot = sns.distplot(plot_data, color = color_vals[i], norm_hist = False, kde = False, label = (str(size_cutoffs[i-1]) + ' <= n < ' + str(size_cutoffs[i])), bins = 20)
            plot_data = list(results.loc[results['total_protein_num'] >= size_cutoffs[len(size_cutoffs) - 1]][show_change])
            distplot = sns.distplot(plot_data, color = color_vals[len(color_vals) - 1], norm_hist = False, kde = False, label = 'n >= ' + str(size_cutoffs[len(size_cutoffs)-1]))
            plt.legend()
        else:
            plot_data = list(results[show_change])
            distplot = sns.distplot(plot_data, color = color_vals[0], norm_hist = False, kde = False, label = 'No size cutoffs', bins = 20)
            plt.legend()
    
    if show_complex is not None:
        num_changed = int(results.loc[results['complex_name'] == show_complex][show_change])
        plt.axvline(num_changed, 0,0.75, color = 'black')
    
    if show_change == 'percent_increased':
        distplot.set_title('Percent of proteins in complex significantly increased in tumor')
    else:
        distplot.set_title('Percent of proteins in complex significantly decreased in tumor')
        
    plt.legend().set_title('Total number of proteins in complex')

'''
pancancer_boxplot
-----------------
Create boxplots of tumor vs normal for all proteins in a given complex, or for a list of proteins, for all cancer types

'''
def pancancer_boxplot(complex_name = None, protein_list = None, show_ttest = True):
    
    if complex_name:
        protein_list = subunitNames[complex_name]
    else:
        complex_name = 'Unique Protein List'
        
    proteomics_dict = get_pancancer_proteomics(protein_list)
    if show_ttest:
        print("t test results for tumor vs. normal:")
        pancancer_ttest(complex_name, protein_list)
        
    for key, val in proteomics_dict.items():
        a4_dims = (11.7, 8.27)
        fig, ax = plt.subplots(figsize = a4_dims)
        plot_data = pd.melt(val, id_vars = ['Sample_Tumor_Normal'], var_name = 'Protein', value_name = 'Proteomics')
        boxplt = sns.boxplot(data = plot_data, x = 'Protein', y = 'Proteomics', hue = 'Sample_Tumor_Normal', showfliers = False)
        if len(protein_list) < 10: 
            boxplt = sns.stripplot(data = plot_data, x = 'Protein', y = 'Proteomics', hue = 'Sample_Tumor_Normal', dodge = True, jitter = True, color = '.3')
        else:
            boxplt.set(xticklabels = [])
        boxplt.set_title(key, fontsize = '20')
        boxplt.set_ylabel('Proteomics', fontsize = '15')
        boxplt.set_xlabel(complex_name, fontsize = '15')
        handles, labels = ax.get_legend_handles_labels()
        l = plt.legend(handles[0:2], labels[0:2])
        plt.show()

def pancancer_ttest(complex_name, protein_list):
    """For a particular protein, test for a significant difference in expression between tumor and normal cells, across all cancers."""
    proteomics_dict = get_pancancer_proteomics(protein_list)
    print(complex_name + '\n')
    for key, val in proteomics_dict.items():
        print(key + ' cancer:')
        tumor = val.loc[val['Sample_Tumor_Normal'] == 'Tumor']
        normal = val.loc[val['Sample_Tumor_Normal'] == 'Normal']
        for protein in protein_list:
            if len(tumor[protein].dropna()) > 2 and len(normal[protein].dropna()) > 2:
                pval = scipy.stats.ttest_ind(tumor[protein].dropna(), normal[protein].dropna())[1]
                print(protein + ': ' + str(pval))
        print('\n')

def get_pancancer_proteomics(protein_list):
    """For each cancer type, get the proteomics dataframe, joined with the Sample_Tumor_Normal column.
    
    Returns:
    dict: Keys are cancer types, values are proteomics dataframes joined with Sample_Tumor_Normal column.
    """
    proteomics_dict = {}
    
    with contextlib.redirect_stdout(None):
        colon_proteomics = colon.join_metadata_to_omics('clinical', 'proteomics', ['Sample_Tumor_Normal'], protein_list)
        colon_proteomics.columns = colon_proteomics.columns.str.replace('_proteomics', '')
        colon_proteomics = colon_proteomics.loc[:,~colon_proteomics.columns.duplicated()]
        proteomics_dict['Colon'] = colon_proteomics

        en_proteomics = en.join_metadata_to_omics('clinical', 'proteomics', ['Proteomics_Tumor_Normal'], protein_list)
        en_proteomics.columns = en_proteomics.columns.str.replace('Proteomics_Tumor_Normal', 'Sample_Tumor_Normal')
        en_proteomics.columns = en_proteomics.columns.str.replace('_proteomics', '')
        en_proteomics = en_proteomics.loc[:,~en_proteomics.columns.duplicated()]
        en_proteomics.loc[en_proteomics['Sample_Tumor_Normal'] == 'Adjacent_normal', 'Sample_Tumor_Normal'] = 'Normal'
        en_proteomics = en_proteomics.loc[(en_proteomics['Sample_Tumor_Normal'] == 'Tumor') | (en_proteomics['Sample_Tumor_Normal'] == 'Normal')]
        proteomics_dict['Endometrial'] = en_proteomics

        renal_proteomics = renal.join_metadata_to_omics('clinical', 'proteomics', ['Sample_Tumor_Normal'], protein_list)
        renal_proteomics.columns = renal_proteomics.columns.str.replace('_proteomics', '')
        renal_proteomics = renal_proteomics.loc[:,~renal_proteomics.columns.duplicated()]
        proteomics_dict['Renal'] = renal_proteomics

        ov_proteomics = ov.join_metadata_to_omics('clinical', 'proteomics', ['Sample_Tumor_Normal'], protein_list)
        ov_proteomics.columns = ov_proteomics.columns.str.replace('_proteomics', '')
        ov_proteomics = ov_proteomics.loc[:,~ov_proteomics.columns.duplicated()]
        proteomics_dict['Ovarian'] = ov_proteomics
    
    return proteomics_dict

def pancancer_analyze_complex(complex_name, include_borderline = False):
    """Determine which proteins in a complex have significantly different expression between tumor and normal samples, for all cancer types."""

    protein_list = subunitNames[complex_name]
    alpha = 0.05 / len(protein_list)
    if include_borderline:
        print('Borderline changes included; alpha set to 0.05\n')
        alpha = 0.05
    
    increased_all = []
    decreased_all = []
    no_change_all = []
    mixed_results = []
    missing_data = []
    
    proteomics_dict = get_pancancer_proteomics(protein_list)
    
    for protein in protein_list:
        num_increased = 0
        num_decreased = 0
        num_no_change = 0
        no_data = False
        for key, val in proteomics_dict.items():
            tumor = val.loc[val['Sample_Tumor_Normal'] == 'Tumor'][protein].dropna()
            normal = val.loc[val['Sample_Tumor_Normal'] == 'Normal'][protein].dropna()
            if len(tumor) > 2 and len(normal) > 2:
                pval = scipy.stats.ttest_ind(tumor, normal)[1]
                if pval <= alpha:
                    if np.mean(tumor) > np.mean(normal):
                        num_increased += 1
                    else:
                        num_decreased += 1
                else:
                    num_no_change += 1
            else:
                no_data = True
                missing_data.append(protein)
                break
                
        if not no_data:     
            if num_decreased + num_no_change == 0:
                increased_all.append(protein)

            elif num_increased + num_no_change == 0:
                decreased_all.append(protein)

            elif num_decreased + num_increased == 0:
                no_change_all.append(protein)

            else:
                mixed_results.append(protein)
            
    print('Increased in all cancers: ')
    print(increased_all)
    print('Decreased in all cancers: ')
    print(decreased_all)
    print('No change in all cancers: ')
    print(no_change_all)
    print('Mixed results: ')
    print(mixed_results)
    print('Missing proteomics data in at least one cancer: ')
    print(missing_data)
    
    return {'increased': increased_all, 'decreased': decreased_all, 'no change': no_change_all, 'mixed': mixed_results, 'missing': missing_data}

def pancancer_plot_proteins(complex_name = None, protein_list = None, show_ttest = False):
    """For a particular complex or a list of proteins, plot the difference between tumor and normal expression across all cancers, with a separate plot for each protein."""
    if complex_name:
        protein_list = subunitNames[complex_name]
    else:
        complex_name = 'Unique Protein List'
        
    if show_ttest:
        pancancer_ttest(complex_name, protein_list)
        
    proteomics_dict = get_pancancer_proteomics(protein_list)
    
    for protein in protein_list:
        plot_data = pd.DataFrame(columns = [protein, 'Sample_Tumor_Normal', 'Cancer_Type'])
        for key, val in proteomics_dict.items():
            cancer_data = val[[protein, 'Sample_Tumor_Normal']].copy()
            cancer_data['Cancer_Type'] = key
            plot_data = plot_data.append(cancer_data)
            
        a4_dims = (11.7, 8.27)
        fig, ax = plt.subplots(figsize = a4_dims)
        boxplt = sns.boxplot(data = plot_data, x = 'Cancer_Type', y = protein, hue = 'Sample_Tumor_Normal', showfliers = False)
        boxplt = sns.stripplot(data = plot_data, x = 'Cancer_Type', y = protein, hue = 'Sample_Tumor_Normal', dodge = True, jitter = True, color = '.3')
        boxplt.set_title(protein, fontsize = '20')
        boxplt.set_ylabel('Proteomics', fontsize = '15')
        boxplt.set_xlabel('Cancer Type', fontsize = '15')
        handles, labels = ax.get_legend_handles_labels()
        l = plt.legend(handles[0:2], labels[0:2])
        plt.show()

def make_separated_boxplot(protlist1, protlist2, cancer_type, show_data = False):
    plot_data = get_pancancer_proteomics(protlist1 + protlist2)[cancer_type]
    plot_data[''] = np.nan
    plot_data = pd.melt(plot_data, id_vars = ['Sample_Tumor_Normal'], var_name = 'Protein', value_name = 'Proteomics')
    plot_order = protlist1 + [''] + protlist2
    
    a4_dims = (11.7, 8.27)
    fig, ax = plt.subplots(figsize = a4_dims)
    boxplt = sns.boxplot(data = plot_data, x = 'Protein', y = 'Proteomics', hue = 'Sample_Tumor_Normal', order = plot_order, showfliers = False)
    if show_data:
        boxplt = sns.stripplot(data = plot_data, x = 'Protein', y = 'Proteomics', hue = 'Sample_Tumor_Normal', dodge = True, jitter = True, color = '0.3', order = plot_order)
        handles, labels = ax.get_legend_handles_labels()
        l = plt.legend(handles[0:2], labels[0:2])
    boxplt.set_title(cancer_type + ' Cancer', fontsize = '20')
    boxplt.set_ylabel('Proteomics', fontsize = '15')
    boxplt.set_xlabel('')

def plot_phosphosites(site_list, cancer_type, show_ttest = True):
    if cancer_type == 'renal':
        clin = renal.get_clinical()[['Sample_Tumor_Normal']]
        phosdf = renal.get_phosphoproteomics()[site_list]
        plotdf = phosdf.join(clin)
        plot_data = pd.melt(plotdf, id_vars = ['Sample_Tumor_Normal'], var_name = 'Site', value_name = 'Phosphorylation')
    else:
        print('Error: unimplemented cancer type')
        return
    
    if show_ttest:
        for site in site_list:
            tumor = plotdf.loc[plotdf['Sample_Tumor_Normal'] == 'Tumor'][site].dropna()
            normal = plotdf.loc[plotdf['Sample_Tumor_Normal'] == 'Normal'][site].dropna()
            pval = scipy.stats.ttest_ind(tumor, normal)[1]
            print('T-test tumor vs normal site ' + site + ': ' + str(pval))
    
    a4_dims = (11.7, 8.27)
    fig, ax = plt.subplots(figsize = a4_dims)
    boxplt = sns.boxplot(data = plot_data, x = 'Site', y = 'Phosphorylation', hue = 'Sample_Tumor_Normal', showfliers = False)
    boxplt = sns.stripplot(data = plot_data, x = 'Site', y = 'Phosphorylation', hue = 'Sample_Tumor_Normal', dodge = True, jitter = True, color = '0.3')
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:2], labels[0:2])
    boxplt.set_title(cancer_type + ' cancer', fontsize = '20')
    boxplt.set_ylabel('Phosphorylation', fontsize = '15')
    boxplt.set_xlabel('Phosphosite')

# Thanks StackOverflow! 
# https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows

def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    #df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df
