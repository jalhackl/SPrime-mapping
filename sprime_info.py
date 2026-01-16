import pandas as pd

def compute_sprime_info(df, phased_evalution=True, add_single_columns_at_end = True):
    """
    Process a DataFrame with interval data and compute:
    - individual_haplotype column
    - interval_length
    - total length per individual_haplotype
    - total length per individual_haplotype per segment
    - segments per individual
    - number of segments per individual
    
    Parameters:
        df (pd.DataFrame): Input dataframe with columns 'individual', 'haplotype', 'start', 'end', 'segment'
        phased_evaluatoin: boolean, whether 'individual_haplotype' column is generated with individual and haplotype information.
        Because SPrime is used with phase data, this will be usually set to True. Default: True
        add_single_columns_at_end: boolean, whether at the end the single 'individual' and 'haplotype' columns are added again. Usually not needed, but could be required for analysis. Default: True
    Returns:
        dict: Contains processed DataFrames:
            - df: original dataframe with 'individual_haplotype' and 'interval_length'
            - total_length_ind_hap: total length per individual_haplotype
            - total_length_ind_hap_seg: total length per individual_haplotype per segment
            - segments_per_ind: segments per individual with count

    """
    
    df = df.copy()
    
    # Create individual_haplotype column
    if 'individual_haplotype' not in df.columns:
        if 'haplotype' in df.columns:
            df['individual_haplotype'] = df['individual'].astype(str) + '_' + df['haplotype'].astype(str)
        else:
            df['individual_haplotype'] = df['individual']
            df['haplotype'] = 0

    if phased_evalution:
        ind_column = "individual_haplotype" 
    else:
        ind_column = "individual"
        
    # 2. Create interval_length column
    df['interval_length'] = df['end'] - df['start']
    
    # Total length per individual/haplotype
    total_length_ind_hap = df.groupby(ind_column)['interval_length'] \
                             .sum() \
                             .reset_index(name='total_length')
    
    # Total length per individual/haplotype per segment
    total_length_ind_hap_seg = df.groupby([ind_column, 'segment'])['interval_length'] \
                                 .sum() \
                                 .reset_index(name='total_length')
    
    # Segments per individual/haplotype
    segments_per_ind = df.groupby(ind_column)['segment'] \
                         .apply(lambda x: sorted(x.unique())) \
                         .reset_index(name='segments')
    
    # Number of segments per individual/haplotype
    segments_per_ind['nr_segments'] = segments_per_ind['segments'].apply(len)

    if add_single_columns_at_end:
        ind_map = (
            df[['individual_haplotype', 'individual', 'haplotype']]
            .drop_duplicates()
        )

        total_length_ind_hap = total_length_ind_hap.merge(
            ind_map, on='individual_haplotype', how='left'
        )

        total_length_ind_hap_seg = total_length_ind_hap_seg.merge(
            ind_map, on='individual_haplotype', how='left'
        )

        segments_per_ind = segments_per_ind.merge(
            ind_map, on='individual_haplotype', how='left'
        )
        
            
    return {
        'df': df,
        'total_length_ind_hap': total_length_ind_hap,
        'total_length_ind_hap_seg': total_length_ind_hap_seg,
        'segments_per_ind': segments_per_ind
    }
