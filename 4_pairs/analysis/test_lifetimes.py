import pandas as pd
import numpy as np

def extract_lifetimes(df):
    """
    Extract lifetimes from a DataFrame where each column contains binary time series data.
    
    Parameters:
    df (pd.DataFrame): DataFrame with shape (labels x time) containing binary data (0s and 1s)
    
    Returns:
    list: List of arrays, where each array contains the lifetimes (consecutive 1s) for each column
    
    Example:
    If a column has data [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1], 
    it returns [3, 2] representing two lifetimes of length 3 and 2.
    """
    lifetimes_per_column = []
    
    for column in df.columns:
        column_data = df[column].values
        lifetimes = []
        current_lifetime = 0
        
        for value in column_data:
            if value == 1:
                current_lifetime += 1
            else:
                if current_lifetime > 0:
                    lifetimes.append(current_lifetime)
                    current_lifetime = 0
        
        # Don't forget the last lifetime if the sequence ends with 1s
        if current_lifetime > 0:
            lifetimes.append(current_lifetime)
        
        lifetimes_per_column.append(np.array(lifetimes))
    
    return lifetimes_per_column

# Alternative more compact version using groupby
def extract_lifetimes_vectorized(df):
    """
    Vectorized version using pandas groupby for better performance.
    """
    lifetimes_per_column = []
    
    for column in df.columns:
        series = df[column]
        # Create groups where consecutive identical values are grouped together
        groups = (series != series.shift()).cumsum()
        
        # Filter for groups where the value is 1 and get their sizes
        lifetimes = series.groupby(groups).size()[series.groupby(groups).first() == 1].values
        lifetimes_per_column.append(lifetimes)
    
    return lifetimes_per_column

# Example usage:
if __name__ == "__main__":
    # Create sample data
    data = {
        'col1': [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1],
        'col2': [1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        'col3': [0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1]
    }
    
    df = pd.DataFrame(data)
    print("Sample DataFrame:")
    print(df)
    print()
    
    # Extract lifetimes
    lifetimes = extract_lifetimes(df)
    
    print("Lifetimes per column:")
    for i, col_lifetimes in enumerate(lifetimes):
        print(f"Column {df.columns[i]}: {col_lifetimes}")
    
    print("\nUsing vectorized version:")
    lifetimes_vec = extract_lifetimes_vectorized(df)
    for i, col_lifetimes in enumerate(lifetimes_vec):
        print(f"Column {df.columns[i]}: {col_lifetimes}")