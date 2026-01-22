import pandas as pd
import numpy as np

# Create test data: 424 rows (212 per chain), with resid and some lifetime columns
data = {
    'resid': list(range(1, 213)) + list(range(1, 213)),  # Chain A: 1-212, Chain B: 1-212
    'lifetime_1': np.random.rand(424) * 100,
    'lifetime_2': np.random.rand(424) * 100,
    'lifetime_3': np.random.rand(424) * 100
}

df = pd.DataFrame(data)
print("Original DataFrame shape:", df.shape)
print("\nFirst 5 rows (Chain A):")
print(df.head())
print("\nRows 212-217 (Chain B start):")
print(df.iloc[211:217])

# Apply transformation
half = len(df) // 2
df1 = df.loc[0:half-1,:].copy()
df2 = df.loc[half:len(df)-1,:].copy()

print("\ndf1 shape:", df1.shape)
print("df2 shape:", df2.shape)

# Sum corresponding rows from both halves
sum_ns = (df1.iloc[:,1:].values + df2.iloc[:,1:].values).sum(axis=1)
print("\nsum_ns shape:", sum_ns.shape)
print("First 5 sum_ns values:", sum_ns[:5])

# Create result DataFrame
result = df1.iloc[:,0].copy().to_frame()
result['sum_ns'] = sum_ns

print("\nResult shape:", result.shape)
print("\nFirst 5 rows of result:")
print(result.head())

# Duplicate for symmetry
df_concat = pd.concat([result, result], ignore_index=True)
print("\nFinal concatenated shape:", df_concat.shape)
print("\nFirst 5 rows:")
print(df_concat.head())
print("\nRows 212-217:")
print(df_concat.iloc[211:217])
