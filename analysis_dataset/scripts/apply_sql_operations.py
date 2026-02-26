#!/usr/bin/env python3
"""
Author
------
Rubi Zarmiento-Garcia

Applies SQL Operations to CSV files

Example
-----
Command line::

    python apply_sql_migrations.py \\
        --input_csv /path/to/input.csv \\
        --output_csv /path/to/output.csv \\
        --operations /path/to/operation1.sql /path/to/operation2.sql

Required Arguments
------------------
--input_csv : str
    Path to input CSV file to transform
--output_csv : str
    Path to output CSV file for results
--operations : list of str
    One or more SQL operation files to apply (applied in sorted order)
"""

import argparse
import sqlite3
import pandas as pd
from pathlib import Path
import re


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Apply SQL operations to database CSV files."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Input CSV file path (database to modify).",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        required=True,
        help="Output CSV file path (modified database).",
    )
    parser.add_argument(
        "--operations",
        type=str,
        nargs='+',
        required=True,
        help="One or more SQL operation files to apply.",
    )
    return parser.parse_args()


def parse_sql_file(sql_file):
    """
    Parse SQL file and extract individual statements.
    
    Handles SQL comments (both single-line -- and multi-line /* */) and
    splits statements by semicolon delimiters. Empty statements are filtered out.

    Parameters
    ----------
    sql_file : str
        Path to SQL file to parse

    Returns
    -------
    list of str
        Individual SQL statements with comments removed and whitespace stripped

    Notes
    -----
    - Single-line comments (--) are removed
    - Multi-line comments (/* */) are removed
    - Statements are split by semicolon (;)
    - Empty or whitespace-only statements are filtered out
    """
    with open(sql_file, 'r') as f:
        content = f.read()
    
    # Remove single-line comments (-- style)
    content = re.sub(r'--.*$', '', content, flags=re.MULTILINE)
    
    # Remove multi-line comments (/* */ style)
    content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
    
    # Split by semicolon to get individual statements
    statements = [stmt.strip() for stmt in content.split(';')]
    
    # Filter out empty statements
    statements = [stmt for stmt in statements if stmt]
    
    return statements


def apply_operations_to_csv(input_csv, output_csv, operation_files):
    """
    Apply SQL operations from files to a CSV database.
    
    Loads a CSV file into an in-memory SQLite database, applies SQL operations
    from one or more SQL files, and writes the transformed data back to CSV.

    Parameters
    ----------
    input_csv : str
        Path to input CSV file to transform
    output_csv : str
        Path to output CSV file for results
    operation_files : list of str
        Paths to SQL operation files to apply (processed in sorted order)

    Returns
    -------
    pd.DataFrame
        Transformed DataFrame after all operations applied

    Raises
    ------
    Exception
        If SQL statement execution fails or file operations fail

    Notes
    -----
    - CSV is loaded into SQLite table named 'data'
    - Operation files are sorted alphabetically before application
    - Duplicate column errors are silently ignored (idempotent operations)
    - Output directories are created automatically if they don't exist
    - Connection is properly closed after operations complete
    
    Examples
    --------
    >>> result = apply_operations_to_csv(
    ...     'input.csv',
    ...     'output.csv',
    ...     ['ops/01_add_cols.sql', 'ops/02_update.sql']
    ... )
    """
    # Read CSV into pandas
    df = pd.read_csv(input_csv, keep_default_na=False, na_values=[''])
    
    # Create in-memory SQLite database
    conn = sqlite3.connect(':memory:')
    
    # Load dataframe into SQLite
    df.to_sql('data', conn, index=False, if_exists='replace')
    
    # Apply each operation file
    for operation_file in sorted(operation_files):
        # Parse SQL statements from file
        statements = parse_sql_file(operation_file)
        
        # Execute each statement
        for i, statement in enumerate(statements, 1):
            try:
                # Execute the statement
                conn.execute(statement)
                conn.commit()
                
            except Exception as e:
                print(f"ERROR: {e}")
                raise
    
    # Read back the modified data
    result_df = pd.read_sql_query("SELECT * FROM data", conn)
    
    # Close connection
    conn.close()
    
    # Write output
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_csv, index=False)
    print(f"Output file: {output_csv}")
    
    return result_df


def main():
    """
    Main function to apply SQL operations to CSV files.

    Workflow
    --------
    1. Parse command-line arguments
    2. Load input CSV to in-memory SQLite database
    3. Apply SQL operations from files in sorted order
    4. Write transformed data to output CSV
    5. Print output file location

    Notes
    -----
    Script exits with non-zero status if any operation fails.
    """
    args = parse_arguments()
    
    # Apply operations
    apply_operations_to_csv(
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        operation_files=args.operations
    )


if __name__ == '__main__':
    main()
