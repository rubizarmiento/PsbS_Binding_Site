#!/usr/bin/env python3
"""
Apply SQL migrations to database CSV files.
Reads .sql files and applies transformations using SQLite.
"""

import argparse
import sqlite3
import pandas as pd
from pathlib import Path
import re


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Apply SQL migrations to database CSV files."
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
        help="Output CSV file path (modified database). If not provided, modifies in place.",
    )
    parser.add_argument(
        "--migrations",
        type=str,
        nargs='+',
        required=True,
        help="One or more SQL migration files to apply.",
    )
    parser.add_argument(
        "--verbose",
        action='store_true',
        help="Print SQL statements as they are executed.",
    )
    return parser.parse_args()


def parse_sql_file(sql_file):
    """
    Parse SQL file and return individual statements.
    Handles comments and multi-line statements.
    
    Parameters:
    -----------
    sql_file : str
        Path to SQL file
    
    Returns:
    --------
    list of str
        Individual SQL statements
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


def apply_migrations_to_csv(input_csv, output_csv, migration_files, verbose=False):
    """
    Apply SQL migrations from files to a CSV.
    
    Parameters:
    -----------
    input_csv : str
        Path to input CSV file
    output_csv : str
        Path to output CSV file
    migration_files : list of str
        Paths to SQL migration files
    verbose : bool
        Whether to print SQL statements
    """
    # Read CSV into pandas
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"  Loaded {len(df)} rows, {len(df.columns)} columns")
    
    # Create in-memory SQLite database
    conn = sqlite3.connect(':memory:')
    
    # Load dataframe into SQLite
    df.to_sql('data', conn, index=False, if_exists='replace')
    print(f"  Loaded into SQLite table 'data'")
    
    # Apply each migration file
    for migration_file in sorted(migration_files):
        print(f"\nApplying migration: {Path(migration_file).name}")
        
        # Parse SQL statements from file
        statements = parse_sql_file(migration_file)
        
        # Execute each statement
        for i, statement in enumerate(statements, 1):
            if verbose:
                print(f"\n  Statement {i}:")
                print(f"  {statement[:100]}..." if len(statement) > 100 else f"  {statement}")
            
            try:
                cursor = conn.execute(statement)
                
                # For ALTER TABLE or UPDATE, show rows affected
                if statement.upper().startswith('UPDATE'):
                    rows_affected = cursor.rowcount
                    print(f"    ✓ Updated {rows_affected} rows")
                elif statement.upper().startswith('ALTER'):
                    print(f"    ✓ Column added")
                else:
                    print(f"    ✓ Executed")
                    
                conn.commit()
                
            except sqlite3.OperationalError as e:
                error_msg = str(e).lower()
                if 'duplicate column' in error_msg:
                    print(f"    ⚠ Warning: Column already exists, skipping...")
                elif 'no such column' in error_msg:
                    print(f"    ✗ Error: {e}")
                    print(f"    Tip: Make sure the column exists in your data")
                    raise
                else:
                    print(f"    ✗ Error: {e}")
                    raise
            except Exception as e:
                print(f"    ✗ Error: {e}")
                raise
    
    # Read back the modified data
    print(f"\nReading modified data from SQLite...")
    result_df = pd.read_sql_query("SELECT * FROM data", conn)
    print(f"  Result: {len(result_df)} rows, {len(result_df.columns)} columns")
    
    # Close connection
    conn.close()
    
    # Write output
    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_csv, index=False)
    print(f"\n✓ Saved to {output_csv}")
    
    # Show new columns
    new_columns = set(result_df.columns) - set(df.columns)
    if new_columns:
        print(f"  New columns added: {', '.join(sorted(new_columns))}")
    
    return result_df


def main():
    args = parse_arguments()
    
    # Determine output path
    output_csv = args.output_csv if args.output_csv else args.input_csv
    
    if output_csv == args.input_csv:
        print("⚠ Warning: Modifying file in place!")
    
    # Apply migrations
    apply_migrations_to_csv(
        input_csv=args.input_csv,
        output_csv=output_csv,
        migration_files=args.migrations,
        verbose=args.verbose
    )
    
    print("\n✓ Done!")


if __name__ == '__main__':
    main()
