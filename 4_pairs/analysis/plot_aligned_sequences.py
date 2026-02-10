import pandas as pd
import random
import matplotlib.pyplot as plt
import numpy as np
from cmap import Colormap

"""
Script to plot aligned sequences with lifetimes heatmap.

Arguments:
--fasta: Path to the aligned FASTA file.
--files: List of CSV files with lifetimes data.
--headers: List of headers corresponding to the CSV files.
--output: Output path for the heatmap image.
--vmin: Minimum value for color normalization.
--vmax: Maximum value for color normalization.
--cmap: Colormap for the heatmap.
--split_every: Number of columns after which to split the table.
--max_resids: Maximum number of residues to display.

Usage example:
python plot_aligned_sequences.py --fasta aligned_sequences.fasta --files lifetimes1.csv lifetimes2.csv --headers header1 header2 --output heatmap.png --vmin 0 --vmax 10000 --cmap Reds --split_every 110


"""

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Plot aligned sequences with lifetimes heatmap.")
    parser.add_argument("--fasta", type=str, required=True, help="Path to the aligned FASTA file.")
    parser.add_argument("--files", type=str, nargs='+', required=True, help="List of CSV files with lifetimes data.")
    parser.add_argument("--headers", type=str, nargs='+', required=True, help="List of headers corresponding to the CSV files.")
    parser.add_argument("--output", type=str, required=True, help="Output path for the heatmap image.")
    parser.add_argument("--vmin", type=float, default=0, help="Minimum value for color normalization.")
    parser.add_argument("--vmax", type=float, default=10000, help="Maximum value for color normalization.")
    parser.add_argument("--cmap", type=str, default='Reds', help="Colormap for the heatmap.")
    parser.add_argument("--split_every", type=int, default=110, help="Number of columns after which to split the table.")
    parser.add_argument("--max_resids", type=int, default=None, help="Maximum number of residues to display.")
    parser.add_argument("--log_transform", action='store_true', help="Apply log10 transformation to lifetime values.")
    parser.add_argument("--edge_color", type=str, default='black', help="Color of table cell edges (default: black).")
    parser.add_argument("--edge_linewidth", type=float, default=1.0, help="Line width of table cell edges (default: 1.0).")
    parser.add_argument("--row_height", type=float, default=2.0, help="Row height scaling factor (default: 2.0).")
    parser.add_argument("--col_width", type=float, default=0.5, help="Column width scaling factor (default: 0.5).")
    parser.add_argument("--show_every_resid", type=int, default=10, help="Show residue ID every N residues (default: 10).")
    parser.add_argument("--font_size", type=float, default=6, help="Font size for table cells (default: 8).")
    parser.add_argument("--row_label_width", type=float, default=0.15, help="Width of row label column (default: 0.15).")
    parser.add_argument("--barplot_output", type=str, default=None, help="Output path for barplot showing vmin/vmax values.")
    args = parser.parse_args()
    return args


def csvs_to_dfs(files_arr, headers=None, log_transform=False):
    dfs = [pd.read_csv(x) for x in files_arr]
    
    # Apply log transformation if requested
    if log_transform:
        for df in dfs:
            if 'sum_ns' in df.columns:
                df['sum_ns'] = np.log10(df['sum_ns'] + 1)
    
    dict_lifetimes = {
    "headers": headers,
    "lifetimes": dfs
    }
    return dict_lifetimes

def fasta_to_dict(file):
    headers = []
    sequences = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                #Remove ">"
                
                headers.append(line[1:])
                sequences.append('')
            elif line:
                if sequences:
                    sequences[-1] += line
    
    #Sequences string to arr
    sequences = [list(x) for x in sequences]

    dict_fasta = {
    "headers": headers,
    "sequences": sequences
    }
    return dict_fasta, headers, sequences

def sequences_to_resids_with_spaces(dict_fasta):
    headers = dict_fasta["headers"]
    sequences = dict_fasta["sequences"]
    arr_type = [[0 if char == "." else 1 for char in sublist] for sublist in sequences]
    arr_resids_and_spaces_fasta = []
    for sublist in arr_type:
        counter = 0
        arr_sublist = []
        for type_res in sublist:
            if type_res != 0:
                counter += 1
                arr_sublist.append(counter)
            else:
                arr_sublist.append(0)
        arr_resids_and_spaces_fasta.append(arr_sublist)
    arr_resids_and_spaces_fasta = np.array(arr_resids_and_spaces_fasta)
    dict_resids_and_spaces_fasta = {
        "headers": headers,
        "sequences": arr_resids_and_spaces_fasta
    }
    
    return dict_resids_and_spaces_fasta

def create_heatmap_table(row_labels, data_2d, vmin, vmax, cmap='viridis', loc='center', figsize=(8, 0.8), split_every=None, color_data=None, label_data=None, max_resids=None, edge_color='black', edge_linewidth=1.0, row_height=0.1, col_width=0.5, font_size=8, row_label_width=0.15):
    """
    Create a colored table where each cell is colored based on its value.
    
    Parameters:
    -----------
    row_labels : list
        Labels for each row
    data_2d : list of lists
        2D array of values (rows x columns) - used as default for colors and labels
    vmin : float
        Minimum value for color normalization
    vmax : float
        Maximum value for color normalization
    cmap : str
        Matplotlib colormap name (default: 'viridis')
    loc : str
        Table location (default: 'center')
    figsize : tuple
        Figure size (width, height)
    split_every : int or None
        If provided, split the table into multiple sub-tables every n columns (default: None)
    color_data : list of lists or None
        Optional 2D array to use for coloring cells (default: uses data_2d)
    label_data : list of lists or None
        Optional 2D array to use for cell text labels (default: uses data_2d)
    max_resids : int or None
        If provided, limit the number of columns to display to this value (default: None)
    
    Returns:
    --------
    fig, ax, table(s) : matplotlib figure, axes, and table object(s)
    """
    # Convert to numpy array for easier manipulation
    data_array = np.array(data_2d)
    n_rows, n_cols = data_array.shape
    
    # Limit columns if max_resids is specified
    if max_resids is not None and max_resids < n_cols:
        data_array = data_array[:, :max_resids]
        n_cols = max_resids
    
    # Use provided color_data or default to data_2d
    if color_data is None:
        color_array = data_array
    else:
        color_array = np.array(color_data)
        if max_resids is not None and max_resids < color_array.shape[1]:
            color_array = color_array[:, :max_resids]
    
    # Use provided label_data or default to data_2d
    if label_data is None:
        label_array = data_array
    else:
        label_array = np.array(label_data)
        if max_resids is not None and max_resids < label_array.shape[1]:
            label_array = label_array[:, :max_resids]
    
    # Normalize color values and get colors for each cell
    normalized = (color_array - vmin) / (vmax - vmin)
    
    # Use cmap library Colormap instead of plt.colormaps
    cm = Colormap(cmap)
    colors = cm(normalized)
    
    # If no split, create single table
    if split_every is None or split_every >= n_cols:
        fig, ax = plt.subplots(figsize=figsize)
        cell_colors = colors.tolist()
        # Override first row colors to white
        cell_colors[0] = [[1.0, 1.0, 1.0, 1.0]] * len(cell_colors[0])
        cell_text = label_array.tolist()
        
        the_table = ax.table(cellText=cell_text,
                             rowLabels=row_labels,
                             cellColours=cell_colors,
                             loc=loc,
                             cellLoc='center')
        
        # Apply edge styling, font size, and font family
        for cell in the_table.get_celld().values():
            cell.set_edgecolor(edge_color)
            cell.set_linewidth(edge_linewidth)
            cell.set_fontsize(font_size)
            cell.get_text().set_fontfamily('sans-serif')
        
        the_table.scale(1, row_height)
        ax.axis('off')
        plt.tight_layout()
        
        return fig, ax, the_table
    
    # Calculate number of sub-tables needed
    n_tables = int(np.ceil(n_cols / split_every))
    
    # Create figure with multiple subplots arranged vertically (one per row)
    fig, axes = plt.subplots(n_tables, 1, figsize=(figsize[0], figsize[1] * n_tables))
    if n_tables == 1:
        axes = [axes]
    
    # Remove spacing to give each subplot equal space
    plt.subplots_adjust(hspace=0, left=0, right=1, top=1, bottom=0)
    
    tables = []
    
    for i, ax in enumerate(axes):
        # Calculate column range for this sub-table
        col_start = i * split_every
        col_end = min((i + 1) * split_every, n_cols)
        n_cols_sub = col_end - col_start
        
        # Extract data, colors, and labels for this sub-table
        sub_colors = colors[:, col_start:col_end, :].tolist()
        # Override first row colors to white
        sub_colors[0] = [[1.0, 1.0, 1.0, 1.0]] * len(sub_colors[0])
        sub_labels = label_array[:, col_start:col_end].tolist()
        
        # Pad last table with empty columns if needed to match split_every width
        if n_cols_sub < split_every:
            padding_needed = split_every - n_cols_sub
            for row_idx in range(len(sub_labels)):
                sub_labels[row_idx].extend([''] * padding_needed)
                sub_colors[row_idx].extend([[1.0, 1.0, 1.0, 1.0]] * padding_needed)
            n_cols_sub = split_every
        
        # Create the sub-table without row labels
        the_table = ax.table(cellText=sub_labels,
                             cellColours=sub_colors,
                             bbox=[0, 0, 1, 1],
                             cellLoc='center')
        
        # Disable auto font size to maintain consistency
        the_table.auto_set_font_size(False)
        
        # Apply edge styling, font size, and font family
        for key, cell in the_table.get_celld().items():
            cell.set_edgecolor(edge_color)
            cell.set_linewidth(edge_linewidth)
            cell.set_fontsize(font_size)
            cell.get_text().set_fontfamily('sans-serif')
        
        ax.axis('off')
        
        tables.append(the_table)
    
    # Don't use tight_layout in split mode to maintain consistent sizing
    
    return fig, axes, tables

def create_barplot(vmin, vmax, output_path):
    """
    Create a barplot showing vmin and vmax values.
    
    Parameters:
    -----------
    vmin : float
        Minimum value for color normalization.
    vmax : float
        Maximum value for color normalization.
    output_path : str
        Path to save the barplot.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    
    labels = ['vmin', 'vmax']
    values = [vmin, vmax]
    colors = ['#3498db', '#e74c3c']
    
    bars = ax.bar(labels, values, color=colors, edgecolor='black', linewidth=1.5)
    
    # Add value labels on top of bars
    for bar, value in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{value:.2f}',
                ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    ax.set_ylabel('Value', fontsize=12, fontfamily='sans-serif')
    ax.set_title('Color Scale Range', fontsize=14, fontweight='bold', fontfamily='sans-serif')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Set y-axis to start from 0
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Barplot saved: {output_path}")

def get_spaced_lifetimes(dict_fasta,dict_lifetimes,dict_resids_and_spaces_fasta, show_every_resid=10):
    """
    Process lifetime data to match FASTA sequence alignment.
    Output order follows the input headers order (--headers argument).
    """
    nested_lifetimes_values = []
    headers = dict_lifetimes["headers"]
    lifetimes_spaced_nested = []
    sequences_spaced_nested = []
    resid_spaced_nested = []

    for i in range(len(headers)):
        header = dict_lifetimes["headers"][i]
        index_header = dict_resids_and_spaces_fasta["headers"].index(header)
        resids = dict_resids_and_spaces_fasta["sequences"][index_header]
        
        resnames = dict_fasta["sequences"][index_header]
        df_resids = list(dict_lifetimes["lifetimes"][i]["resid"])    
        df_lifetimes = list(dict_lifetimes["lifetimes"][i]["sum_ns"])    
        lifetimes_spaced = []
        resid_spaced = []
        
        # Find the last non-zero resid for this sequence
        last_resid = max([r for r in resids if r != 0]) if any(r != 0 for r in resids) else 0
        
        for e in range(len(resids)):
            resid = resids[e]
            # Show resid if: it's 1, divisible by show_every_resid, or the last resid
            if resid != 0 and (resid == 1 or resid % show_every_resid == 0 or resid == last_resid):
                resid_spaced.append(str(resid))
            else:
                resid_spaced.append('')
            if resid == 0:
                lifetimes_spaced.append(0)
            else:
                # Check if resid exists in df_resids, if not assign 0 (no lifetime data)
                if resid in df_resids:
                    lifetimes_spaced.append(df_lifetimes[df_resids.index(resid)])
                else:
                    lifetimes_spaced.append(0)
        lifetimes_spaced_nested.append(lifetimes_spaced)
        sequences_spaced_nested.append(resnames)
        resid_spaced_nested.append(resid_spaced)
    
    return lifetimes_spaced_nested, sequences_spaced_nested, resid_spaced_nested

def main():
    args = parse_arguments()

    dfs = csvs_to_dfs(args.files, headers=args.headers, log_transform=args.log_transform)

    # Transform vmin and vmax to log10 scale if log_transform is enabled
    if args.log_transform:
        args.vmin = np.log10(args.vmin + 1)
        args.vmax = np.log10(args.vmax + 1)

    dict_fasta, headers, sequences = fasta_to_dict(args.fasta)
    dict_resids_and_spaces_fasta = sequences_to_resids_with_spaces(dict_fasta)

    lifetimes_spaced_nested, sequences_spaced_nested, resid_spaced_nested = get_spaced_lifetimes(dict_fasta, dfs, dict_resids_and_spaces_fasta, show_every_resid=args.show_every_resid)

    # Calculate and print min/max values of the lifetime data
    lifetimes_array = np.array(lifetimes_spaced_nested)
    data_min = np.min(lifetimes_array)
    data_max = np.max(lifetimes_array)
    print(f"Data range: min={data_min:.2f}, max={data_max:.2f}")
    print(f"Color scale: vmin={args.vmin}, vmax={args.vmax}")

    # Prepend resid row to the data (first row will be residue IDs)
    # For color data, use zeros for the resid row (will be overridden to white anyway)
    zero_row = [[0] * len(resid_spaced_nested[0])]
    lifetimes_with_resid = zero_row + lifetimes_spaced_nested
    sequences_with_resid = [resid_spaced_nested[0]] + sequences_spaced_nested
    row_labels_with_resid = ['Resid'] + args.headers

    # Use different data for colors vs labels
    create_heatmap_table(
        row_labels=row_labels_with_resid, 
        data_2d=lifetimes_with_resid,  # Used as default
        color_data=lifetimes_with_resid,  # Data for colormap
        label_data=sequences_with_resid,  # Different data for cell text
        vmin=args.vmin, 
        vmax=args.vmax, 
        cmap=args.cmap,
        split_every=args.split_every,
        max_resids=args.max_resids,
        edge_color=args.edge_color,
        edge_linewidth=args.edge_linewidth,
        row_height=args.row_height,
        col_width=args.col_width,
        font_size=args.font_size,
        row_label_width=args.row_label_width
    )

    # Save figure
    plt.savefig(args.output, dpi=600)
    print(f"Output: {args.output}")
    
    # Create barplot with automatic naming if barplot_output is provided, or always if True
    if args.barplot_output:
        # Generate barplot filename from main output filename
        import os
        base_name = os.path.splitext(args.output)[0]
        ext = os.path.splitext(args.output)[1]
        barplot_path = f"{base_name}_barplot{ext}"
        create_barplot(args.vmin, args.vmax, barplot_path)

if __name__ == "__main__":
    main()