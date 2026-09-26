#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import os
from dataclasses import dataclass, field
import seaborn as sns
from pathlib import Path

#whereami
ROOT = Path(__file__).resolve().parent

#handle input
def parse_args(argv):
    print("Parsing input...")
    parser = argparse.ArgumentParser(description="Generate plots from compiled host element results")
    parser.add_argument("input_file",metavar="INPUT_FILE", help="Path to Edwards Excel file, or a another excel file containing indexed elements and their proportions")
    parser.add_argument("output_dir", metavar="OUTPUT_DIR", type=Path, help="Directory where plots will be saved")

    #optionals for MGE heatmap (defaults follow Edwards design)
    parser.add_argument("--sheet_mge", metavar="SHEET_NAME", default="Host_Element_Prevalence", type=str, required=False, help="Name of the sheet in the input Excel file containing MGE data")
    parser.add_argument("--mge_cols", metavar="MGE_COLS_IDENTIFIER", default="__Proportional_Called", type=str, required=False, help="Identifier to identify all MGE columns in the input file")
    parser.add_argument("--index_id_mge", metavar="INDEX_IDENTIFIER", default="Element", type=str, required=False, help="Identifier for the MGE index column in the input file")

    #optionals for gene heatmap (defaults follow Edwards design)
    parser.add_argument("--sheet_gene", metavar="SHEET_NAME", default="Proportional_Data", type=str, required=False, help="Name of the sheet in the input Excel file containing gene data")
    parser.add_argument("--gene_cols_gene", metavar="GENE_COLS_IDENTIFIER", default="__Proportional_Called", type=str, required=False, help="Identifier to identify all gene columns in the input file")
    parser.add_argument("--index_id_gene", metavar="INDEX_IDENTIFIER", default="Element_Gene", type=str, required=False, help="Identifier for the gene index column in the input file")
    return parser.parse_args(argv)

# store stuff to modify heatmaps
@dataclass(frozen=True)
class HeatmapConfig:
    # figure and axis configuration
    figsize: tuple[float, float]    = (16, 12)
    title: str                      = "Heatmap"
    xlabel: str                     = "Host"
    ylabel: str                     = "MGE"

    # heatmap specific configuration
    vmin: float                     = 0.0
    vmax: float                     = 1.0
    cmap: str                       = "Reds"
    annot: bool                     = True
    fmt: str                        = ".2f"
    cbar_kws: dict[str, str]        = field(default_factory=dict)
    linewidths: float               = 0.5

# read and store the mmseqs data generated from Edwards script (reads all sheets from the Excel file)
@dataclass  
class mmseqs2Results:
    mmseqs_sheets: dict[str, pd.DataFrame] = field(default_factory=dict)

    @classmethod
    def from_excel(cls, input_path: str, **kwargs):
        try:
            #reads ALL sheets from Edwards Excel file (by name)
            sheets = pd.read_excel(input_path, sheet_name=None, **kwargs)
            return cls(mmseqs_sheets=sheets)
        except ValueError as e:
            raise ValueError(f"Error reading Excel file {input_path!r}: {e}")

    def get_named_sheet(self, sheet_name: str) -> pd.DataFrame:
        try:
            return self.mmseqs_sheets[sheet_name]
        except KeyError as e:
            raise KeyError(f"Error: sheet {sheet_name} does not exists. {e}")

#prepare MGE data and write plots
class HostElementPlotter:
    def __init__(self, mmseqs_data_sheet: pd.DataFrame, config: HeatmapConfig, output_filename: str = "MGE_heatmap.png"):
        self.mmseqs_data_sheet: pd.DataFrame = mmseqs_data_sheet
        self.config = config
        self.output_filename = output_filename
    
    #prepare heatmap data for MGEs
    def prepare_MGE_heatmap_data(self,  mmseqs_data_sheet: pd.DataFrame, 
                                        col_identifier: str, 
                                        index_col_header: str,
                                        clean_col_identifier: bool = True) -> pd.DataFrame:
        #prepares heatmap data:
        # Cols: columns that contain the col_identifier (proportions) 
        # Index: set the index to the specified column header (MGE names)
        proportional_cols: list[str] = [col for col in mmseqs_data_sheet.columns if col_identifier in col]
        heatmap_data = pd.DataFrame(mmseqs_data_sheet[proportional_cols])
        if clean_col_identifier:
            heatmap_data.columns = [col.replace(col_identifier, "") for col in proportional_cols]
        heatmap_data.index = mmseqs_data_sheet[index_col_header]
        return heatmap_data

    def write_MGE_heatmap(self, heatmap_data: pd.DataFrame, output_dir: Path):
        config = self.config

        figure = plt.figure(figsize=config.figsize)
        axis_obj = sns.heatmap(heatmap_data, 
                                vmin=config.vmin,
                                vmax=config.vmax,
                                cmap=config.cmap, 
                                annot=config.annot, 
                                fmt=config.fmt, 
                                linewidths=config.linewidths, 
                                cbar_kws=config.cbar_kws)
        axis_obj.set_title(config.title, fontsize=16)
        axis_obj.set_xlabel(config.xlabel, fontsize=14)
        axis_obj.set_ylabel(config.ylabel, fontsize=14)
        axis_obj.set_xticklabels(axis_obj.get_xticklabels(), rotation=0, fontsize=16)
        axis_obj.set_yticklabels(axis_obj.get_yticklabels(), rotation=0, fontsize=16)
        axis_obj.figure.tight_layout()
        output_dir.mkdir(parents=True, exist_ok=True)
        axis_obj.figure.savefig(str(output_dir / self.output_filename))
        plt.close(figure)

#prepare gene data and write plots
class GenePlotter:
    def __init__(self, mmseqs_data_sheet: pd.DataFrame, config: HeatmapConfig, output_filename: str = "gene_heatmap.png"):
        self.mmseqs_data_sheet: pd.DataFrame = mmseqs_data_sheet
        self.config = config
        self.output_filename = output_filename

    # prepare heatmap data for genes
    def prepare_gene_heatmap_data(self,
                                  mmseqs_data_sheet: pd.DataFrame,
                                  col_identifier: str,
                                  index_col_header: str,
                                  clean_col_identifier: bool = True) -> pd.DataFrame:
        proportional_cols: list[str] = [col for col in mmseqs_data_sheet.columns if col_identifier in col]
        heatmap_data = mmseqs_data_sheet[proportional_cols].copy()
        if clean_col_identifier:
            heatmap_data.columns = [col.replace(col_identifier, "") for col in proportional_cols]
        heatmap_data.index = mmseqs_data_sheet[index_col_header]
        return heatmap_data

    # reads data and groups genes, based on index prefix
    def get_gene_groups(self, 
                        index_col_header: str,
                        sep: str = "_") -> dict[str, pd.DataFrame]:
        
        gene_groups: dict[str, pd.DataFrame] = {}
        unique_prefixes = set()

        #get unique prefixes from the index column
        for row in self.mmseqs_data_sheet[index_col_header]:
            prefix = row.split(sep)[0]
            unique_prefixes.add(prefix + sep)
        print(f"{len(unique_prefixes)} Unique prefixes: {unique_prefixes}")

        #group rows by unique prefix (each prefix gets its own DataFrame)
        for prefix in unique_prefixes:
            mask_vector = self.mmseqs_data_sheet[index_col_header].str.startswith(prefix)
            gene_groups[prefix] = self.mmseqs_data_sheet.loc[mask_vector]
            print(f"Prefix: {prefix}, Genes: {len(gene_groups[prefix])}")

        return gene_groups

    def write_gene_plot(self, 
                        gene_data: pd.DataFrame, 
                        output_dir: Path, 
                        output_filename: str):
        output_path = output_dir / output_filename
        config = self.config
        

        figure = plt.figure(figsize=config.figsize)
        axis_obj = sns.heatmap(gene_data, 
                                vmin=config.vmin,
                                vmax=config.vmax,
                                cmap=config.cmap, 
                                annot=config.annot, 
                                fmt=config.fmt, 
                                linewidths=config.linewidths, 
                                cbar_kws=config.cbar_kws)
        axis_obj.set_title(config.title, fontsize=16)
        axis_obj.set_xlabel(config.xlabel, fontsize=14)
        axis_obj.set_ylabel(config.ylabel, fontsize=14)
        axis_obj.set_xticklabels(axis_obj.get_xticklabels(), rotation=0, fontsize=16)
        axis_obj.set_yticklabels(axis_obj.get_yticklabels(), rotation=0, fontsize=16)
        axis_obj.figure.tight_layout()
        output_dir.mkdir(parents=True, exist_ok=True)
        axis_obj.figure.savefig(str(output_path))
        plt.close(figure)

def main():
    #args
    args = parse_args(sys.argv[1:])

    #get data
    print(f"Reading: {args.input_file}")
    mmseqs2_data = mmseqs2Results.from_excel(args.input_file)

    #________make MGE heatmap________
    mge_proportions = mmseqs2_data.get_named_sheet(args.sheet_mge)
    MGE_config = HeatmapConfig(title="MGE Heatmap", 
                               ylabel="MGEs", 
                               cbar_kws={"label": "MGE proportion"})
    plotter = HostElementPlotter(mmseqs_data_sheet=mge_proportions, 
                                 config=MGE_config)
    MGE_heatmap_data = plotter.prepare_MGE_heatmap_data(mge_proportions, 
                                                        args.mge_cols, 
                                                        args.index_id_mge)

    print(
        f"Generating MGE heatmap ({MGE_heatmap_data.shape[0]} elements, "
        f"{MGE_heatmap_data.shape[1]} hosts)..."
    )

    plotter.write_MGE_heatmap(heatmap_data=MGE_heatmap_data, 
                              output_dir=args.output_dir)
    print(f"Saved MGE heatmap: {args.output_dir / plotter.output_filename}")

    #________Make gene heatmaps________
    print(f"Preparing gene heatmaps '{args.sheet_gene}'...")

    gene_proportions = mmseqs2_data.get_named_sheet(args.sheet_gene)
    GENE_config = HeatmapConfig(title="Gene Heatmap", 
                                ylabel="Genes", 
                                cbar_kws={"label": "Gene proportion"})
    gene_plotter = GenePlotter(mmseqs_data_sheet=gene_proportions, 
                               config=GENE_config)
    gene_groups = gene_plotter.get_gene_groups(index_col_header=args.index_id_gene)

    #create output directory for gene heatmaps
    gene_specific_output_dir = args.output_dir / "gene_specific"
    gene_specific_output_dir.mkdir(parents=True, exist_ok=True)

    # fenerate and save gene heatmaps for each gene group
    for prefix, gene_data in gene_groups.items():
        output_filename = f"{prefix}{gene_plotter.output_filename}"
        print(f"Generating gene heatmap for '{prefix}' ({len(gene_data)} genes)...")
        gene_heatmap_data = gene_plotter.prepare_gene_heatmap_data(gene_data, 
                                                                   args.gene_cols_gene, 
                                                                   args.index_id_gene)
        gene_plotter.write_gene_plot(gene_data=gene_heatmap_data, 
                                     output_dir=gene_specific_output_dir,
                                     output_filename=output_filename)
    print(f"Saved gene heatmapa: {gene_specific_output_dir}")

if __name__ == "__main__":
    main()