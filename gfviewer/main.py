from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
import math
import os
import argparse
from reportlab.lib.units import cm
import PyPDF2
from PyPDF2 import PdfReader, PdfWriter
from Bio.Graphics import BasicChromosome
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


color_code = {0: (round(230 / 255, 2), round(25 / 255, 2), round(75 / 255, 2)),
              1: (round(60 / 255, 2), round(180 / 255, 2), round(75 / 255, 2)),
              2: (round(255 / 255, 2), round(225 / 255, 2), round(25 / 255, 2)),
              3: (round(0 / 255, 2), round(130 / 255, 2), round(200 / 255, 2)),
              4: (round(245 / 255, 2), round(130 / 255, 2), round(48 / 255, 2)),
              5: (round(145 / 255, 2), round(30 / 255, 2), round(180 / 255, 2)),
              6: (round(70 / 255, 2), round(240 / 255, 2), round(240 / 255, 2)),
              7: (round(240 / 255, 2), round(50 / 255, 2), round(230 / 255, 2)),
              8: (round(210 / 255, 2), round(245 / 255, 2), round(60 / 255, 2)),
              9: (round(250 / 255, 2), round(190 / 255, 2), round(212 / 255, 2)),
              10: (round(0 / 255, 2), round(128 / 255, 2), round(128 / 255, 2)),
              11: (round(220 / 255, 2), round(190 / 255, 2), round(255 / 255, 2)),
              12: (round(170 / 255, 2), round(110 / 255, 2), round(40 / 255, 2)),
              13: (round(255 / 255, 2), round(250 / 255, 2), round(200 / 255, 2)),
              14: (round(128 / 255, 2), round(0 / 255, 2), round(0 / 255, 2)),
              15: (round(170 / 255, 2), round(255 / 255, 2), round(195 / 255, 2)),
              16: (round(128 / 255, 2), round(128 / 255, 2), round(0 / 255, 2)),
              17: (round(255 / 255, 2), round(215 / 255, 2), round(180 / 255, 2)),
              18: (round(0 / 255, 2), round(0 / 255, 2), round(128 / 255, 2)),
              19: (round(128 / 255, 2), round(128 / 255, 2), round(128 / 255, 2)),
              20: (round(0 / 255, 2), round(0 / 255, 2), round(0 / 255, 2))}

color_index = {0: 'red',
              1: 'green',
              2: 'yellow',
              3: 'blue',
              4: 'orange',
              5: 'purple',
              6: 'cyan',
              7: 'magenta',
              8: 'lime',
              9: 'pink',
              10: 'teal',
              11: 'lavender',
              12: 'brown',
              13: 'beige',
              14: 'maroon',
              15: 'mint',
              16: 'olive',
              17: 'apricot',
              18: 'navy',
              19: 'grey',
              20: 'black'}

max_gfs = 19
centro_color = 19


def print_color_guide():
    """Prints the default color scheme."""
    print("Default Color Scheme for Plotting:\n")
    for index, color in color_index.items():
        print(f"  {index + 1}: {color}")
    print("\nUse this as a reference for custom color assignments.")


def get_fasta_sequence_lengths(fasta_file):
    """Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequence lengths as values."""
    sequence_lengths = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_lengths[record.id] = len(record.seq)

    return sequence_lengths


def get_sequence_lengths(in_file):
    seq_lens = {}

    with open(in_file, "r") as file:
        for line in file:
            tokens = line.strip().split(',')
            seq_lens[tokens[0]] = int(tokens[1])

    return seq_lens


def get_color_maps(in_file):
    color_map = {}

    with open(in_file, "r") as file:
        for line in file:
            tokens = line.strip().split(',')
            if len(tokens) not in [2, 4]:
                print(
                    f"Error: Color file is in invalid format. Requires <gf_id>,<color_code> or <gf_id>,<[0,1]>,<[0,1]>,<[0,1]> per line")
                exit(1)
            elif len(tokens) == 2:
                color_map[tokens[0]] = color_code[int(tokens[1]) - 1]
            elif len(tokens) == 4:
                color_map[tokens[0]] = (float(tokens[1]), float(tokens[2]), float(tokens[3]))

    return color_map


def detect_and_read_file(file_path):
    """Detects the file type (CSV, TSV, XLSX) and reads it into a pandas DataFrame."""
    _, file_extension = os.path.splitext(file_path)

    if file_extension.lower() == ".csv":
        df = pd.read_csv(file_path)
    elif file_extension.lower() == ".tsv":
        df = pd.read_csv(file_path, sep="\t")
    elif file_extension.lower() == ".xlsx":
        df = pd.read_excel(file_path, sheet_name=0)  # Reads the first sheet by default
    else:
        raise ValueError("Unsupported file format. Please provide a CSV, TSV, or XLSX file.")

    return df


def convert_to_gb(data, gb_path, seq_lengths):
    records = {}
    gfs = []

    for seq in seq_lengths:
        records[seq] = SeqRecord(Seq("N" * seq_lengths[seq]), id=seq, annotations={"molecule_type": "DNA"})

    for _, row in data.iterrows():
        if str(row["gene_family"]) == "centromere" and str(row["strand"]) == "0":
            feature_1 = SeqFeature(FeatureLocation(row["start"], row["end"]), type=str(row["gene_family"]), strand=1)
            feature_1.qualifiers["gene"] = [str(row["gene_id"]).capitalize()]
            records[str(row["chromosome"])].features.append(feature_1)

            feature_2 = SeqFeature(FeatureLocation(row["start"], row["end"]), type=str(row["gene_family"]), strand=-1)
            feature_2.qualifiers["gene"] = [str(row["gene_id"]).capitalize()]
            records[str(row["chromosome"])].features.append(feature_2)
        else:
            feature = SeqFeature(FeatureLocation(row["start"], row["end"]), type=str(row["gene_family"]),
                                 strand=1 if row["strand"] == "+" else -1)

            # Add gene qualifiers (ID and Family)
            feature.qualifiers["gene"] = [str(row["gene_id"])]  # Storing gene ID

            records[str(row["chromosome"])].features.append(feature)

            if str(row["gene_family"]) not in gfs:
                gfs.append(str(row["gene_family"]))

    for seq in records:
        out_file = gb_path + "/" + seq + ".gb"
        with open(out_file, "w") as out:
            SeqIO.write(records[seq], out, "genbank")

        print(f"GB file successfully created: {out_file}")

    return gfs


def combine_legend_with_plot(legend_file, gf_file):
    # Open both PDFs
    legend = PdfReader(legend_file)
    gf = PdfReader(gf_file)
    writer = PdfWriter()

    # Overlay the legend PDF onto the first page
    legend_page = legend.pages[0]  # Modify first page
    gf_page = gf.pages[0]
    legend_page.merge_page(gf_page)
    writer.add_page(legend_page)

    # Save the final output
    with open(gf_file, "wb") as output:
        writer.write(output)


def merge_pdfs(out_path, total_pages):
    merger = PyPDF2.PdfMerger()

    # Iterate through input PDFs and add their pages to the merger
    for p in range(1, total_pages + 1):
        with open(out_path + "/" + str(p) + ".pdf", "rb") as pdf:
            merger.append(pdf)

    # Write the merged PDF to the output file
    with open(out_path + "/mgf.pdf", "wb") as output:
        merger.write(output)

    for p in range(1, total_pages + 1):
        os.remove(out_path + "/" + str(p) + ".pdf")


def plot_legends(color_map, out_file, location, orientation, legends_per_page, centromeres, nrows):
    if legends_per_page:
        fig, ax = plt.subplots(figsize=(11.69, 8.27))
    else:
        fig, ax = plt.subplots()
    color_list = []
    gfs = list(color_map.keys())
    for gf in gfs:
        patch = mpatches.Patch(color=color_map[gf], label=gf)
        color_list.append(patch)
    if centromeres:
        centromere = mpatches.Patch(color=color_code[centro_color], label='Centromere')  # grey
        color_list.append(centromere)
    if legends_per_page:
        if location == 'lower':
            ax.legend(handles=color_list, loc="lower center", ncol=int(math.ceil(len(gfs) / nrows)))
        elif location == 'upper':
            ax.legend(handles=color_list, loc="upper center", ncol=int(math.ceil(len(gfs) / nrows)))
        elif location == 'left':
            ax.legend(handles=color_list, loc="center left")
        elif location == 'right':
            ax.legend(handles=color_list, loc="center right")
    else:
        if orientation == 'horizontal':
            ax.legend(handles=color_list, loc="center", ncol=int(math.ceil(len(gfs) / nrows)))
        elif orientation == 'vertical':
            ax.legend(handles=color_list, loc="center")

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")


def plot_gene_families(args):
    plot_path = args.output_directory + "/plot"
    os.makedirs(plot_path, exist_ok=True)

    _, file_extension = os.path.splitext(args.genome_file)
    if file_extension.lower() in [".fasta", ".fa", ".fna"]:
        seq_lengths = get_fasta_sequence_lengths(args.genome_file)
    elif file_extension.lower() == ".txt":
        seq_lengths = get_sequence_lengths(args.genome_file)
    else:
        raise ValueError("Unsupported file format. Please provide a FASTA or TXT file.")

    chrs = list(seq_lengths.keys())
    max_len = seq_lengths[max(seq_lengths, key=seq_lengths.get)]  # Could compute this from the entries dict

    df = detect_and_read_file(args.data_file)
    gb_path = args.output_directory + "/tmp"
    os.makedirs(gb_path, exist_ok=True)
    gfs = convert_to_gb(df, gb_path, seq_lengths)

    if args.color_map_file:
        color_map = get_color_maps(args.color_map_file)
        gfs = list(color_map.keys())
    else:
        color_map = {}
        index = 0
        for gf in gfs:
            color_map[gf] = color_code[index]
            index += 1

    if len(gfs) > max_gfs:
        print(f"Error: Number of gene families exceeds limit.")
        exit(1)

    legend_file = plot_path + "/legends.pdf"
    plot_legends(color_map, legend_file, args.legend_location, args.legend_orientation, args.legends_per_page,
                 args.centromeres, args.number_of_rows_in_legends)

    current_chr = 0
    total_pages = int(math.ceil(len(chrs) / args.number_of_chromosomes_per_page))
    for p in range(1, total_pages + 1):
        chr_diagram = BasicChromosome.Organism()
        chr_diagram.page_size = (29.693 * cm, 21.006 * cm)  # A4 landscape

        current_iteration = args.number_of_chromosomes_per_page
        if p == total_pages and len(chrs) % args.number_of_chromosomes_per_page != 0:
            current_iteration = len(chrs) % args.number_of_chromosomes_per_page
        current_range = current_chr + current_iteration
        while current_chr < current_range:
            ch = chrs[current_chr]
            input_gb_file = gb_path + "/" + ch + ".gb"
            record = SeqIO.read(input_gb_file, "genbank")
            length = len(record)
            features = []
            for f in record.features:
                if f.type in gfs:
                    f.qualifiers["color"] = [color_map[f.type]]
                    features.append(f)
                elif f.type == "centromere" and args.centromeres:
                    f.qualifiers["color"] = [color_code[centro_color]]
                    features.append(f)

            cur_chromosome = BasicChromosome.Chromosome(ch.capitalize())
            # Set the scale to the MAXIMUM length plus the two telomeres in bp,
            # want the same scale used on all chromosomes so they can be
            # compared to each other
            cur_chromosome.scale_num = max_len + 2 * args.telomere_length

            # Add an opening telomere
            start = BasicChromosome.TelomereSegment()
            start.scale = args.telomere_length
            cur_chromosome.add(start)

            # Add a body - again using bp as the scale length here.
            body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
            body.scale = length
            cur_chromosome.add(body)

            # Add a closing telomere
            end = BasicChromosome.TelomereSegment(inverted=True)
            end.scale = args.telomere_length
            cur_chromosome.add(end)

            # This chromosome is done
            chr_diagram.add(cur_chromosome)
            current_chr += 1

        current_file = plot_path + "/" + str(p) + ".pdf"
        chr_diagram.draw(current_file, "")

        if args.legends_per_page:
            combine_legend_with_plot(legend_file, current_file)

    if args.concatenate_pages:
        merge_pdfs(plot_path, total_pages)

    if args.legends_per_page:
        os.remove(legend_file)

    print(f"Multi-gene families are successfully plotted inside: {plot_path}")


def parse_arguments():
    # Create a parser
    parser = argparse.ArgumentParser(description="Process input arguments for the GFViewer tool.")

    # Define arguments for the user parameters
    parser.add_argument('-d', '--data_file', type=str, required=True,
                        help='A file (.xlsx/.csv/.tsv) containing gene family and location information for each gene')
    parser.add_argument('-g', '--genome_file', type=str, required=True,
                        help='A fasta (.fasta/.fna/.fa) file containing the genome sequence or a text (.txt) file containing the chromosome ids with their lengths; <seq_id>,<seq_length> per line')
    parser.add_argument('-o', '--output_directory', type=str, required=True,
                        help='Path to the output directory')
    parser.add_argument('-c', '--color_map_file', type=str,
                        help='A text (.txt) file containing the gene family (gf) ids with their color codes; <gf_id>,<color_code> or <gf_id>,<[0,1]>,<[0,1]>,<[0,1]> per line')
    parser.add_argument('-l', '--legend_location', type=str, default='lower',
                        help='Specify the location (upper/lower/left/right) of legends only when adding them to each page (default: lower)')
    parser.add_argument('-or', '--legend_orientation', type=str, default='horizontal',
                        help='Specify the orientation (horizontal/vertical) of legends only when plotting legends separately (default: horizontal)')
    parser.add_argument('-t', '--telomere_length', type=int, default=10000,
                        help='The length of telomeres in bp used in the plot (default: 10000)')
    parser.add_argument('-p', '--number_of_chromosomes_per_page', type=int, default=3,
                        help='Number of chromosomes to be plotted per page (default: 3)')
    parser.add_argument('-r', '--number_of_rows_in_legends', type=int, default=2,
                        help='Number of rows in the legends (default: 2)')

    # Define optional flags
    parser.add_argument('-cen', '--centromeres', action='store_true',
                        help='Plot centromeres of the chromosomes along with multi gene families')
    parser.add_argument('-lpp', '--legends_per_page', action='store_true',
                        help='Plot legends per page in the PDF')
    parser.add_argument('-conc', '--concatenate_pages', action='store_true',
                        help='Concatenate the pages into a single PDF file')

    # Parse the arguments
    args = parser.parse_args()

    # Ensure the required files exist and output directory is valid
    if not args.data_file:
        print(f"Error: Data file is required.")
        exit(1)
    if not args.genome_file:
        print(f"Error: Genome file is required.")
        exit(1)
    if not args.output_directory:
        print(f"Error: Output directory is required.")
        exit(1)

    if not os.path.isfile(args.data_file):
        print(f"Error: Data file '{args.data_file}' does not exist.")
        exit(1)

    if not os.path.isfile(args.genome_file):
        print(f"Error: Genome file '{args.genome_file}' does not exist.")
        exit(1)

    if args.color_map_file and not os.path.isfile(args.color_map_file):
        print(f"Warning: Color map file '{args.color_map_file}' does not exist. Using default colors.")

    if args.legend_location and args.legend_location not in ['upper', 'lower', 'left', 'right']:
        print(f"Warning: Invalid legends location '{args.legend_location}'. Using default location.")
        args.legend_location = 'lower'

    if args.legend_orientation and args.legend_orientation not in ['horizontal', 'vertical']:
        print(f"Warning: Invalid legends orientation '{args.legend_orientation}'. Using default orientation.")
        args.legend_orientation = 'horizontal'

    return args


def main():
    # Parse the arguments
    args = parse_arguments()

    # Here you can access each argument, for example:
    print("Data File:", args.data_file)
    print("Genome File:", args.genome_file)
    print("Output Directory:", args.output_directory)
    print("Telomere Length:", args.telomere_length)
    print("Chromosomes Per Page:", args.number_of_chromosomes_per_page)
    print("Legends location (only when adding to each page):", args.legend_location)
    print("Legends orientation (only when plotting legends separately):", args.legend_orientation)
    print("Rows in Legends (only when upper/lower location or horizontal orientation):", args.number_of_rows_in_legends)

    # Create output directory and the required sub-directories
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
        print(f"Output directory '{args.output_directory}' is created.")

    # Process flags
    if args.centromeres:
        print("Centromeres of the chromosomes will be plotted along with multi-gene families")

    if args.legends_per_page:
        print("Legends will be plotted per page.")

    if args.concatenate_pages:
        print("Pages will be concatenated into a single PDF.")

    plot_gene_families(args)


if __name__ == "__main__":
    main()
