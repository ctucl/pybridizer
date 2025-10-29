"""
Pybridizer - A tool for designing HCR (Hybridization Chain Reaction) probes.

This module provides functionality to:
1. Process FASTA sequences
2. Create oligonucleotide pairs
3. Perform BLAST searches
4. Filter and optimize probe pairs
5. Add HCR hairpin sequences

Dependencies:
- Biopython
- pandas
- biothings_client
- BLAST+ command line tools

Author: Chintan Trivedi, UCL

"""

import sys
from subprocess import Popen, PIPE
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, MeltingTemp
from Bio import Align
from Bio.Seq import Seq
from biothings_client import get_client
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def read_fasta(file_path):
    """
    Reads the first record from a FASTA file with error handling.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        tuple: (ID, description, sequence) if successful.
               (None, None, None) if an error occurs.
    """
    try:
        records = list(SeqIO.parse(file_path, 'fasta'))

        if not records:
            print(f"Error: No FASTA records found in '{file_path}'. The file might be empty.", 
                  file=sys.stderr)
            return None, None, None

        first_record = records[0]
        ID = first_record.id
        desc = first_record.description
        sequence = first_record.seq

        print('Transcript ID: ', ID)
        print('Transcript Description: ', desc)
        print('Transcript Sequence Length: ', len(sequence), 'nt')

        return ID, desc, sequence

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.", file=sys.stderr)
        return None, None, None
    except ValueError as e:
        print(f"Error: Could not parse '{file_path}' as FASTA. Details: {e}", file=sys.stderr)
        return None, None, None
    except Exception as e:
        print(f"An unexpected error occurred while processing '{file_path}': {e}", 
              file=sys.stderr)
        return None, None, None


def create_oligos(sequence, oligo_length=25, gap_length=2, frame_start_position=0):
    """
    Designs a set of forward and reverse oligos from an input sequence.

    Args:
        sequence (str): The input nucleotide sequence.
        oligo_length (int, optional): The length of each oligo. Defaults to 25.
        gap_length (int, optional): The gap between oligos. Defaults to 2.
        frame_start_position (int, optional): Starting position. Defaults to 0.

    Returns:
        list: Generated oligos as SeqRecord objects.
    """
    if not all(isinstance(x, int) for x in [oligo_length, gap_length, frame_start_position]):
        raise TypeError("Oligo length, gap length, and frame start position must be integers.")

    if oligo_length <= 0:
        raise ValueError("Oligo length must be positive.")
    if gap_length < 0:
        raise ValueError("Gap length must be non-negative.")
    if frame_start_position < 0 or frame_start_position >= len(sequence):
        raise ValueError("Frame start position is out of range.")

    oligos_all = []
    n = frame_start_position
    k = frame_start_position + oligo_length + gap_length
    j = len(sequence)

    while (k + oligo_length - gap_length) <= j:
        e_rna = sequence[n:n + oligo_length]
        o_rna = sequence[k:k + oligo_length]
        EV = e_rna.reverse_complement()
        OD = o_rna.reverse_complement()
        oligos_all.extend([EV, OD])
        n += (oligo_length + gap_length) * 2
        k += (oligo_length + gap_length) * 2

    print('Oligos tiled along the transcript sequence:')
    print(oligos_all)
    return oligos_all


def blast_oligos(oligos_all, dbname, taxid=None):
    """
    Performs BLAST search for oligonucleotides.

    Args:
        oligos_all (list): List of oligonucleotide sequences
        dbname (str): BLAST database name
        taxid (str/int, optional): Taxonomy ID for filtering. If None, no taxid filter is applied.

    Returns:
        pandas.DataFrame: Filtered BLAST results
    """
    if not oligos_all:
        raise ValueError("No oligonucleotides provided for BLAST search")

    records = []
    for (index, sequence) in enumerate(oligos_all):
        records.append(SeqRecord(sequence, id=str(index+1)))
    
    try:
        SeqIO.write(records, 'oligos.faa', 'fasta')
    except IOError as e:
        raise IOError(f"Failed to write temporary FASTA file: {e}")

    # Base BLAST command
    cmd_list = [
        'blastn',
        '-db', dbname,
        '-query', 'oligos.faa',
        '-outfmt', '6',
        '-out', 'result_tab.txt',
        '-task', 'blastn-short',
        '-evalue', '100',
        '-strand', 'minus'
    ]

    # Add taxid parameter only if provided
    if taxid is not None:
        cmd_list.extend(['-taxids', str(taxid)])

    print('Running batch BLAST on oligos...')
    try:
        blast_process = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        stdout, stderr = blast_process.communicate()
        return_code = blast_process.returncode

        if return_code != 0:
            raise IOError(f"BLAST Error (code {return_code}): {stderr.decode()}")
        
        print('BLAST completed successfully')
        
        blast_result = pd.read_csv(
            'result_tab.txt', 
            sep="\t", 
            header=None,
            names=[
                'qid', 'sseqid', 'pident', 'length', 'mismatch',
                'gapopen', 'qstart', 'qend', 'sstart', 'send',
                'evalue', 'bitscore'
            ]
        )

        blast_result = blast_result[
            ~blast_result['sseqid'].str.startswith(('XM', 'XR')) &
            (blast_result['length'] > 0.6 * blast_result['length'].max())
        ]

        return blast_result

    except FileNotFoundError:
        raise FileNotFoundError("BLAST executable not found. Please ensure BLAST+ is installed.")
    except Exception as e:
        raise IOError(f"BLAST execution failed: {str(e)}")



def fetch_geneIDs(blast_result, taxid='zebrafish'):
    """
    Fetches gene IDs for BLAST results using biothings client.

    Args:
        blast_result (pandas.DataFrame): BLAST results
        taxid (str): Taxonomy identifier (default: 'zebrafish')

    Returns:
        pandas.DataFrame: BLAST results with gene IDs
    """
    if blast_result.empty:
        raise ValueError("Empty BLAST results provided")
    
    if 'sseqid' not in blast_result.columns:
        raise ValueError("Invalid BLAST results format: missing 'sseqid' column")

    IDs = blast_result['sseqid'].tolist()
    IDs = ','.join(IDs)

    try:
        gene_client = get_client('gene')
        gids = gene_client.querymany(
            IDs,
            scopes='refseq.rna',
            fields='entrezgene',
            species=taxid,
            returnall=True
        )

        if not gids['out']:
            raise ValueError(f"No gene IDs found for taxonomy {taxid}")

        genes = gids['out']
        df = pd.DataFrame(genes)
        df = df.set_index(blast_result.index)
        geneid = df['entrezgene']
        
        blast_data = blast_result.copy()
        blast_data.insert(1, "geneid", geneid, allow_duplicates=True)
        
        return blast_data

    except Exception as e:
        raise ValueError(f"Failed to fetch gene IDs: {str(e)}")


def adjacency_filter(dframe):
    """
    Filters dataframe to keep only adjacent query IDs.

    Args:
        dframe (pandas.DataFrame): Input dataframe with 'qid' column

    Returns:
        pandas.DataFrame: Filtered dataframe
    """
    dataset = dframe.copy()
    uqids = dataset.qid.unique()
    
    j = -1 
    for i in range(len(uqids)-1):
        if not (uqids[i]-j==1 or uqids[i+1]-uqids[i]==1):
            indexNames = dataset[dataset['qid'] == uqids[i]].index
            dataset.drop(indexNames, inplace=True)
        j = uqids[i]

    if uqids[-1]-uqids[-2] != 1:
        indexNames = dataset[dataset['qid'] == uqids[-1]].index
        dataset.drop(indexNames, inplace=True)
            
    return dataset



def filter_and_rank(oligos_all, blast_data, GC_range=[0.37, 0.85], Tm_range=[47, 85]):
    """
    Filters and ranks oligo pairs based on GC content and melting temperature.

    Args:
        oligos_all (list): List of oligonucleotide sequences
        blast_data (pandas.DataFrame): BLAST results (can contain either 'geneid' or 'sseqid' column)
        GC_range (list): Acceptable GC content range [min, max]
        Tm_range (list): Acceptable melting temperature range [min, max]

    Returns:
        pandas.DataFrame: Filtered and ranked probe pairs
    """
    if not isinstance(GC_range, list) or len(GC_range) != 2:
        raise ValueError("GC_range must be a list of [min, max] values")
    if not isinstance(Tm_range, list) or len(Tm_range) != 2:
        raise ValueError("Tm_range must be a list of [min, max] values")
    if GC_range[0] >= GC_range[1] or Tm_range[0] >= Tm_range[1]:
        raise ValueError("Range minimum must be less than maximum")

    # Determine which column to use for matching
    id_column = 'geneid' if 'geneid' in blast_data.columns else 'sseqid'
    
    data_new = blast_data.copy()
    print('Filtering oligos based on GC content and Melting Temperature...')

    # Filter by GC content and melting temperature
    GCcon = [gc_fraction(oligo) for oligo in oligos_all]
    MT = [MeltingTemp.Tm_GC(oligo, strict=False) for oligo in oligos_all]

    for idx, (gc, mt) in enumerate(zip(GCcon, MT)):
        if ((gc <= GC_range[0] or gc >= GC_range[1]) or 
            (mt <= Tm_range[0] or mt >= Tm_range[1])):
            indexNames = data_new[data_new['qid'] == idx+1].index
            data_new.drop(indexNames, inplace=True)

    # Get unique query IDs after initial filtering
    uqids = sorted(data_new.qid.unique())
    
    print('Filtering adjacent oligos with same transcript hits...')
    valid_pairs = []
    i = 0
    
    while i < len(uqids) - 1:
        current_id = uqids[i]
        next_id = uqids[i + 1]
        
        # Check if IDs are consecutive and represent a valid pair (even/odd)
        if next_id - current_id == 1 and current_id % 2 == 1:
            # Get hits for both oligos using the determined column
            hits1 = set(data_new.loc[data_new['qid'] == current_id, id_column])
            hits2 = set(data_new.loc[data_new['qid'] == next_id, id_column])
            
            # Check for shared off-target hits
            if len(hits1 & hits2) <= 1:  # Allow only target hit overlap
                valid_pairs.append((current_id, next_id))
            i += 2
        else:
            i += 1
    
    # Create final datasheet with valid pairs only
    pairs_data = []
    for pair in valid_pairs:
        oligo1_pos, oligo2_pos = pair
        score = ((len(data_new.loc[data_new['qid'] == oligo1_pos, id_column].unique()) +
                 len(data_new.loc[data_new['qid'] == oligo2_pos, id_column].unique())) / 2)
        
        pairs_data.append({
            'Oligo1_Position': oligo1_pos,
            'Oligo2_Position': oligo2_pos,
            'Oligo1_Sequence': str(oligos_all[oligo1_pos-1]),
            'Oligo2_Sequence': str(oligos_all[oligo2_pos-1]),
            'Score (average hits)': score
        })

    if not pairs_data:
        print("Warning: No valid probe pairs found after filtering!")
        return pd.DataFrame()

    probe_datasheet = pd.DataFrame(pairs_data)
    probe_datasheet = probe_datasheet.sort_values('Score (average hits)').reset_index(drop=True)

    print(f'Generated {len(probe_datasheet)} non-overlapping probe pairs')
    return probe_datasheet



def add_hairpin(probe_datasheet, hairpin):
    """
    Adds hairpin sequences to probes.

    Args:
        probe_datasheet (pandas.DataFrame): Probe information
        hairpin (str): Hairpin type ('B1'-'B5')

    Returns:
        pandas.DataFrame: Updated probe datasheet with hairpin sequences
    """
    print('Generating HCR probes...')

    hairpin_sequences = {
        'B1': ('gAggAgggCAgCAAACggAA', 'TAgAAgAgTCTTCCTTTACg'),
        'B2': ('CCTCgTAAATCCTCATCAAA', 'AAATCATCCAgTAAACCgCC'),
        'B3': ('gTCCCTgCCTCTATATCTTT', 'TTCCACTCAACTTTAACCCg'),
        'B4': ('CCTCAACCTACCTCCAACAA', 'ATTCTCACCATATTCgCTTC'),
        'B5': ('CTCACTCCCAATCTCTATAA', 'AACTACCCTACAAATCCAAT'),
    }

    if hairpin not in hairpin_sequences:
        raise ValueError(f"Invalid hairpin type '{hairpin}'. Must be one of: {', '.join(hairpin_sequences.keys())}")

    I_even, I_odd = hairpin_sequences[hairpin]

    probe_datasheet['HCRprobe1'] = [I_even.upper() + seq for seq in probe_datasheet['Oligo1_Sequence']]
    probe_datasheet['HCRprobe2'] = [seq + I_odd.upper() for seq in probe_datasheet['Oligo2_Sequence']]

    print('HCR probes designed')
    return probe_datasheet



def plot_sequence_alignment(sequence, sequence_id, probe_datasheet, top_n=None, figsize=(15, 8), save_path=None, output_fasta=None, return_records=True):
    """
    Creates a sequence-level visualization of probe alignments with probe pairs on same rows.
    Also optionally writes a FASTA alignment (padded to target length) of the target + probe sequences.

    Args:
        sequence (Bio.Seq): Input sequence
        sequence_id (str): ID of the input sequence
        probe_datasheet (pd.DataFrame): DataFrame containing probe information
        top_n (int, optional): Number of top-ranked pairs to plot. If None, plots all pairs
        figsize (tuple): Figure dimensions (width, height)
        save_path (str, optional): Path to save the plot image
        output_fasta (str, optional): Path to save FASTA alignment. If provided, writes padded FASTA.
        return_records (bool, optional): If True, return list of SeqRecord objects along with fig
    Returns:
        matplotlib.figure.Figure: The generated figure
        (optional) list[SeqRecord]: list of SeqRecord objects when return_records is True
    """
    # Setup aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1

    # Sort probe_datasheet by score and select top_n
    sorted_probes = probe_datasheet.sort_values('Score (average hits)', ascending=True)

    if top_n is not None:
        if not isinstance(top_n, int) or top_n <= 0:
            raise ValueError("top_n must be a positive integer")
        if top_n > len(sorted_probes):
            print(f"Warning: Requested {top_n} pairs but only {len(sorted_probes)} pairs available")
            top_n = len(sorted_probes)
        rows_to_plot = sorted_probes.head(top_n)
    else:
        rows_to_plot = sorted_probes

    # Setup figure with white background
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    fig, ax = plt.subplots(figsize=figsize)

    # Get tab20b colors for probes
    cmap = plt.cm.tab20b
    probe_colors = cmap(np.linspace(0, 1, len(rows_to_plot)))

    # Create custom colormap with white background
    custom_colors = np.vstack(([1, 1, 1, 1],  # white background
                             [0.8, 0.8, 0.8, 1],  # grey target
                             probe_colors))
    custom_cmap = mcolors.ListedColormap(custom_colors)

    # Process sequences and create visualization array
    target_seq = str(sequence)
    vis_array = np.zeros((len(rows_to_plot) + 1, len(sequence)))  # +1 for target sequence

    # Set target sequence row to 1 (will be grey)
    vis_array[0] = 1

    # Prepare records for FASTA output (start with target)
    records = [SeqRecord(Seq(target_seq), id=sequence_id, description="Target sequence")]

    # Process each probe pair
    for idx, (_, row) in enumerate(rows_to_plot.iterrows()):
        # Get probe sequences (use same orientation as plotted)
        oligo1_seq = str(Seq(row['Oligo1_Sequence']).reverse_complement())
        oligo2_seq = str(Seq(row['Oligo2_Sequence']).reverse_complement())

        # Get alignments
        alignment1 = aligner.align(target_seq, oligo1_seq)[0]
        alignment2 = aligner.align(target_seq, oligo2_seq)[0]

        # Get positions (start on target)
        pos1 = int(alignment1.coordinates[0][0])
        pos2 = int(alignment2.coordinates[0][0])

        # Fill probe positions with index+2 (since 0 is white, 1 is grey)
        for i in range(len(oligo1_seq)):
            if pos1 + i < len(sequence):
                vis_array[idx + 1][pos1 + i] = idx + 2
        for i in range(len(oligo2_seq)):
            if pos2 + i < len(sequence):
                vis_array[idx + 1][pos2 + i] = idx + 2

        # Prepare padded sequences for FASTA output
        pad1 = '-' * pos1 + oligo1_seq + '-' * (len(target_seq) - pos1 - len(oligo1_seq))
        pad2 = '-' * pos2 + oligo2_seq + '-' * (len(target_seq) - pos2 - len(oligo2_seq))

        # Add position information to descriptions
        oligo1_end = pos1 + len(oligo1_seq)
        oligo2_end = pos2 + len(oligo2_seq)
        
        records.append(SeqRecord(Seq(pad1), 
                               id=f"Pair_{idx+1}_Oligo1", 
                               description=f"Position: {pos1}-{oligo1_end}"))
        records.append(SeqRecord(Seq(pad2), 
                               id=f"Pair_{idx+1}_Oligo2", 
                               description=f"Position: {pos2}-{oligo2_end}"))

    # Create plot
    im = ax.imshow(vis_array, aspect='auto',
                   cmap=custom_cmap,
                   interpolation='none')

    # Set labels
    ax.set_yticks(range(len(rows_to_plot) + 1))
    labels = [sequence_id] + [f'Pair {i+1}' for i in range(len(rows_to_plot))]
    ax.set_yticklabels(labels)

    # Set x-axis ticks at sequence_length/20 intervals (avoid zero division)
    seq_len = max(1, len(sequence))
    tick_interval = max(1, seq_len // 20)
    tick_positions = list(range(0, seq_len, tick_interval))
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_positions)
    ax.set_xlabel('Sequence Position (nt)')

    # Remove grid
    ax.grid(False)

    # Adjust layout
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')

    if output_fasta:
        try:
            SeqIO.write(records, output_fasta, "fasta")
            print(f"Alignment FASTA saved to {output_fasta}")
        except Exception as e:
            print(f"Failed to write FASTA output: {e}", file=sys.stderr)

    if return_records:
        return fig, records

    return fig

