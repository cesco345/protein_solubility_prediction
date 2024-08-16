import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

def setup_logging(verbose):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(filename="solubility_prediction.log", level=level, format="%(asctime)s - %(levelname)s - %(message)s")
    console = logging.StreamHandler()
    console.setLevel(level)
    logging.getLogger("").addHandler(console)

def validate_files(*file_paths):
    for file_path in file_paths:
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        if not path.is_file():
            raise IsADirectoryError(f"Expected a file, got a directory: {file_path}")
        logging.info(f"Validated file: {file_path}")

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import logging

def parse_fasta(fasta_file):
    compositions = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq)
            analysis = ProteinAnalysis(seq)
            aa_count = analysis.count_amino_acids()
            
            # Define the disorder propensity scale
            disorder_scale = {
                'A': 0.06, 'R': -0.18, 'N': -0.01, 'D': -0.05, 'C': 0.02, 'Q': -0.01, 'E': -0.01, 
                'G': 0.15, 'H': -0.06, 'I': 0.08, 'L': 0.09, 'K': -0.1, 'M': 0.08, 'F': 0.11, 
                'P': 0.04, 'S': 0.05, 'T': 0.04, 'W': 0.1, 'Y': 0.05, 'V': 0.07
            }
            
            composition = {
                'ID': record.id,
                'K-R': aa_count.get('K', 0) - aa_count.get('R', 0),
                'D-E': aa_count.get('D', 0) - aa_count.get('E', 0),
                'naa': len(seq),
                'KyteDoo': analysis.gravy(),
                'abs-charge': abs(analysis.charge_at_pH(7.0)),
                'pI': analysis.isoelectric_point(),
                'FoldIndex': 2.785 * analysis.gravy() - (abs(analysis.charge_at_pH(7.0)) / len(seq)) - 1.151,
                'disorder': sum(analysis.protein_scale(window=7, param_dict=disorder_scale))/len(seq),
                'entropy': analysis.molar_extinction_coefficient()[0]
            }
            compositions.append(composition)
        logging.info(f"Parsed {len(compositions)} sequences from FASTA file")
        return pd.DataFrame(compositions)
    except Exception as e:
        logging.error(f"Error in parse_fasta: {str(e)}")
        raise


def generate_profiles(compositions):
    profiles = []
    for _, row in compositions.iterrows():
        for prop in ['K-R', 'D-E', 'naa', 'KyteDoo', 'abs-charge', 'pI', 'FoldIndex', 'disorder', 'entropy']:
            profiles.append({
                'ID': row['ID'],
                'Property': prop,
                'Value': row[prop]
            })
    return pd.DataFrame(profiles)
