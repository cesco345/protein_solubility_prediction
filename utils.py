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

def parse_fasta(fasta_file):
    compositions = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        analysis = ProteinAnalysis(seq)
        composition = {
            'ID': record.id,
            'K-R': analysis.count_amino_acids()['K'] - analysis.count_amino_acids()['R'],
            'D-E': analysis.count_amino_acids()['D'] - analysis.count_amino_acids()['E'],
            'naa': len(seq),
            'KyteDoo': analysis.gravy(),
            'abs-charge': abs(analysis.charge_at_pH(7.0)),
            'pI': analysis.isoelectric_point(),
            'FoldIndex': 2.785 * analysis.gravy() - (abs(analysis.charge_at_pH(7.0)) / len(seq)) - 1.151,
            'disorder': sum(analysis.protein_scale(window=7, param_dict='ED'))/len(seq),
            'entropy': analysis.molar_extinction_coefficient()[0]
        }
        compositions.append(composition)
    return pd.DataFrame(compositions)

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
