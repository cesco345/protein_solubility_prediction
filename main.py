import argparse
import sys
import traceback
import pandas as pd
import logging
from config import Config
from prediction import predict_solubility
from utils import setup_logging, validate_files, parse_fasta, generate_profiles

def main():
    parser = argparse.ArgumentParser(description="Predict protein solubility from sequence data.")
    parser.add_argument("input_file", help="Path to the input file (FASTA or composition CSV)")
    parser.add_argument("output_file", help="Path to the output predictions file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("-n", "--num_features", type=int, default=10, help="Number of features to use for prediction")
    parser.add_argument("--fasta", action="store_true", help="Input file is in FASTA format")
    args = parser.parse_args()

    config = Config()
    setup_logging(args.verbose)

    try:
        validate_files(args.input_file, config.SEQ_REFERENCE_DATA_FILE)
        
        if args.fasta:
            logging.info("Processing FASTA input")
            compositions = parse_fasta(args.input_file)
            logging.info(f"Compositions shape: {compositions.shape}")
            logging.info(f"Compositions columns: {compositions.columns}")
            logging.info(f"Compositions head:\n{compositions.head()}")
            
            profiles = generate_profiles(compositions)
            logging.info(f"Profiles shape: {profiles.shape}")
            logging.info(f"Profiles columns: {profiles.columns}")
            logging.info(f"Profiles head:\n{profiles.head()}")
        else:
            logging.info("Processing CSV input")
            compositions = pd.read_csv(args.input_file)
            profiles = pd.read_csv(args.input_file.replace('composition', 'profiles'))

        predict_solubility(compositions, profiles, args.output_file, config, args.num_features)
        logging.info(f"Predictions written to {args.output_file}")
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        logging.error("Traceback:")
        logging.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
