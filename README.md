# Protein Solubility Prediction

This project predicts protein solubility based on sequence composition and profile data. A Python-based protein solubility prediction tool, building on the Protein-sol algorithm developed by the Warwicker group at the University of Manchester.

## Installation

1. Clone this repository:

git clone https://github.com/cesco345/protein_solubility_prediction.git
cd protein_solubility_prediction

2. Install the required packages:
pip install -r requirements.txt

## Usage

Run the prediction script with:
python main.py <input_file> <output_file> [options]

### Arguments:

- `input_file`: Path to the input file (FASTA or composition CSV)
- `output_file`: Path to the output predictions file

### Options:

- `-v, --verbose`: Enable verbose output
- `-n NUM_FEATURES, --num_features NUM_FEATURES`: Number of features to use for prediction (default: 10)
- `--fasta`: Specify if the input file is in FASTA format

### Examples:

For FASTA input:
python main.py data/example_input.fasta predictions.txt --fasta -v -n 8

For CSV input:
python main.py data/example_composition.csv predictions.txt -v -n 8

## Input File Formats

1. FASTA file: Standard FASTA format with protein sequences
2. Composition file: CSV file with columns including 'ID', 'K-R', 'D-E', 'naa', 'KyteDoo', 'abs-charge', etc.

## Output Format

The output file is a CSV with the following columns:
- ID: Sequence identifier
- percent-sol: Predicted solubility percentage
- scaled-sol: Scaled solubility (0-1)
- population-sol: Population average solubility
- pI: Isoelectric point

## Configuration

Adjust parameters in `src/config.py` to modify file paths, column names, and other settings.

## Logging

Logs are written to `solubility_prediction.log`. Use the `-v` option for verbose logging.

## License

This project is licensed under the MIT License.



