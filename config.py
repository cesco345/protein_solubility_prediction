import os

class Config:
    # File paths
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(BASE_DIR, "data")
    SEQ_REFERENCE_DATA_FILE = os.path.join(DATA_DIR, "seq_reference_data.txt")

    # Column names
    ID_COLUMN = "ID"
    ORF_ID_COLUMN = " ORF-ID"
    PROPERTY_COLUMN = "Property"
    VALUE_COLUMN = "Value"

    # Prediction parameters
    LOW_SOL_KEY = "LOW"
    TOP_SOL_KEY = "TOP"
    POP_SOL_KEY = "POP"
    ZDF_KEY = "ZDF"
    AVG_KEY = "AVG"

    # Output format
    OUTPUT_HEADER = "SEQUENCE PREDICTIONS,ID,percent-sol,scaled-sol,population-sol,pI"
    OUTPUT_FORMAT = "SEQUENCE PREDICTIONS,{},{:.3f},{:.3f},{:.3f},{:.2f}"

    # Logging
    LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    LOG_FILE = "solubility_prediction.log"
