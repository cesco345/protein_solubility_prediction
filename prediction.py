import logging
import numpy as np
import pandas as pd
from typing import Dict, List



def read_reference_data(file_path: str) -> Dict:
    logging.info(f"Reading reference data from {file_path}")
    ref_data = {}
    current_section = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split(',')
            if parts[0] in ['POP', 'TOP', 'LOW', 'ZDF']:
                current_section = parts[0]
                if current_section not in ref_data:
                    ref_data[current_section] = []
            if current_section:
                ref_data[current_section].append(parts)

    for section in ['POP', 'TOP', 'LOW']:
        if section in ref_data:
            df = pd.DataFrame(ref_data[section])
            ref_data[section] = {
                'AVG': float(ref_data[section][0][3]),
                'data': df
            }

    if 'ZDF' in ref_data:
        headers = ref_data['ZDF'][0][4:]
        zscore_diff = [float(x) for x in ref_data['ZDF'][1][4:]]
        use_for_prob = ref_data['ZDF'][2][4:]
        ref_data['ZDF'] = {
            'headers': headers,
            'zscore_diff': zscore_diff,
            'use_for_prob': use_for_prob
        }

    logging.info("Reference data loaded successfully")
    return ref_data

def calculate_prediction(row: pd.Series, headers: List[str], weights: np.ndarray, low_values: np.ndarray, top_values: np.ndarray) -> float:
    feature_values = row[headers].astype(float).values
    scaled_values = np.clip((feature_values - low_values) / (top_values - low_values), 0, 1)
    prediction = np.sum(weights * scaled_values) / np.sum(np.abs(weights))
    return prediction

def predict_solubility(compositions: pd.DataFrame, profiles: pd.DataFrame, output_file: str, config: 'Config', num_features: int):
    logging.info(f"Starting solubility prediction with {num_features} features")
    ref_data = read_reference_data(config.SEQ_REFERENCE_DATA_FILE)
    
    logging.info(f"Composition shape: {compositions.shape}")
    logging.info(f"Profiles shape: {profiles.shape}")

    compositions = compositions.rename(columns={config.ORF_ID_COLUMN: config.ID_COLUMN})
    profiles_pivot = profiles.pivot(index=config.ID_COLUMN, columns=config.PROPERTY_COLUMN, values=config.VALUE_COLUMN)
    profiles_pivot.reset_index(inplace=True)

    data = pd.merge(compositions, profiles_pivot, on=config.ID_COLUMN, how='inner')

    features_to_use = [i for i, use in enumerate(ref_data[config.ZDF_KEY]['use_for_prob']) if use == 'y'][:num_features]
    weights = np.array([ref_data[config.ZDF_KEY]['zscore_diff'][i] for i in features_to_use])
    headers = [ref_data[config.ZDF_KEY]['headers'][i] for i in features_to_use]

    low_values = ref_data[config.LOW_SOL_KEY]['data'].iloc[2, 4:].astype(float).values[features_to_use]
    top_values = ref_data[config.TOP_SOL_KEY]['data'].iloc[2, 4:].astype(float).values[features_to_use]

    logging.info(f"Features used for prediction: {headers}")
    logging.info(f"Feature weights: {weights}")

    data['prediction'] = data.apply(lambda row: calculate_prediction(row, headers, weights, low_values, top_values), axis=1)

    low_sol = ref_data[config.LOW_SOL_KEY][config.AVG_KEY]
    top_sol = ref_data[config.TOP_SOL_KEY][config.AVG_KEY]
    pop_sol = ref_data[config.POP_SOL_KEY][config.AVG_KEY]

    with open(output_file, 'w') as out_handle:
        out_handle.write(f"{config.OUTPUT_HEADER}\n")
        for _, row in data.iterrows():
            adjusted_prediction = (row['prediction'] + 1) / 2
            percent_sol = low_sol + (top_sol - low_sol) * adjusted_prediction
            scaled_sol = (percent_sol - low_sol) / (top_sol - low_sol)
            out_handle.write(f"{config.OUTPUT_FORMAT.format(row[config.ID_COLUMN], percent_sol, scaled_sol, pop_sol, row['pI'])}\n")

    logging.info(f"Prediction Statistics:")
    logging.info(f"Number of predictions: {len(data)}")
    logging.info(f"Mean raw prediction: {data['prediction'].mean():.3f}")
    logging.info(f"Min raw prediction: {data['prediction'].min():.3f}")
    logging.info(f"Max raw prediction: {data['prediction'].max():.3f}")
    
    adjusted_predictions = (data['prediction'] + 1) / 2
    logging.info(f"Mean adjusted prediction: {adjusted_predictions.mean():.3f}")
    logging.info(f"Min adjusted prediction: {adjusted_predictions.min():.3f}")
    logging.info(f"Max adjusted prediction: {adjusted_predictions.max():.3f}")
