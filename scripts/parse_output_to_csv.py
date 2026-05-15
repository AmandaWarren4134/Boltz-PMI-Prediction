import os
import argparse
import json
import csv
import numpy as np
from pathlib import Path

# PROJECT_ROOT = Path(os.environ["PROJECT_ROOT"])
# SCRIPT_DIR = Path(os.environ["SCRIPT_DIR"])
# YAML_DIR = Path(os.environ["YAML_DIR"])
OUTPUT_DIR = Path(os.environ["OUTPUT_DIR"])
PROCESSED_DIR = Path(os.environ["PROCESSED_DIR"])

def get_predictions_from_output(output_folder: Path, verbose: bool = False):
    confidence_list = []
    num_missing_confidence_files = 0
    num_missing_affinity_files = 0

    # Each subfolder corresponds to the prediction
    # associated with a single YAML file.
    # The name follows the pattern:
    # boltz_results_{protein_name}__{ligand_name}
    for subfolder in output_folder.iterdir():
        if not subfolder.is_dir():
            continue

        stem = subfolder.name.split("boltz_results_")[1]

        if "__" in stem:
            protein_name = stem.split("__")[0] 
            ligand_name = stem.split("__")[1]
        else:
            protein_name, ligand_name = stem, "UNKNOWN"

        pred_folder = subfolder / "predictions" / stem
        confidence_file = pred_folder / f"confidence_{stem}_model_0.json"
        affinity_file = pred_folder / f"affinity_{stem}.json"

        if not confidence_file.exists():
            if verbose:
                print(f"\tMissing confidence: {stem}")
            num_missing_confidence_files += 1

        if not affinity_file.exists():
            if verbose:
                print(f"\tMissing affinity: {stem}")
            num_missing_affinity_files += 1

        if not confidence_file.exists() or not affinity_file.exists():
            continue

        with open(confidence_file, "r") as f:
            confidence = json.load(f)

        with open(affinity_file, "r") as f:
            affinity = json.load(f)

        protein_ligand_data = {'protein': protein_name, 'ligand': ligand_name}
        prediction_result = {**protein_ligand_data, **confidence, **affinity}
        confidence_list.append(prediction_result)

    if verbose:
        print(f"\tMissing confidence files: {num_missing_confidence_files}")
        print(f"\tMissing affinity files: {num_missing_affinity_files}")

    add_additional_affinity_units(confidence_list)

    return confidence_list


def add_additional_affinity_units(confidence_list: list):
    # Add micromolar affinity, molar affinity, and -log10(molar affinity)

    for item in confidence_list:
        affinity_pred_value = item.get("affinity_pred_value")

        micromolar_affinity_pred_value = 10 ** affinity_pred_value
        molar_affinity_pred_value = micromolar_affinity_pred_value / 1e6
        neg_log_molar_affinity_pred_value = -np.log10(molar_affinity_pred_value)

        item["micromolar_affinity_pred_value"] = micromolar_affinity_pred_value
        item["molar_affinity_pred_value"] = molar_affinity_pred_value
        item["neg_log_molar_affinity_pred_value"] = neg_log_molar_affinity_pred_value


def save_data_as_csv(job_name: str, data: list):
    if not PROCESSED_DIR.exists():
        PROCESSED_DIR.mkdir(parents=True)

    csv_file = PROCESSED_DIR / f"confidence_predictions_{job_name}.csv"

    if not data:
        print(f"No data for {job_name}. Skipping CSV creation.")
        return

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

    print(f"Saved confidence predictions to {csv_file.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output_folder", type=str, nargs="+", help="Name(s) of the output folder(s) to process. If not provided, --process_all will be used")
    ap.add_argument("--job_name", type=str, help="Name for the output CSV")
    ap.add_argument("--process_all", action="store_true", help="Process all output folders in the outputs directory")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    
    folders = []

    if args.process_all:
        print("Processing all output folders...")
        for folder_path in OUTPUT_DIR.iterdir():
            if folder_path.is_dir():
                folders.append(folder_path)
    elif args.output_folder:
        for folder_name in args.output_folder:
            folder_path = OUTPUT_DIR / folder_name
            if folder_path.is_dir():
                folders.append(folder_path)
            else:
                print(f"Output folder not found: {folder_name}")

    for folder_path in folders:
        print(f"Processing output folder: {folder_path.name}")
        confidence_list = get_predictions_from_output(folder_path, verbose=args.verbose)

        if not confidence_list:
            print(f"No confidence data found for output folder: {folder_path.name}")
            continue

        if args.job_name and len(folders) == 1:
            job_name = args.job_name
        else:            
            job_name = folder_path.name
        
        save_data_as_csv(job_name, confidence_list)


if __name__ == "__main__":
    main()
