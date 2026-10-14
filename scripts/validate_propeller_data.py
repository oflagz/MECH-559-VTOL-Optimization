"""Usage: python -m scripts.validate_propeller_data POLARS.csv TESTS.csv"""
import argparse
from src.rotor.polar_table import PolarTable
from src.rotor.validation import validate_csv, summary

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("polars",help="CSV with reynolds,alpha_deg,cl,cd")
    parser.add_argument("measurements",help="CSV with propeller measurements and geometry")
    args = parser.parse_args()
    results = validate_csv(args.measurements,PolarTable(args.polars))
    for row in results:
        print(row)
    print("Summary:",summary(results))

if __name__ == "__main__":
    main()
