#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Generate random quadrature data.")
    parser.add_argument("-q", "--quadrature", type=int, required=True,
                        help="Number of quadrature points.")
    parser.add_argument("-o", "--output_file", type=str, required=True,
                        help="Output CSV file name.")

    args = parser.parse_args()

    Q = args.quadrature
    filename = args.output_file

    np.random.seed(42)

    data = {
        'dXdx_init': np.random.random(Q * 9),
        'dudX': np.random.random(Q * 9),
        'ddudX': np.random.random(Q * 9),
    }

    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)

    print(f"Data saved to '{filename}' with {Q} quadrature points.")

if __name__ == "__main__":
    main()
