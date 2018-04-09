
#!/usr/bin/env python
"""
Given a TBA coefficients file or a TBA significance file, maps the 
motif names to gene names
"""

### imports ###
import argparse
import numpy as np
import os
import pandas as pd
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a TBA coefficients file or \
    a TBA significance file, maps the motif names to gene names')
    parser.add_argument("feature_path",
        help="path to a standardized feature set created by create_features.py",
        type = str)
    
    parser.add_argument
    parser.add_argument("motif_files",
        help="list of motif files",
        type=str,
        nargs="+")
    
    
    # parse arguments
    args = parser.parse_args()

    feature_path = args.feature_path
    label_path = args.label_path
    output_path = args.output_path
