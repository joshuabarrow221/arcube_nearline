#!/usr/bin/env python3
"""Compatibility entry point; the central lifetime module owns all plotting."""
import argparse
from lifetime import update_json

def main(input_file, output_file):
    update_json(input_file, [], output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--output_file', required=True)
    main(**vars(parser.parse_args()))
