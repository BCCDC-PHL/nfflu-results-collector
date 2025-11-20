#!/usr/bin/env python3
import os, sys
import argparse
import logging
import pandas as pd

from nfflu_results_collector.collector import Nfflu_Results_Collector


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--analysis-dir', required=True)
    parser.add_argument('-o', '--output-summary', required=True, help='Path to output summary CSV file')
    parser.add_argument('-O', '--output-mixture', help='Path to output mixture report CSV file')
    parser.add_argument('-s', '--output-symlinks', help='Path to output symlinks directory')
    parser.add_argument('-a', '--auto-nfflu', action='store_true', help='Run on auto-nfflu results structure')
    parser.add_argument('--log-level', default='INFO', help='Logging level (default: INFO)')

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    collector = Nfflu_Results_Collector({'auto-nfflu': args.auto_nfflu})
    collector.collect_run_summary(args.analysis_dir, args.output_summary)
    if args.output_mixture:
        collector.collect_mixture_report(args.analysis_dir, args.output_mixture)
    if args.output_symlinks:
        collector.symlink_consensus_fastas(args.analysis_dir, args.output_symlinks)

if __name__ == "__main__":
    sys.exit(main())