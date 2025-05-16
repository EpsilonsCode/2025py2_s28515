#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever with filtering, CSV reporting and visualization
"""
import time

from Bio import Entrez
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys


class NCBIRetriever:
    def __init__(self, email, api_key=None):
        """Initialize with NCBI credentials."""
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid, min_len=None, max_len=None):
        """Search for records with optional length filtering."""
        search_term = f"txid{taxid}[Organism]"
        if min_len:
            search_term += f" AND {min_len}:{max_len or ''}[SLEN]"

        handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
        results = Entrez.read(handle)
        return results if int(results["Count"]) > 0 else None

    def fetch_records(self, webenv, query_key, count, batch_size=10):
        """Fetch records in batches."""
        records = []
        for start in range(0, count, batch_size):
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key
            )
            records.extend(handle.read().split('//\n')[:-1])
            time.sleep(0.34)  # NCBI rate limit
        return records

    def parse_records(self, records):
        """Extract relevant data from GenBank records."""
        data = []
        for record in records:
            lines = record.split('\n')
            accession = length = description = None
            for line in lines:
                if line.startswith('ACCESSION'):
                    accession = line.split()[1]
                elif line.startswith('DEFINITION'):
                    description = line[12:].replace('\n', ' ').strip()
                elif line.startswith('LOCUS'):
                    length = int(line.split()[2])
            if accession:
                data.append([accession, length, description])
        return data

    def generate_csv(self, data, filename):
        """Generate CSV report from parsed data."""
        df = pd.DataFrame(data, columns=['Accession', 'Length', 'Description'])
        df.to_csv(filename, index=False)
        return df

    def plot_lengths(self, df, filename):
        """Generate length distribution plot."""
        df = df.sort_values('Length', ascending=False)
        plt.figure(figsize=(12, 6))
        plt.plot(df['Accession'], df['Length'], 'o-')
        plt.xticks(rotation=90)
        plt.xlabel('Accession Number')
        plt.ylabel('Sequence Length')
        plt.title('Sequence Length Distribution')
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()


def main():
    print("NCBI GenBank Data Retriever")
    print("---------------------------")

    # Get input from user
    email = input("Enter your NCBI registered email: ")
    api_key = input("Enter your NCBI API key (optional, press Enter to skip): ")
    taxid = input("Enter Taxonomic ID: ")
    min_len = input("Enter minimum sequence length (optional, press Enter to skip): ")
    max_len = input("Enter maximum sequence length (optional, press Enter to skip): ")
    output = input("Enter output base name (default: 'output'): ") or 'output'

    # Convert length inputs to integers or None
    min_len = int(min_len) if min_len else None
    max_len = int(max_len) if max_len else None

    retriever = NCBIRetriever(email, api_key)

    print(f"Searching for taxID: {taxid}...")
    results = retriever.search_taxid(taxid, min_len, max_len)

    if not results:
        sys.exit("No records found")

    count = int(results["Count"])
    print(f"Found {count} records")

    records = retriever.fetch_records(results["WebEnv"], results["QueryKey"], min(count, 100))
    data = retriever.parse_records(records)

    csv_file = f"{output}.csv"
    df = retriever.generate_csv(data, csv_file)
    print(f"Saved report to {csv_file}")

    plot_file = f"{output}.png"
    retriever.plot_lengths(df, plot_file)
    print(f"Saved plot to {plot_file}")


if __name__ == "__main__":
    main()
    