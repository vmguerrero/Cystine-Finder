import argparse
import requests
import xml.etree.ElementTree as ET
from Bio import SeqIO
from io import StringIO
import csv
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import sys

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

icon_path = resource_path('Icon.ico')

def read_peptides(file_path):
    peptides = []
    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            peptides.append({"Protein": row["Protein"], "Peptide": row["Peptide"]})
    return peptides

def get_protein_data(Protein):
    xml_url = f"https://rest.uniprot.org/uniprotkb/{Protein}.xml"
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{Protein}.fasta"

    xml_response = requests.get(xml_url)
    xml_content = xml_response.content

    try:
        root = ET.fromstring(xml_content)
    except ET.ParseError as e:
        raise Exception(f"Error parsing {Protein} xml") from e

    fasta_response = requests.get(fasta_url)
    fasta_io = StringIO(fasta_response.text)
    record = SeqIO.read(fasta_io, "fasta")
    sequence = str(record.seq)

    return root, sequence

def parse_disulfide_bonds(root):
    namespace = {'uniprot': 'http://uniprot.org/uniprot'}
    disulfide_bonds = []
    cysteine_positions = {}
    for feature in root.findall('.//uniprot:feature[@type="disulfide bond"]', namespaces=namespace):
        begin = feature.find('.//uniprot:begin', namespaces=namespace)
        end = feature.find('.//uniprot:end', namespaces=namespace)
        position = feature.find('.//uniprot:position', namespaces=namespace)
        if begin is not None and end is not None:
            try:
                begin_pos = int(begin.get('position'))
                end_pos = int(end.get('position'))
                disulfide_bonds.append((begin_pos, end_pos))
                cysteine_positions[begin_pos] = end_pos
                cysteine_positions[end_pos] = begin_pos
            except AttributeError as e:
                print(f"Error in feature: {ET.tostring(feature, encoding='unicode')}")
                raise e
        elif position is not None:
            try:
                pos = int(position.get('position'))
                disulfide_bonds.append((pos, 0))
                cysteine_positions[pos] = 0
            except AttributeError as e:
                print(f"Error in feature: {ET.tostring(feature, encoding='unicode')}")
                raise e
        else:
            print(f"Element 'begin', 'end' or 'position' not found for feature: {ET.tostring(feature, encoding='unicode')}")
    return disulfide_bonds, cysteine_positions

def search_disulfide_bonds(protein_sequence, peptides, cysteine_positions):
    results = []
    for peptide_info in peptides:
        peptide = peptide_info["Peptide"]
        peptide_positions = []
        start_pos = protein_sequence.find(peptide)
        if start_pos == -1:
            results.append((peptide, [], [], [], peptide_info["Protein"]))
            continue
        
        peptide_cysteine_positions = []
        disulfide_partner_positions = []
        for i, aa in enumerate(peptide):
            if aa == 'C':
                peptide_abs_pos = start_pos + i + 1
                peptide_cysteine_positions.append(i + 1)
                if peptide_abs_pos in cysteine_positions:
                    peptide_positions.append(peptide_abs_pos)
                    disulfide_partner_positions.append(cysteine_positions[peptide_abs_pos])
                else:
                    disulfide_partner_positions.append(None)
        
        results.append((peptide, peptide_cysteine_positions, peptide_positions, disulfide_partner_positions, peptide_info["Protein"]))
    return results

def write_results(output_file, results):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(["Peptide", "CysteinePositionsInPeptide", "CysteinePositionsInProtein", "DisulfidePartnerPositionsInProtein", "Protein"])
        for result in results:
            writer.writerow(result)

def main(peptide_file, output_file):
    peptides = read_peptides(peptide_file)
    all_results = []

    proteins = set([peptide["Protein"] for peptide in peptides])
    for Protein in proteins:
        root, sequence = get_protein_data(Protein)
        disulfide_bonds, cysteine_positions = parse_disulfide_bonds(root)
        protein_peptides = [peptide for peptide in peptides if peptide["Protein"] == Protein]
        peptide_disulfide_info = search_disulfide_bonds(sequence, protein_peptides, cysteine_positions)
        
        for result in peptide_disulfide_info:
            all_results.append(result)

    write_results(output_file, all_results)

    update_results_tree(all_results)

def update_results_tree(results):
    for i in treeview.get_children():
        treeview.delete(i)
    for result in results:
        treeview.insert("", "end", values=result)

def browse_file():
    file_path = filedialog.askopenfilename()
    if file_path:
        entry_input_file.delete(0, tk.END)
        entry_input_file.insert(0, file_path)

def run_process():
    input_file = entry_input_file.get()
    if input_file:
        output_file = os.path.splitext(input_file)[0] + "_cystines.tsv"
        try:
            main(input_file, output_file)
            messagebox.showinfo("Success", f"Output file: {output_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))
    else:
        messagebox.showwarning("Error", "Please, provide an input file path")

app = tk.Tk()
app.title("Cystine Finder. Software for Locating Disulfide Bonds in Protein Sequences. Developed by Victor M. Guerrero-Sanchez")

frame = ttk.Frame(app, padding="10")
frame.grid(row=0, column=0)

label_input_file = ttk.Label(frame, text="Input file:")
label_input_file.grid(row=0, column=0, pady=5, padx=5)
entry_input_file = ttk.Entry(frame, width=50)
entry_input_file.grid(row=0, column=1, pady=5, padx=5)
button_browse_input = ttk.Button(frame, text="Browse", command=browse_file)
button_browse_input.grid(row=0, column=2, pady=5, padx=5)

button_run = ttk.Button(app, text="Run", command=run_process)
button_run.grid(row=1, column=0, pady=20)

columns = ["Peptide", "CysteinePositionsInPeptide", "CysteinePositionsInProtein", "DisulfidePartnerPositionsInProtein", "Protein"]
treeview = ttk.Treeview(app, columns=columns, show="headings")
for col in columns:
    treeview.heading(col, text=col)
treeview.grid(row=2, column=0, padx=10, pady=10, sticky=(tk.W, tk.E, tk.N, tk.S))
app.iconbitmap(icon_path)

app.mainloop()
