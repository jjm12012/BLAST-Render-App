import streamlit as st
import os
import subprocess
import pandas as pd
from Bio import SeqIO

# Set up the Streamlit UI
st.title("BLAST Analysis Web App")
st.write("Upload .ab1 files and a reference sequence (.txt) to process BLAST alignments.")

# Create necessary directories
base_dir = "blast_data"
fasta_dir = os.path.join(base_dir, "FASTA_Files")
blast_output_dir = os.path.join(base_dir, "BLAST_Results")
blast_db_dir = os.path.join(base_dir, "BLAST_DB")

for folder in [fasta_dir, blast_output_dir, blast_db_dir]:
    os.makedirs(folder, exist_ok=True)

# File uploader
uploaded_ab1_files = st.file_uploader("Upload .ab1 files", type=["ab1"], accept_multiple_files=True)
uploaded_ref_file = st.file_uploader("Upload reference sequence (.txt)", type=["txt"])

if st.button("Process Files"):
    if not uploaded_ab1_files or not uploaded_ref_file:
        st.error("Please upload both .ab1 files and a reference sequence.")
    else:
        # Save uploaded files
        for uploaded_file in uploaded_ab1_files:
            file_path = os.path.join(base_dir, uploaded_file.name)
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

        ref_file_path = os.path.join(base_dir, uploaded_ref_file.name)
        with open(ref_file_path, "wb") as f:
            f.write(uploaded_ref_file.getbuffer())

        # Convert reference sequence to FASTA
        reference_fasta = os.path.join(base_dir, "reference.fasta")
        with open(ref_file_path, "r") as ref_txt, open(reference_fasta, "w") as ref_fasta:
            lines = ref_txt.readlines()
            ref_fasta.write(">Reference_Sequence\n")
            ref_fasta.writelines(lines)

        # Convert .ab1 to FASTA
        for filename in os.listdir(base_dir):
            if filename.lower().endswith(".ab1"):
                input_path = os.path.join(base_dir, filename)
                output_path = os.path.join(fasta_dir, filename.replace(".ab1", ".fasta"))
                try:
                    record = SeqIO.read(input_path, "abi")
                    with open(output_path, "w") as fasta_file:
                        SeqIO.write(record, fasta_file, "fasta")
                except Exception as e:
                    st.error(f"Error converting {filename}: {e}")

        # Create BLAST database
        blast_db_path = os.path.join(blast_db_dir, "reference_db")
        makeblastdb_cmd = f'makeblastdb -in "{reference_fasta}" -dbtype nucl -out "{blast_db_path}"'
        try:
            subprocess.run(makeblastdb_cmd, shell=True, check=True)
            st.success("BLAST database created successfully.")
        except subprocess.CalledProcessError as e:
            st.error(f"Error creating BLAST database: {e}")

        # Run BLASTN
        results = []
        for filename in os.listdir(fasta_dir):
            if filename.endswith(".fasta"):
                query_fasta = os.path.join(fasta_dir, filename)
                output_file = os.path.join(blast_output_dir, filename.replace(".fasta", "_blast_results.txt"))

                blastn_cmd = f'blastn -query "{query_fasta}" -db "{blast_db_path}" -out "{output_file}" -outfmt 7'
                try:
                    subprocess.run(blastn_cmd, shell=True, check=True)
                    st.success(f"BLASTN completed for {filename}. Results saved.")
                except subprocess.CalledProcessError as e:
                    st.error(f"Error running BLASTN for {filename}: {e}")

                # Process BLAST results for summary
                with open(output_file, "r") as blast_file:
                    lines = blast_file.readlines()
                cleaned_results = [line for line in lines if not line.startswith("#") and line.strip()]
                
                if cleaned_results:
                    results.append([filename] + cleaned_results[0].strip().split("\t")[:6])

        # Display results as a table
        if results:
            df = pd.DataFrame(results, columns=["Query File", "Query ID", "Subject ID", "% Identity", "Alignment Length", "Mismatches", "Gap Opens"])
            st.write("### BLAST Summary Results")
            st.dataframe(df)

