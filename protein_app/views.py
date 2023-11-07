import asyncio
import io
import time
import concurrent.futures
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd
from collections import Counter
from django.http import HttpResponse
from django.shortcuts import render
from .models import ProteinSequence
import logomaker
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor


matplotlib.use('Agg')  # Use the Agg backend for Matplotlib

def process_sequences(sequences_chunk):
    max_sequence_length = max(len(seq) for seq in sequences_chunk)
    columns = sorted(set(letter for seq in sequences_chunk for letter in seq if len(seq) >= max_sequence_length))

    df = pd.DataFrame(columns=columns)

    for col_idx in range(max_sequence_length):
        column_data = [seq[col_idx] if col_idx < len(seq) else None for seq in sequences_chunk]

        letter_counts = dict(Counter(column_data))
        letter_counts = {k: v for k, v in letter_counts.items() if k is not None}
        filtered_counts = {letter: count for letter, count in letter_counts.items() if count > 2}

        for letter, count in filtered_counts.items():
            df.at[col_idx, letter] = count

    df = df.fillna(0)
    df = df.divide(df.sum(axis=0), axis=1)
    return df

def process_sequences_and_generate_logo(sequences_group, center_letter, batch_size):
    results = []

    for i in range(0, len(sequences_group), batch_size):
        batch_sequences = sequences_group[i:i + batch_size]
        final_df = process_sequences(batch_sequences)
        final_df = final_df.apply(pd.to_numeric, errors='coerce').fillna(0)

        ww_logo = logomaker.Logo(final_df, color_scheme='NajafabadiEtAl2017', vpad=0.1, width=0.8)
        ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
        ww_logo.highlight_position(p=7, color='gold', alpha=0.5)

        for spine in ww_logo.ax.spines.values():
            spine.set_visible(False)

        ww_logo.ax.set_xlim([-0.5, 14.5])

        buffer = io.BytesIO()
        plt.savefig(buffer, format='svg', bbox_inches='tight', pad_inches=0)
        plt.close()

        svg_content = buffer.getvalue().decode('utf-8')
        buffer.close()

        results.append((svg_content, center_letter))

    return results

def analyze_protein(request):
    start_time = time.time()
    
    if request.method == 'POST':
        sequence = request.POST['protein_sequence']
        sequences = sequence.strip().split('\n')

        sequences_s, sequences_t, sequences_y = [], [], []

        for seq in sequences:
            center_letter = seq[7] if len(seq) > 7 else None
            if center_letter == 'S':
                sequences_s.append(seq)
            elif center_letter == 'T':
                sequences_t.append(seq)
            elif center_letter == 'Y':
                sequences_y.append(seq)

        chunk_size = 100
        batch_size = 1000  # Adjust the batch size based on available memory and performance
        graphs = []

        # Use ProcessPoolExecutor for parallel processing
        with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            tasks = []
            
            for sequences_group, center_letter in [(sequences_s, 'S'), (sequences_t, 'T'), (sequences_y, 'Y')]:
                task = executor.submit(process_sequences_and_generate_logo, sequences_group, center_letter, batch_size)
                tasks.append(task)
                
            for task in tasks:
                results = task.result()
                graphs.extend(results)

        # Save the protein sequence to the database
        protein_sequence = ProteinSequence(sequence=sequence)
        protein_sequence.save()

        context = {
            'graphs': graphs,
            'protein_sequence': sequence,
        }
        end_time = time.time()
        execution_time = end_time - start_time
        print("Execution time:", execution_time)

        return render(request, 'protein_app/results.html', context)

    return render(request, 'protein_app/input_form.html')