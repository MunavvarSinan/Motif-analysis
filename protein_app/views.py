from django.shortcuts import render
from django.http import HttpResponse
from .models import ProteinSequence
from collections import Counter
import pandas as pd
import numpy as np
import logomaker
import multiprocessing
from functools import partial
import io
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend for Matplotlib


def process_sequences(sequences_chunk):
    max_sequence_length = max(len(seq) for seq in sequences_chunk)
    columns = sorted(set(letter for seq in sequences_chunk for letter in seq if len(
        seq) >= max_sequence_length))

    df = pd.DataFrame(columns=columns)

    for col_idx in range(max_sequence_length):
        # Extract column data from valid sequences
        column_data = [seq[col_idx] if col_idx < len(
            seq) else None for seq in sequences_chunk]

        letter_counts = dict(Counter(column_data))

        # Filter out None values (sequences of different lengths)
        letter_counts = {k: v for k,
                         v in letter_counts.items() if k is not None}

        # Filter out letters with a count less than or equal to 2
        filtered_counts = {letter: count for letter,
                           count in letter_counts.items() if count > 2}

        # Update the DataFrame
        for letter, count in filtered_counts.items():
            df.at[col_idx, letter] = count

    # Fill NaN values with zeros
    df = df.fillna(0)

    # Normalize frequencies
    df = df.divide(df.sum(axis=0), axis=1)

    return df


def analyze_protein(request):
    if request.method == 'POST':
        sequence = request.POST['protein_sequence']

        # Split the sequence into lines if needed
        sequences = sequence.strip().split('\n')
        sequences_s = []
        sequences_t = []
        sequences_y = []

        # Split sequences into three groups based on center letter
        for seq in sequences:
            center_letter = seq[7] if len(seq) > 7 else None
            if center_letter == 'S':
                sequences_s.append(seq)
            elif center_letter == 'T':
                sequences_t.append(seq)
            elif center_letter == 'Y':
                sequences_y.append(seq)

        # Process sequences and generate graphs for each group
        graphs = []

        for sequences_group, center_letter in [(sequences_s, 'S'), (sequences_t, 'T'), (sequences_y, 'Y')]:
            # Example code:
            chunk_size = 100
            batch_size = 1000  # Number of sequences to process in each batch

            # Create a thread pool for parallel processing
            with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
                futures = []
                results = []

                # Batch sequences for processing
                for i in range(0, len(sequences_group), batch_size):
                    batch_sequences = sequences_group[i:i + batch_size]

                    # Submit the batch for processing
                    future = executor.submit(process_sequences, batch_sequences)
                    futures.append(future)

                # Collect results from futures
                for future in as_completed(futures):  # Use as_completed
                    result = future.result()
                    results.append(result)

            # Concatenate the list of DataFrames into the final DataFrame
            # Concatenate horizontally (along columns)
            final_df = pd.concat(results, axis=1)

            # Set the 'pos' index to match your provided DataFrame
            final_df.index.name = 'pos'

            # Normalize frequencies
            final_df = final_df.divide(final_df.sum(axis=0), axis=1)
            final_df = final_df.replace([np.inf, -np.inf, np.nan], 0)
            print(final_df)

            ww_logo = logomaker.Logo(final_df,
                                     color_scheme='NajafabadiEtAl2017',
                                     vpad=.1,
                                     width=.8)

            # Style using Logo methods
            # Change spacing and rotation
            ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
            ww_logo.highlight_position(p=7, color='gold', alpha=.5)

            # Remove the box around the graph
            for spine in ww_logo.ax.spines.values():
                spine.set_visible(False)

            # Set the x-axis limits to 1-15 without resizing the figure
            ww_logo.ax.set_xlim([-0.5, 14.5])  # Adjusted limits without resizing

            # Create an in-memory buffer for the Matplotlib figure
            buffer = io.BytesIO()

            # Save the Matplotlib figure as SVG in the buffer
            plt.savefig(buffer, format='svg', bbox_inches='tight',
                        pad_inches=0)  # Remove padding
            plt.close()

            # Get the SVG content from the buffer
            svg_content = buffer.getvalue().decode('utf-8')
            buffer.close()

            # Append the graph and center letter to the list
            graphs.append((svg_content, center_letter))

        # Save the protein sequence and analysis results to the database
        protein_sequence = ProteinSequence(sequence=sequence)
        protein_sequence.save()

        # Pass the graphs and protein sequence to the template
        context = {
            'graphs': graphs,
            'protein_sequence': sequence,
        }

        return render(request, 'protein_app/results.html', context)

    return render(request, 'protein_app/input_form.html')



# def analyze_protein(request):
#     if request.method == 'POST':
#         sequence = request.POST['protein_sequence']

#         # Split the sequence into lines if needed
#         sequences = sequence.strip().split('\n')
#         sequences_s = []
#         sequences_t = []
#         sequences_y = []

#         # Split sequences into three groups based on center letter
#         for seq in sequences:
#             center_letter = seq[7] if len(seq) > 7 else None
#             if center_letter == 'S':
#                 sequences_s.append(seq)
#             elif center_letter == 'T':
#                 sequences_t.append(seq)
#             elif center_letter == 'Y':
#                 sequences_y.append(seq)

#         # Print sequences for each group
#         print("Sequences with 'S' in the center:")
#         for seq in sequences_s:
#             print(seq)

#         print("Sequences with 'T' in the center:")
#         for seq in sequences_t:
#             print(seq)

#         print("Sequences with 'Y' in the center:")
#         for seq in sequences_y:
#             print(seq)

#         # Example code:
#         chunk_size = 100
#         batch_size = 1000  # Number of sequences to process in each batch

#         # Create a thread pool for parallel processing
#         with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
#             futures = []
#             results = []

#             # Batch sequences for processing
#             for i in range(0, len(sequences), batch_size):
#                 batch_sequences = sequences[i:i + batch_size]

#                 # Submit the batch for processing
#                 future = executor.submit(process_sequences, batch_sequences)
#                 futures.append(future)

#             # Collect results from futures
#             for future in as_completed(futures):  # Use as_completed
#                 result = future.result()
#                 results.append(result)
#         # Concatenate the list of DataFrames into the final DataFrame
#         # Concatenate horizontally (along columns)
#         final_df = pd.concat(results, axis=1)

#         # Set the 'pos' index to match your provided DataFrame
#         final_df.index.name = 'pos'

#         # Normalize frequencies
#         final_df = final_df.divide(final_df.sum(axis=0), axis=1)
#         final_df = final_df.replace([np.inf, -np.inf, np.nan], 0)
#         print(final_df)

#         ww_logo = logomaker.Logo(final_df,
#                                  color_scheme='NajafabadiEtAl2017',
#                                  vpad=.1,
#                                  width=.8)

#         # Style using Logo methods
#         # Change spacing and rotation
#         ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
#         ww_logo.highlight_position(p=7, color='gold', alpha=.5)

#         # Remove the box around the graph
#         for spine in ww_logo.ax.spines.values():
#             spine.set_visible(False)

#         # Set the x-axis limits to 1-15 without resizing the figure
#         ww_logo.ax.set_xlim([-0.5, 14.5])  # Adjusted limits without resizing

#         # Create an in-memory buffer for the Matplotlib figure
#         buffer = io.BytesIO()

#         # Save the Matplotlib figure as SVG in the buffer
#         plt.savefig(buffer, format='svg', bbox_inches='tight',
#                     pad_inches=0)  # Remove padding
#         plt.close()

#         # Get the SVG content from the buffer
#         svg_content = buffer.getvalue().decode('utf-8')
#         buffer.close()

#         # Save the protein sequence and analysis results to the database
#         protein_sequence = ProteinSequence(sequence=sequence)
#         protein_sequence.save()

#         # Pass the SVG content to the template
#         context = {
#             'svg_content': svg_content,
#             'protein_sequence': sequence,
#         }

#         return render(request, 'protein_app/results.html', context)

#     return render(request, 'protein_app/input_form.html')
