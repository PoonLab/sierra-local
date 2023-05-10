import tempfile
import os

def get_input_sequences(handle, return_dict=False):
    """
    Parse open file as FASTA, return a list of sequences.

    :param handle:  Open stream to FASTA file in read mode
    :param return_dict:  Option to return a dictionary of sequences keyed by header
    :return:  if return_dict is True, returns dictionary; else a list of sequences
    """
    headers = []
    sequences = []

    header = ''
    sequence = ''
    for line in handle:
        if line.startswith('$'):
            # skip comment
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                headers.append(header)
                sequences.append(sequence)
                sequence = ''  # reset
            header = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()

    # handle last entry
    headers.append(header)
    sequences.append(sequence)

    return dict(zip(headers, sequences)) if return_dict else sequences



