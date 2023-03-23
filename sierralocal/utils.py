def get_input_sequences(handle, return_dict=False):
    """
    Parse open file as FASTA, return a list of sequences.
    @param handle: _io.TextIOWrapper, Open stream to FASTA file in read mode
    @param return_dict: bool, Option to return a dictionary of sequences keyed
    by header
    @return: dict or list, if return_dict is True, returns dictionary; else a list
    of sequences
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
