# DNA-translation-and-logistics-curve-
DNA translation and logistics curve 
def translate_dna_to_protein(dna_seq):
    codon_table = {
        "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
        "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
        "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
        "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",                 
        "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
        "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
        "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
        "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
        "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
        "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
        "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
        "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
        "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
        "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
        "TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
        "TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
    }
    
    protein = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein += codon_table.get(codon.upper(), "X")  # "X" for unknown codons
    return protein

import numpy as np
import pandas as pd

def logistic_growth(time, K=1.0, r=1.0, N0=0.01, lag_length=5, exp_length=20):
    growth = np.zeros_like(time)
    for i, t in enumerate(time):
        if t < lag_length:
            growth[i] = N0
        elif lag_length <= t < lag_length + exp_length:
            dt = t - lag_length
            growth[i] = K / (1 + ((K - N0) / N0) * np.exp(-r * dt))
        else:
            growth[i] = K
    return growth

def generate_growth_df(n_curves=100, time_points=100):
    data = []
    time = np.linspace(0, 50, time_points)

    for i in range(n_curves):
        K = np.random.uniform(0.8, 1.2)
        r = np.random.uniform(0.2, 1.0)
        lag = np.random.uniform(2, 10)
        exp = np.random.uniform(10, 30)
        N0 = np.random.uniform(0.001, 0.05)
        growth = logistic_growth(time, K=K, r=r, N0=N0, lag_length=lag, exp_length=exp)
        data.append(growth)

    df = pd.DataFrame(data, columns=[f"T{int(t)}" for t in time])
    return df

def time_to_80_percent_max(time, growth_curve):
    threshold = 0.8 * max(growth_curve)
    for t, pop in zip(time, growth_curve):
        if pop >= threshold:
            return t
    return None  # in case threshold never reached

def hamming_distance(str1, str2):
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len, "_")
    str2 = str2.ljust(max_len, "_")
    
    return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))

# Example usage:
slack_username = "Amarachi"
twitter_handle = "Laura_Olarinde"
distance = hamming_distance(slack_username, twitter_handle)

# Generate growth curves
df = generate_growth_df()

# Get one example curve to calculate time to 80%
time = np.linspace(0, 50, 100)
sample_curve = df.iloc[0].values
time_80 = time_to_80_percent_max(time, sample_curve)

print("Time to 80% max:", time_80)
print("Hamming Distance:", distance)
