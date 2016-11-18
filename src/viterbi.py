from data import parse_fna
from math import log, exp
import sys
import json

def run_viterbi(gene, params):


    a = log(params['transitions']['begin'][0]) + log(params['emissions']['state_1'][raw_gene[0]])
    b = log(params['transitions']['begin'][1]) + log(params['emissions']['state_2'][raw_gene[0]])
    hmm_model = [(a, b)]

    for i in range(1, len(raw_gene)):
        prev_a, prev_b = hmm_model[i-1]

        a_stay = prev_a + log(params['transitions']['state_1'][0]) + log(params['emissions']['state_1'][raw_gene[i]])
        a_trans = prev_b + log(params['transitions']['state_2'][0]) + log(params['emissions']['state_1'][raw_gene[i]])
        a = max(a_stay, a_trans)

        b_stay = prev_b + log(params['transitions']['state_2'][1]) + log(params['emissions']['state_2'][raw_gene[i]])
        b_trans = prev_a + log(params['transitions']['state_1'][1]) + log(params['emissions']['state_2'][raw_gene[i]])
        b = max(b_stay, b_trans)

        hmm_model.append((a, b))

    return hmm_model

def run_em(hmm, gene):
    s1_stay_count = 0
    s1_trans_count = 0
    s2_stay_count = 0
    s2_trans_count = 0
    s1_emission_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
    s2_emission_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
    prev_not_cpg = hmm[0][0] > hmm[0][1]

    for i in range(1, len(hmm)):
        a, b = hmm[i]

        # not cpg island
        if a >= b:
            if prev_not_cpg:
                s1_stay_count += 1
            else:
                s2_trans_count += 1
                prev_not_cpg = True
            s1_emission_counts[gene[i]] += 1
        # cpg island
        else:
            if prev_not_cpg:
                s1_trans_count += 1
                prev_not_cpg = False
            else:
                s2_stay_count += 1
            s2_emission_counts[gene[i]] += 1

    s1_sum_emissions = 0
    for key in s1_emission_counts:
        s1_sum_emissions += s1_emission_counts[key]

    s2_sum_emissions = 0
    for key in s2_emission_counts:
        s2_sum_emissions += s2_emission_counts[key]

    params = {
        'emissions': {
            'state_1': {
                'A': s1_emission_counts['A'] / s1_sum_emissions,
                'G': s1_emission_counts['G'] / s1_sum_emissions,
                'C': s1_emission_counts['C'] / s1_sum_emissions,
                'T': s1_emission_counts['T'] / s1_sum_emissions
            },

            'state_2': {
                'A': s2_emission_counts['A'] / s2_sum_emissions,
                'G': s2_emission_counts['G'] / s2_sum_emissions,
                'C': s2_emission_counts['C'] / s2_sum_emissions,
                'T': s2_emission_counts['T'] / s2_sum_emissions
            }
        },

        'transitions': {
            'begin': [0.9999, 0.0001],
            'state_1': [s1_stay_count / (s1_stay_count + s1_trans_count), s1_trans_count / (s1_stay_count + s1_trans_count)],
            'state_2': [s2_trans_count / (s2_stay_count + s2_trans_count), s2_stay_count / (s2_stay_count + s2_trans_count)]
        }
    }
    return params


# state 1 represents normal nucleotide distributions
# state 2 represents cpg islands
def find_cpg_islands(hmm, k):
    in_cpg = False
    ret = []
    i = 0
    start = None
    while k != 0 and i < len(hmm):
        a, b = hmm[i]
        # not in a cpg island
        if a >= b:
            # exiting cpg island
            if in_cpg:
                in_cpg = False
                ret.append((start + 1, i))
                k -= 1

        # in cpg island
        else:
            # entering cpg island
            if not in_cpg:
                start = i
                in_cpg = True
        i += 1
    return ret

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Needs a params filename as an argument")
        exit()

    # parse parameters for viterbi
    with open(sys.argv[1], 'r') as params_file:
        params = json.loads(params_file.read())

    raw_gene = parse_fna()

    for i in range(9):
        hmm = run_viterbi(raw_gene, params)

        print('Transition parameters:')
        print('State 1: %f %f' % (params['transitions']['state_1'][0], params['transitions']['state_1'][1]))
        print('State 2: %f %f' % (params['transitions']['state_2'][0], params['transitions']['state_2'][1]))

        print('Emission parameters:')
        print('CPG Island A\tG\tC\tT')
        A = params['emissions']['state_1']['A']
        G = params['emissions']['state_1']['G']
        C = params['emissions']['state_1']['C']
        T = params['emissions']['state_1']['T']
        print('State_1:   %f\t%f\t%f\t%f' % (A, G, C, T))
        A = params['emissions']['state_2']['A']
        G = params['emissions']['state_2']['G']
        C = params['emissions']['state_2']['C']
        T = params['emissions']['state_2']['T']
        print('State_2:   %f\t%f\t%f\t%f' % (A, G, C, T))
        print('Log probability of most likely path: %f' % max(hmm[-1][1], hmm[-1][0]))
        cpg_islands = find_cpg_islands(hmm, 5)
        print('Number of hits: %d' % len(cpg_islands))
        for start, end in cpg_islands:
            print('\t%d, %d' % (start, end))
        params = run_em(hmm, raw_gene)



