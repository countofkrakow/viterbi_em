from data import parse_fna
from math import log, exp
import sys
import json

def run_viterbi(gene, params, k):
    raw_gene = parse_fna()

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

    

    return

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Needs a params filename as an argument")
        exit()

    # parse parameters for viterbi
    with open(sys.argv[1], 'r') as params_file:
        params = json.loads(params_file.read())


