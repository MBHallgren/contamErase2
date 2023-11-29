import os
import sys

def eval_rmlst(validates_genes, unique_genes):
    conditions = {
        'scenario_1': {
            'total_varified_genes': 53,
            'total_unique_genes': 1
        },
        'scenario_2': {
            'total_varified_genes': 52,
            'total_unique_genes': 3
        },
        'scenario_3': {
            'total_varified_genes': 50,
            'total_unique_genes': 5
        },
        'scenario_4': {
            'total_varified_genes': 40,
            'total_unique_genes': 10
        }

    }

    for scenario in conditions:
        if validates_genes >= conditions[scenario]['total_varified_genes'] and unique_genes >= conditions[scenario]['total_unique_genes']:
            return True
    return False

path = '/home/people/malhal/contamErase/data/nanopore/intra/paper/results/'
path_list = os.listdir(path)
results_dict = {}

for item in path_list:
    print (item)
    if not item.endswith('.py') and not item.endswith('one'):
        new_path = path + item
        results_dict[item] = {}
        if os.path.exists(new_path + '/rmlst_staphylococcus_aureus/rmlst_results.txt'):
            with open(new_path + '/rmlst_staphylococcus_aureus/rmlst_results.txt', 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#N'):
                        name = line.split(' ')[0]
                        results_dict[item][name] = {}
                    if line.startswith('#Validated genes:'):
                        results_dict[item][name]['validated_genes'] = line.split(' ')[-1]
                    if line.startswith('#Validated unique genes:'):
                        results_dict[item][name]['unique_genes'] = line.split(' ')[-1]
                    if line.startswith('#Unvalidated genes:'):
                        results_dict[item][name]['unvalidated_genes'] = line.split(' ')[-1]
                        results_dict[item][name]['color'] = 'red'
                if len(results_dict[item]) == 0:
                    results_dict[item] = 'NA'
                elif len(results_dict[item]) > 1:
                    found_rmlst = False
                    top_unique_genes = 0
                    best_template = ''
                    best_hit = {}
                    for name in results_dict[item]:
                        if eval_rmlst(int(results_dict[item][name]['validated_genes']), int(results_dict[item][name]['unique_genes'])):
                            best_hit['color'] = 'green'
                            found_rmlst = True
                            if int(results_dict[item][name]['unique_genes']) > top_unique_genes:
                                top_unique_genes = int(results_dict[item][name]['unique_genes'])
                                best_hit = results_dict[item][name]
                                best_hit['color'] = 'green'
                                best_template = name

                        else:
                            if not found_rmlst:
                                if int(results_dict[item][name]['unique_genes']) > top_unique_genes:
                                    top_unique_genes = int(results_dict[item][name]['unique_genes'])
                                    best_hit = results_dict[item][name]
                                    best_template = name
                    results_dict[item] = {}
                    results_dict[item][best_template] = best_hit
        else:
            results_dict[item] = 'NA'

print (results_dict)

outfile = 'csv/nanopore.csv'.format(item)
with open (outfile, 'w') as f:
    print ('id', 'name', 'validated_genes', 'unique_genes', 'unvalidated_genes', 'color', sep=',', file=f)
    for id in path_list:
        if id in results_dict:
            if results_dict[id] == 'NA':
                print (id, 'NA', '0', '0', '0', 'black', sep=',', file=f)
            else:
                for name in results_dict[id]:
                    print (id, name, results_dict[id][name]['validated_genes'], results_dict[id][name]['unique_genes'], results_dict[id][name]['unvalidated_genes'], results_dict[id][name]['color'], sep=',', file=f)