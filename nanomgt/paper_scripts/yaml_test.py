import yaml

with open('../intra_species_conditions.yaml', 'r') as yaml_file:
    conditions = yaml.safe_load(yaml_file)
print (conditions)