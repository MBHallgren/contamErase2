import os
import sys

path = '/home/people/plan/kma_nanomlst/seqs/run38/guppy_5.0.11_sup'
files = os.listdir(path)

for item in files:
    name = item.split('.')[0]
    cmd = 'kma -i {} -o {} -t_db /home/people/malhal/papers/rmlst/test/38_db/specie_db -ID 10 -ont -md 10 -matrix -eq 14 -mct 0.5'.format(os.path.join(path, item), name)
    os.system(cmd)