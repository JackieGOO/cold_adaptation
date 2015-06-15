from subprocess import check_output
import os

# set the environmental variable of IUPred
os.environ['IUPred_PATH'] = "/home/yzolotarov/iupred"

res = check_output(["/home/yzolotarov/iupred/iupred", 
	"/home/yzolotarov/iupred/P53_HUMAN.seq", "long"])

res = res.split('\n')
scores = [float(line.split()[2]) for line in res[:-1] if '#' not in line]

disordered = [s for s in scores if s >= 0.5]

print disordered, len(disordered), len(scores)
