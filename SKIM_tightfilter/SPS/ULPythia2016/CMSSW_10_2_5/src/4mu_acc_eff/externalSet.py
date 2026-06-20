import subprocess

externalList = []
prefix = "/eos/user/l/leyao/26JJ/MC_Maker/HelacOnia2016/CMSSW_10_6_20/src/NTUPLE/ggpsi1psi1g_gpt0p8/"
res = subprocess.check_output(["ls", prefix])
folders = res.decode('utf-8').strip().split()
for folder in folders:
    res = subprocess.check_output(["ls", prefix + folder])
    files = res.decode('utf-8').strip().split()
    for file in files:
        externalList.append(prefix + folder + '/' + file)
# print("Find", len(externalList), "external datasets.\n")