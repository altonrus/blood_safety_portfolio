import os, sys, time

n_iters = 10120
n_jobs = 220
iters_per_job = n_iters/n_jobs
N_rep = 4000000

iter_start = 1

for job in range(n_jobs):
  with open("PSA_ZIKV_%03d.sbatch" % (job), 'w') as f:
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("#SBATCH --job-name=PSA_ZIKV_%03d\n" % (job))
    f.write("#SBATCH --time=47:00:00\n")
    f.write("#SBATCH --ntasks=1\n")
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --mem-per-cpu=1G\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --ntasks-per-node=1\n")
    #f.write("#SBATCH --mail-type=END\n")
    #f.write("#SBATCH --mail-user=altonr@stanford.edu\n")
    f.write("\n")
    f.write("module load python/3.6.1\n")
    f.write("module load py-numpy/1.14.3_py36\n")
    f.write("module load py-pandas/0.23.0_py36\n")
    f.write("python3 ZIKV_sim_sherlockPSA.py %d %d %d PSA_ZIKV_%03d.csv" % (iter_start, iters_per_job, N_rep, job))
  
  os.system("sbatch PSA_ZIKV_%03d.sbatch" % (job))
  time.sleep(1)
  iter_start += iters_per_job
