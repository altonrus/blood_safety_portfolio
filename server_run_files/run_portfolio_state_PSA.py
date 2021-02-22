import os, sys, time

n_iters = 10000
n_jobs = 40
iters_per_job = n_iters/n_jobs

iter_start = 1

for job in range(n_jobs):
  with open("PSA_portfolio_%03d.sbatch" % (job), 'w') as f:
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("#SBATCH --job-name=PSA_portfolio_%03d\n" % (job))
    f.write("#SBATCH --time=5:00:00\n")
    f.write("#SBATCH --ntasks=1\n")
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --mem-per-cpu=1G\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --ntasks-per-node=1\n")
    #f.write("#SBATCH --mail-type=END\n")
    #f.write("#SBATCH --mail-user=altonr@stanford.edu\n")
    f.write("\n")
    f.write("module load R\n")
    f.write("Rscript PSA_for_sherlock.R %d %d" % (iter_start, iter_start+iters_per_job -1))
  
  os.system("sbatch PSA_portfolio_%03d.sbatch" % (job))
  iter_start += iters_per_job
