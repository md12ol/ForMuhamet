import subprocess
import itertools
import time

# params list
mut_rates = [0.027, 0.05]
num_islands = [1, 4]

# path to cpp
executable_name = "./cmake-build-release-wsl/pagmo_nsga2"

combinations = list(itertools.product(num_islands, mut_rates))

total_runs = len(combinations)

# starting simu
print(f"=== STARTING AUTOMATED HYPERPARAMETER SWEEP ===")
print(f"Total planned simulations: {total_runs}")
print("-" * 50)

for i, (isl, mut) in enumerate(combinations, 1):
    print(f"[{i}/{total_runs}] Setup: Islands={isl} | Mut={mut}")

    command = [executable_name, str(isl), str(mut)]

    start_time = time.time()

    try:
        subprocess.run(command, check=True)

        duration = time.time() - start_time
        print(f"   -> Done! Duration: {duration/60:.2f} minutes.\n")

    except subprocess.CalledProcessError as e:
        print(f"   -> ERROR: Simulation {i} crashed. Skipping to next...\n")
    except FileNotFoundError:
        print(f"   -> FATAL ERROR: Executable '{executable_name}' not found! Check your path.")
        break

print("=== ALL COMPUTATIONS HAVE BEEN SUCCESSFULLY COMPLETED! ===")