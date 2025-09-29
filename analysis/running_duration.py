import os
import datetime
import struct
import matplotlib.pyplot as plt


def get_runtime(parent_dir = "."):
    results = {}

    for folder in sorted(os.listdir(parent_dir)):
        if folder.startswith("run_") and os.path.isdir(os.path.join(parent_dir, folder)):
            run_path = os.path.join(parent_dir, folder)
            file1 = os.path.join(run_path, "a_e_t0.csv")
            file2 = os.path.join(run_path, "particles.txt")

            if os.path.exists(file1) and os.path.exists(file2):
                # get modification times
                t1 = os.path.getmtime(file1)
                t2 = os.path.getmtime(file2)

                # convert to datetime
                dt1 = datetime.datetime.fromtimestamp(t1)
                dt2 = datetime.datetime.fromtimestamp(t2)

                # runtime difference
                duration = (dt2 - dt1).total_seconds() / 60.0

                # Read particle file
                with open(file2, 'rb') as pfile:
                    try:
                        Np_t0   = struct.unpack('i', pfile.read(4))[0]  # Read int
                        time_t0 = struct.unpack('d', pfile.read(8))[0]  # Read double
                        r_dust  = struct.unpack('d', pfile.read(8))[0]  # Read double
                    except struct.error:
                        continue

                results[folder] = {
                    "start": dt1,
                    "end": dt2,
                    "duration": duration,
                    "N_particle": Np_t0,
                    "r_dust": r_dust
                }
    
    # --- print results ---
    for folder, info in results.items():
        print(f"[{folder}]")
        print(f"  Start     : {info['start']}")
        print(f"  End       : {info['end']}")
        print(f"  Duration  : {info['duration']:.2f} minutes")
        print(f"  N_particle: {info['N_particle']}")
        print(f"  r_dust    : {info['r_dust']:.3e}\n")
                            
    return results

def plot_runtime_vs_rdust(results, output_dir = "."):
    # extract values
    runs = list(results.keys())
    durations = [results[r]["duration"] for r in runs]
    r_dusts = [results[r]["r_dust"] for r in runs]

    # make scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(r_dusts, durations, c="blue", marker="o")

    # axes
    plt.xscale("log")
    plt.xlabel("Dust Radius (m)")
    plt.ylabel("Runtime (minutes)")
    plt.title("Simulation Duration vs dust radius")
    plt.grid(True)

    # annotate points with run names
    for run, x, y in zip(runs, r_dusts, durations):
        plt.annotate(run, (x, y), textcoords="offset points", xytext=(5,5), fontsize=6)

    plt.tight_layout()
   
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "runtime_vs_rdust.png"), dpi=150)
    
    plt.show()


if __name__ == '__main__':
    # parent directory containing run_001, run_002, ...
    parent_dir = "."

    output_dir = "."
    
    run_time = get_runtime(parent_dir)
    #plot_runtime_vs_rdust(run_time, output_dir)
