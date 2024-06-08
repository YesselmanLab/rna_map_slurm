import subprocess
import json
from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class SlurmOptions:
    name: str
    time: str = "12:00:00"
    mem_per_job: str = "2GB"
    cpus_per_task: int = 1
    extra_header_cmds: str = ""


def generate_submit_file(path, jobs):
    f = open(path, "w")
    for j in jobs:
        f.write(f"sbatch {j}\n")
    f.close()


def get_job_header(args: SlurmOptions, job_dir=None):
    if job_dir is not None:
        output_path = f"{job_dir}/{args.name}.out"
        error_path = f"{job_dir}/{args.name}.err"
    else:
        output_path = f"{args.name}.out"
        error_path = f"{args.name}.err"

    return f"""#!/bin/bash
#SBATCH --time={args.time}         
#SBATCH --mem-per-cpu={args.mem_per_job}       
#SBATCH --cpus-per-task={args.cpus_per_task}    
#SBATCH --output={output_path}
#SBATCH --error={error_path}

{args.extra_header_cmds}
"""


def get_user_jobs(user):
    try:
        # Run the squeue command to get jobs of the user in JSON format
        result = subprocess.run(
            [
                "squeue",
                "--user",
                user,
                "--format",
                "%i %P %j %u %t %M %l %D %R",
                "--noheader",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

        # Split the output into lines
        lines = result.stdout.strip().split("\n")

        # Define keys for the job attributes
        keys = [
            "JobID",
            "Partition",
            "Name",
            "User",
            "State",
            "Time",
            "TimeLimit",
            "Nodes",
            "NodeList",
        ]

        # Parse the lines into a list of dictionaries
        jobs = [dict(zip(keys, line.split())) for line in lines]

        return jobs
        # Write the jobs data to a JSON file
        #with open(f"{user}_jobs.json", "w") as f:
        #    json.dump(jobs, f, indent=4)

        #print(f"Job data for user '{user}' has been saved to {user}_jobs.json")

    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def search_jobs(user, search_term):
    try:
        # Read the jobs data from the JSON file
        with open(f"{user}_jobs.json", "r") as f:
            jobs = json.load(f)

        # Filter jobs based on the search term
        filtered_jobs = [
            job for job in jobs if search_term.lower() in json.dumps(job).lower()
        ]

        # Print the filtered jobs
        if filtered_jobs:
            print(f"Found {len(filtered_jobs)} job(s) matching '{search_term}':")
            for job in filtered_jobs:
                print(json.dumps(job, indent=4))
        else:
            print(f"No jobs found matching '{search_term}'")

    except FileNotFoundError:
        print(
            f"No job data found for user '{user}'. Please run the get_user_jobs function first."
        )
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    user = input("Enter the SLURM username: ")
    get_user_jobs(user)

    search_term = input("Enter a search term to filter jobs: ")
    search_jobs(user, search_term)
