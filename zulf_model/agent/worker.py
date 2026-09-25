"""Child-process entry point for background jobs."""
import sys

from .jobs import run_worker

if __name__ == "__main__":
    run_worker(sys.argv[1])
