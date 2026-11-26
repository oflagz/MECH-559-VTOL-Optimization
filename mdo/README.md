From this folder:

    python3 system.py

That runs every pack and blade count, 540 cases. It prints each one as it finishes and writes search_log.txt here. The last lines are the best design, or "no feasible design".

One pack, without the full search:

    python3 -c "import system; print(system.solve_one(10, 3, 3))"

Needs numpy, scipy, and matplotlib.
