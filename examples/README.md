# Example usage

This directory contains python scripts using Conquest2a to do various things, including file writing and plotting (p)DOS and band structures with matplotlib. Re-use or modify as you wish.

In particular, as of CONQUEST v1.4, all output files are written to the same directory as whatever corresponding `Conquest_input` is in. All associated files (`.ion`, coordinates etc.) are assumed to also be in the same directory. Therefore most of these scripts start with resolving an absolute path to some directory to take care of issues.