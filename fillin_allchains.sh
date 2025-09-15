#!/bin/bash

# Fill in chain IDs
while getopts "i:o:h" opt; do 
  case $opt in
  i) input="$OPTARG";; # input file
  o) output="$OPTARG";; # output file
  *) echo "Usage: $0 [-i input file][-o output file]"
    exit 1 ;;
  esac
done


awk '{
  if ($1=="ATOM" || $1=="HETATM") {
    atomnum = substr($0,7,5)
    atomname = substr($0,13,4)
    resname = substr($0,18,3)
    chain   = " "
    resnum  = substr($0,23,4)

    # assign chains based on atom number (adjust cutoffs as needed)
    if (atomnum+0 <= 11750) {
      chain="A"
    } else if ((11751 <= atomnum+0 ) && (atomnum+0 <= 23500)) {
      chain="B"
    } else if ((23501 <= atomnum+0) && (atomnum+0 <= 35250)) {
      chain="C"
    } else {
      chain="D"
    }

    # rebuild line with chain correctly in col 22
    printf("%-6s%5s %-4s %3s %1s%4s%s\n",
           substr($0,1,6), atomnum, atomname, resname, chain, resnum, substr($0,27))
  } else {
    print $0
  }
}' "$input" > "$output"

# Fill in the atom column 
obabel "$output" -O "$output"
# obabel medoid_cluster0_fp11.pdb -O medoid_cluster0_fp11.pdb