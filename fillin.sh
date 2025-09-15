!/bin/bash

input="added_names.pdb"
output="md_out.pdb"

awk '{
  if ($1=="ATOM" || $1=="HETATM") {
    atomnum = substr($0,7,5)
    atomname = substr($0,13,4)
    resname = substr($0,18,3)
    chain   = " "
    resnum  = substr($0,23,4)

    # assign chain A if atom number <= 2699
    if (atomnum+0 <= 2699) {
      chain="A"
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


# #### By TER records ####
# #!/bin/bash

# input="md.pdb"
# output="md_out.pdb"

# chain1="A"
# chain2="D"

# awk -v chain="$chain" '{
#   if ($1=="ATOM" || $1=="HETATM") {
#     atomnum = substr($0,7,5)
#     atomname = substr($0,13,4)
#     resname = substr($0,18,3)
#     resnum  = substr($0,23,4)

#     printf("%-6s%5s %-4s %3s %1s%4s%s\n",
#            substr($0,1,6), atomnum, atomname, resname, chain, resnum, substr($0,27))
#   } else if ($1=="TER") {
#     print $0
#     chain=chain2   # switch after TER
#   } else {
#     print $0
#   }
# }' "$input" > "$output"


