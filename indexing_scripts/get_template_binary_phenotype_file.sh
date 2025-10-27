# This takes as input a graph and writes a template binary phenotype file to stdout
# The template will have all sample names with phenotype 1

GRAPH=$1

printf "FID\tIID\tPHENO\n" 
vg paths -M -x $GRAPH | tail -n +2 | awk '{ if ( $3 == "NO_SAMPLE_NAME") {print $5} else {print $3}}' | sort | uniq | sed -r 's/(.*)/\1\t\1\t1/g'

