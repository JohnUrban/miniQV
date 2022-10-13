grep ^total */*/*/minimap2-endbonus.cigar-QV-table.txt | awk 'OFS="\t" {gsub("/","\t"); print $2,$3,$5,$8}'
