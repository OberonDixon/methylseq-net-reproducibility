#!/bin/bash
set -euo pipefail
HG38_2BIT=$(methylseqnet-repro-path genomes)/hg38.2bit
CPG_ANNOTATION=$(methylseqnet-repro-path annotations)/hg38_cpgIsland.bed
twoBitToFa "$HG38_2BIT" stdout | maskOutFa stdin hard stdout \
  | /clusterfs/nilah/oberon/repos/kent/src/utils/cpgIslandExt/cpg_lh /dev/stdin 2> cpg_lh.err \
    |  awk '{$2 = $2 - 1; width = $3 - $2;  printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n", $1, $2, $3, $5, $6, width, $6, width*$7*0.01, 100.0*2*$6/width, $7, $9);}' \
     | sort -k1,1 -k2,2n > "$CPG_ANNOTATION"