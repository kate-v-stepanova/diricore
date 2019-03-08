#!/bin/bash


source ./diricore_virtualenv/bin/activate

set -e;
set -u;


OUTDIR="./data/output/rpf_5p_density/";
INDIR="./data/output/tophat_out/";
PLOTDIR="./data/output/figures/";
SAMPLENAME_FILE="./data/input/metadata/rpf_density_samplenames.tsv";
CONTRAST_FILE="./data/input/metadata/rpf_density_contrasts.tsv";

species=$1;
projectname=$2;
minreads=$3;

INDEXDATA_FILE="./staticdata/${species}/transcript_data.hdf5";
MAPS_FILE="./staticdata/${species}/codon_regions.width_61.hdf5";
MAPSSTART_FILE="./staticdata/${species}/codon_regions.START_Other_ATG.width_61.hdf5";


###
# { setup
of="${OUTDIR}/${projectname}.txcoord_counts.hdf5";

mkdir ${OUTDIR} || true;
mkdir -p "${PLOTDIR}/rpf_5p_density_plots/" || true;

# map RPFs to transcriptome coordinates
./diricore/bin/map_rpfs_to_transcriptome_positions.py \
    -t "${INDEXDATA_FILE}" \
    -o "${of}" \
    -b <(\
        ls -1 ${INDIR}/*/accepted_hits.hqmapped_dedup.bam | sort -V | while read fn; do
                b=$(basename $(dirname "$fn"));
                b=${b%%.*};

                echo -e "${b}\t${fn}";
            done \
        )
# }

# { generate RPF density shift plots
datafile="${of}";
echo -ne "\
Ala\tGCA,GCC,GCG,GCT
Arg\tCGA,CGC,CGG,CGT,AGA,AGG
Asn\tAAC,AAT
Asp\tGAC,GAT
Cys\tTGC,TGT
Gln\tCAA,CAG
Glu\tGAA,GAG
Gly\tGGA,GGC,GGG,GGT
His\tCAC,CAT
Ile\tATA,ATC,ATT
Leu\tCTA,CTC,CTG,CTT,TTA,TTG
Lys\tAAA,AAG
Met\tATG
Phe\tTTC,TTT
Pro\tCCA,CCC,CCG,CCT
Ser\tTCA,TCC,TCG,TCT,AGC,AGT
Thr\tACA,ACC,ACG,ACT
Trp\tTGG
Tyr\tTAC,TAT
Val\tGTA,GTC,GTG,GTT
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/${projectname}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python ./diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${datafile}" \
        ${MAPS_FILE} \
        ${codons}
done

# RPF density at special codons
# (START/other ATG)
echo -ne "\
ATG_split\tSTART_ATG,Other_ATG
" | while read aa codongroupstr; do
    of="${PLOTDIR}/rpf_5p_density_plots/${projectname}.m${minreads}.${aa}.rpf_5p_density_shift_plot.pdf";
    codons=$(echo $codongroupstr | sed 's/,/ /g');

    python ./diricore/bin/plot_rpf_5p_density.py \
        -c "${CONTRAST_FILE}" \
        -n "${SAMPLENAME_FILE}" \
        -o "${of}" \
        -m ${minreads} \
        "${datafile}" \
        ${MAPSSTART_FILE} \
        ${codons}
done
