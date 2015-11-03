#!/bin/bash


# if you want to use S step you must set this
REAL_TE_INS=$1

INS=$2



####################################################################################################
########################################## evaluate simulation #####################################
####################################################################################################

echo $INS

mkdir evaluate_INS_SIM

INS_NAME=$(basename ${INS/.gff3/})

#TP
intersectBed -a $INS -b $REAL_TE_INS -u > evaluate_INS_SIM/$INS_NAME.true_positives.gff3

#FP
intersectBed -a $INS -b $REAL_TE_INS -v  > evaluate_INS_SIM/$INS_NAME.false_positives.gff3

#detected
intersectBed -a $REAL_TE_INS -b $INS -u > evaluate_INS_SIM/$INS_NAME.detected.gff3

#undetected
intersectBed -a $REAL_TE_INS -b $INS  -v  > evaluate_INS_SIM/$INS_NAME.undetected.gff3

wc -l evaluate_INS_SIM/*

/users/so/phdcourse/jitterbug/jitterbug-code/scripts/plot_gff_annots_2files.py -g evaluate_INS_SIM/$INS_NAME.true_positives.gff3,TP -G evaluate_INS_SIM/$INS_NAME.false_positives.gff3,FP -o evaluate_INS_SIM/$INS_NAME.TP-FP.pdf


