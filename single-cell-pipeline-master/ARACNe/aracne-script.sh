#!bin/bash
ARACNE_BIN=/ifs/scratch/c2b2/ac_lab/CZI/aracne-1.11/aracne-ap.jar
ARACNE_DIR='/ifs/scratch/c2b2/ac_lab/lv2395/ARACNe'
BOOTSTRAP_NUM=200
P_THRESH=1E-8

## Process Arguments
RUN_ID=$1 # a name for this ARACNe run
BASE_DIR=$2 # the base directory where results will be generated / stored
EXP_FILE=$3 # the expression file in .rds format (features X samples)

## Referenced Variables
A_TABLE=${BASE_DIR}${RUN_ID}_ARACNe-table.tsv

## Prep data from .rds to ARACNe Table
qsub -N adp_${RUN_ID} -wd ${BASE_DIR} -l mem=16G -l s_rt=:30: -l h_rt=:30: -b y \
/nfs/apps/R/3.5.1/bin/Rscript ${ARACNE_DIR}/aracne_data-prep.r --rds=${EXP_FILE} --out=${A_TABLE}

## Create sub-directory for each regulator set
REG_SETS=()
for (( c = 4; c <= $#; c++ ))
do
	FILE=${!c}
	REG=${FILE##*/}
	REG=${REG%%.*}
	mkdir -p ${BASE_DIR}/${RUN_ID}_${REG}
	REG_SETS=(${REG_SETS[@]} ${REG})
done

## Calculate MI, generate boostraps, and consolidate for each reg set
HOLD_LIST=''
for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
do
	## Set iteration variables
	d=$((c + 4))
	REG_FILE=${!d}
	REG=${REG_SETS[c]}
	WORK_DIR=${BASE_DIR}${RUN_ID}_${REG}
	RUN_NAME=${RUN_ID}_${REG}
	HOLD_LIST="$HOLD_LIST,aCon_${RUN_NAME}"
	## Threshold calculation
	qsub -N aThresh_${RUN_NAME} -hold_jid adp_${RUN_ID} -wd ${WORK_DIR} -l mem=16G -l s_rt=2:0:0 -l h_rt=2:0:0 -b y \
	java -Xmx5G -jar ${ARACNE_BIN} -e ${A_TABLE} -o ${WORK_DIR} --tfs ${REG_FILE} --pvalue ${P_THRESH} --seed 666 --calculateThreshold
	## Bootstrap generation
	for (( i = 1; i <= $BOOTSTRAP_NUM; i++ ))
	do
		qsub -N abs_${RUN_NAME} -hold_jid aThresh_${RUN_NAME} -wd ${WORK_DIR} -l mem=32G -l s_rt=8:0:0 -l h_rt=8:0:0 -b y \
		java -Xmx5G -jar ${ARACNE_BIN} -e ${A_TABLE} -o ${WORK_DIR} --tfs ${REG_FILE} --pvalue ${P_THRESH} --seed $i
	done
	## Consolidation
	qsub -N aCon_${RUN_NAME} -hold_jid abs_${RUN_NAME} -wd ${WORK_DIR} -l mem=16G -l s_rt=16:0:0 -l h_rt=16:0:0 -b y \
	/nfs/apps/R/3.5.1/bin/Rscript ${ARACNE_DIR}/aracne_consolidate.r ${WORK_DIR} ${A_TABLE} ${REG_FILE} bonferroni 0.01
done 

## Merge and reg process jobs
CAT_LIST=''
for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
do
	REG=${REG_SETS[c]}
	WORK_DIR=${BASE_DIR}${RUN_ID}_${REG}
	CAT_LIST="$CAT_LIST ${WORK_DIR}/finalNetwork_4col.tsv"
done
HOLD_LIST=${HOLD_LIST#?}
qsub -N aMerge_${RUN_ID} -hold_jid ${HOLD_LIST} -wd ${BASE_DIR} -l mem=16G -l s_rt=16:0:0 -l h_rt=16:0:0 -b y \
bash ${ARACNE_DIR}/aracne_merge.sh ${CAT_LIST} ${BASE_DIR}${RUN_ID}_finalNet-merged.tsv
qsub -N arp_${RUN_ID} -hold_jid aMerge_${RUN_ID} -wd ${BASE_DIR} -l mem=16G -l s_rt=16:0:0 -l h_rt=16:0:0 -b y /nfs/apps/R/3.5.1/bin/Rscript \
${ARACNE_DIR}/aracne_regProc.r --a_file=${BASE_DIR}${RUN_ID}_finalNet-merged.tsv --exp_file=${A_TABLE} --out_dir=${BASE_DIR} --out_name=${RUN_ID}

## Wait til all consolidation is complete
# FIN_CHECK=0
# while [ ${FIN_CHECK} = 0 ]
# do
# 	FIN_CHECK=1
# 	for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
# 	do
# 		echo 'Generating networks...'
# 		sleep 1m
# 		REG=${REG_SETS[c]}
# 		WORK_DIR=${BASE_DIR}${RUN_ID}_${REG}
# 		if [[ -f ${WORK_DIR}/finalNetwork_4col.tsv ]] ; then
# 			FIN_CHECK=$((FIN_CHECK*1))
# 		else
# 			FIN_CHECK=$((FIN_CHECK*0))
# 		fi
# 	done
# done

# ## Merge and reg process
# CAT_LIST=''
# for (( c = 0; c < ${#REG_SETS[@]}; c++ ))
# do
# 	REG=${REG_SETS[c]}
# 	WORK_DIR=${BASE_DIR}${RUN_ID}_${REG}
# 	CAT_LIST="$CAT_LIST ${WORK_DIR}/finalNetwork_4col.tsv"
# done
# cat ${CAT_LIST} > ${BASE_DIR}${RUN_ID}_finalNet-merged.tsv
# echo ${A_TABLE}
# qsub -N regProc_${RUN_ID} -wd ${BASE_DIR} -l mem=16G -l s_rt=16:0:0 -l h_rt=16:0:0 -b y /nfs/apps/R/3.5.1/bin/Rscript \
# ${ARACNE_DIR}/aracne_regProc.r --a_file=${BASE_DIR}${RUN_ID}_finalNet-merged.tsv --exp_file=${A_TABLE} --out_dir=${BASE_DIR} --out_name=${RUN_ID}
