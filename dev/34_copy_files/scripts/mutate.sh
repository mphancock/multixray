#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=01:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N mutate
#$ -pe smp 8
#$ -t 1-1
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_219_cctbx

# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "131_1_state/118669" --new "165_J1_i/0"
# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "129_2_state/118663" --new "165_J1_i/1"
# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "137_4_state_1_cond/131215" --new "165_J1_i/2"
# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "153_N8_J1/558377" --new "165_J1_i/3"
# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "146_1_state_2_cond/515356" --new "164_J2_ik/0"
# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "130_2_cond/118668" --new "164_J2_ik/1"
# python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "138_4_state_2_cond/131217" --new "164_J2_ik/2"
python ~/xray/dev/34_copy_files/scripts/mutate.py --orig "154_N8_J2/558391" --new "164_J2_ik/3"
