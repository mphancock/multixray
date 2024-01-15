# Check usage of /wynton/gruop/sali.
beegfs-quota -h -p group -u mhancock

# Check number of files/storage on /wynton/group.
beegfs-ctl --getquota --storagepoolid=11 --uid "mhancock"
