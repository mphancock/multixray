# Check usage of /wynton/group/sali.
beegfs-quota -h -p group -u mhancock

# Check number of files/storage on /wynton/home.
beegfs-ctl --getquota --storagepoolid=11 --uid "mhancock"
