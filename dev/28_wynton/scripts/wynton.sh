# Check usage of /wynton/group/sali.
beegfs-quota -h -p group -u mhancock

# Check number of files/storage on /wynton/home.
beegfs-ctl --getquota --storagepoolid=11 --uid "mhancock"

## monitor i/o on beegfs
beegfs-ctl --userstats --nodetype=storage --names  --perinterval=10 --interval=10 --filter=<YOUR_USERNAME>