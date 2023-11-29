echo -n "Nodes alive: "
qhost | grep -v -e ' - ' -e '-----' -e 'HOSTNAME' | wc -l
qhost | grep -v -e ' - ' -e '-----' -e 'HOSTNAME' | awk 'BEGIN{X=0}{X+=strtonum($3)}END{print "CPUs: " X}'
qhost | grep -v -e ' - ' -e '-----' -e 'HOSTNAME' | awk 'BEGIN{X=0}{X+=strtonum($4)}END{print "Total Load: " X}'
echo "Nodes dead: "
qhost | grep ' - ' | grep -v  -e '-----' -e 'HOSTNAME' -e global | awk '{print $1}'
echo "Running per user: "
#qstat | awk '{print $4}' | sort | grep -v -e '^$' -e user | uniq -c | sort -n
qstat -u '*' -g d | awk '{print $4,$5}' | sort | grep -v -e '^ $' -e user | grep 'r$' | uniq -c | sort -nr
echo "Queued per user: "
qstat -u '*' -g d | awk '{print $4,$5}' | sort | grep -v -e '^ $' -e user | grep 'qw$' | uniq -c | sort -nr
echo "Other per user: "
qstat -u '*' -g d | awk '{print $4,$5}' | sort | grep -v -e '^ $' -e user | grep -v -e 'qw$' -e 'r$' | uniq -c | sort -nr