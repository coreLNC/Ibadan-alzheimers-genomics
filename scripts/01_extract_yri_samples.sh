H="ftp.1000genomes.ebi.ac.uk"
P="vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
curl -s "http://${H}/${P}" | awk '$2=="YRI" {print $1}' > yri_samples.txt && wc -l yri_samples.txt
