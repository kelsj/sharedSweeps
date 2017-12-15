#run script to add local rho to unstd iHS file 
perl makeInputFor.iHS.std.by.freq.rho.bins.pl FIN 22

#run script to calculate stdized iHS using stdization bins file and unstd iHS
perl iHS.std.by.freq.rho.bins.pl FIN 22  